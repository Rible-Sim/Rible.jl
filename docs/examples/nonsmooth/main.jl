using LinearAlgebra
using Statistics
using SparseArrays
using StaticArrays
using Rotations
# using Parameters
import GeometryBasics as GB
using OffsetArrays
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using LaTeXStrings
using RecursiveArrayTools
using Interpolations
using EponymTuples
using CircularArrays
using TypeSortedCollections
using Printf
using CoordinateTransformations
using Meshing
using Test
using Unitful
using Match
using FileIO
using Revise
import TensegrityRobots as TR
import Meshes
cd("examples\\nonsmooth")
include("def.jl"); includet("def.jl")
# includet("plotting.jl")
include("../analysis.jl"); includet("../analysis.jl")
include("../vis.jl"); includet("../vis.jl")
figdir::String = raw"C:\Users\luo22\OneDrive\Papers\Ph.D.Thesis\ns"
# ----------- unibot ---------------
unibot = uni(0.0;
            μ = 0.9,
            e = 0.0,
			z0 = 0.2
)

plot_traj!(unibot;
	# zlims = (0,1),
	showground = true
)

function uni_dynfuncs(bot)
    (;tg) = bot
    function F!(F,q,q̇,t)
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_tensiles!(tg)
        TR.apply_gravity!(tg)
        F .= TR.generate_forces!(tg)
    end
    function Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t)
		∂F∂q̌ .= 0
        ∂F∂q̌̇ .= 0
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_tensiles!(tg)
        TR.build_∂Q̌∂q̌!(∂F∂q̌,tg)
        TR.build_∂Q̌∂q̌̇!(∂F∂q̌̇,tg)
    end
	bars = TR.get_rigidbars(tg)
	bar1 = bars[1]
	rbid = bar1.prop.id

	function prepare_contacts!(contacts,q)
		TR.update_rigids!(tg,q)
		for c in contacts
			gap = bar1.state.rps[c.id][3]
			TR.activate!(c,gap)
		end
		active_contacts = filter(contacts) do c
			c.state.active
		end
		na = length(active_contacts)
		inv_μ_vec = ones(eltype(q),3na)
		n = [0,0,1.0]
		for (i,ac) in enumerate(active_contacts)
			(;id,state) = ac
			state.frame = TR.SpatialFrame(n)
			inv_μ_vec[3(i-1)+1] = 1/ac.μ
		end
        es = [ac.e for ac in active_contacts]
		gaps = [ac.state.gap for ac in active_contacts]
		H = Diagonal(inv_μ_vec)
        active_contacts, gaps, H, es
    end

	function get_directions_and_positions(active_contacts,q)
		na = length(active_contacts)
		TR.update_rigids!(tg,q)
		D = Matrix{eltype(q)}(undef,3na,length(q))
		ŕ = Vector{eltype(q)}(undef,3na)
		for (i,ac) in enumerate(active_contacts)
			(;id,state) = ac
			(;n,t1,t2) = state.frame
            C = bar1.state.cache.Cps[id]
			CT = C*TR.build_T(tg,rbid)
			dm = hcat(n,t1,t2) |> transpose
			D[3(i-1)+1:3(i-1)+3,:] = dm*CT
			ŕ[3(i-1)+1:3(i-1)+3] = dm*bar1.state.rps[id]
		end
		D,ŕ
	end

	@eponymtuple(F!,Jac_F!,prepare_contacts!,get_directions_and_positions)
end

tspan = (0.0,1.5)
h = 1e-3

prob = TR.SimProblem(unibot,uni_dynfuncs)

TR.solve!(prob,TR.Zhong06();tspan,dt=h,ftol=1e-14,maxiters=50,exception=false)

TR.solve!(prob,TR.ZhongCCP();tspan,dt=h,ftol=1e-14,maxiters=50,exception=false,verbose=true)

plot_traj!(unibot)

GM.activate!()
with_theme(theme_pub;
		markersize = 0.4fontsize,
		fontsize = 0.8fontsize,
		figure_padding = (fontsize,0,1.5fontsize,fontsize),
	) do 
	plot_traj!(unibot;
		figsize = (0.9tw,0.7tw),
		xlims = (-0.01,0.15),
		ylims = (-0.15,0.05),
		zlims = (-0.01,0.3),
		gridsize = (2,3),
		attimes = [0,0.179,0.255,0.412,0.694,1.177],
		AxisType=Axis3,
		doslide =false,
		showinfo = false,
		showlabels = false,
        showpoints = false,
		showground = true,
		sup! = (ax,_) -> begin
			# ax.zlabeloffset = 2fontsize
			# ax.xlabeloffset = 1.5fontsizew
			# ax.azimuth = 4.335530633326986
			# ax.elevation = 0.27269908169872403
		end,
		# figname = "unibot",
	)
end

GM.activate!(); plotsave_energy(unibot)
CM.activate!(); plotsave_energy(unibot,"unibot_energy")

unibots = [
	begin
		TR.solve!(
			TR.SimProblem(
				uni(0.0; μ = 0.9, e, z0 = 0.2),
				uni_dynfuncs
			),
			TR.ZhongCCP();
			tspan,dt=h,ftol=1e-14,maxiters=50,exception=false
		)
	end
	for e = [1.0,0.7,0.3,0.0]
]

plot_traj!(unibots[1])
GM.activate!(); plotsave_energy(unibots[3])
function plotsave_velocity_restitution(bots,showlegend,
		figname=nothing;
		ymids = [1.0,0.7,0.3,0.0]
	)
	with_theme(theme_pub;
			resolution = (tw,0.4th),
			figure_padding = (fontsize,fontsize,0,0),

		) do
		fig = Figure()
		grids = [GridLayout(fig[i,j]) for i in 1:length(bots), j = 1:2]
		for (botid,bot) in enumerate(bots)
			ax1 = Axis(grids[botid,1][1,1], xlabel = tlabel, ylabel = L"v~(\mathrm{m/s})")
			ṙp2_mid = get_mid_velocity!(bot,1,2)
			t_mids = get_time_mids(bot)
			lines!(ax1,t_mids,ṙp2_mid[3,:], label="v́ₙ")
			lines!(ax1,t_mids,ṙp2_mid[1,:], label="v́ₜ" )
			if botid == 4 && showlegend
				axislegend(ax1;position=:cb,orientation=:horizontal)
			end
			xlims!(ax1,0,1.5)
			# ylims!(ax1,-6,ymax)
			Label(grids[botid,1][1,1,TopLeft()], "($(alphabet[2botid-1]))")
			ax2 = Axis(grids[botid,2][1,1], xlabel = "Count", ylabel = L"e_{\mathrm{eff}}")
			vz = get_velocity!(bot,1,2)[3,:]

			contacts_traj_voa = VectorOfArray(bot.contacts_traj)
			c1_traj = contacts_traj_voa[2,:]
			# steps = 1:length(c1_traj)
			idx_imp = findall(isimpact, c1_traj)
			@show idx_imp |> length
			e_eff = -vz[idx_imp]./vz[idx_imp.-1]
			scatter!(ax2,e_eff)
			# @show e_eff
			@show mean(e_eff)
			ylims!(ax2,ymids[botid]-0.08,ymids[botid]+0.08)
			# ax2.yticks = [0,0.3,0.7,1.0]
			# ax2.xticks = collect(1:length(e_eff))
			Label(grids[botid,2][1,1,TopLeft()], "($(alphabet[2botid]))")
		end
		rowgap!(fig.layout,-fontsize)
		savefig(fig,figname)
		fig
	end
end
GM.activate!(); plotsave_velocity_restitution(unibots,true)
CM.activate!(); plotsave_velocity_restitution(unibots,true,"unibot_restitution")

unibot_e5 = TR.solve!(
	TR.SimProblem(
		uni(0.0; μ = 0.01, e=0.5, z0 = 0.2),
		uni_dynfuncs
	),
	TR.ZhongCCP();
	tspan,dt=0.5h,ftol=1e-14,maxiters=50,exception=false,verbose=true,
)
GM.activate!(); plot_traj!(unibot_e5)

me = TR.mechanical_energy!(unibot_e5)
me.E |> lines

GM.activate!(); plotsave_friction_direction(
	[unibot_e5],L"e",[0.5]; 
	resolution = (0.9tw,0.28tw),
	cid=2
	)

CM.activate!(); plotsave_friction_direction(
		[unibot_e5],L"e=",[0.5], 
		"unibot_friction_direction"; 
		resolution = (0.9tw,0.28tw)
)

function plotsave_point_traj_vel(bots,figname=nothing)
	with_theme(theme_pub;
		resolution=(1.0tw,0.4tw),
	) do
		fig = Figure()
		nbot = length(bots)
		for (i,bot) in enumerate(bots)
			ax1 = Axis(
				fig[i,1],
				xlabel=L"x~(\mathrm{m})",
				ylabel=L"y~(\mathrm{m})",
				aspect=DataAspect(),
			)
			ax2 = Axis(
				fig[i,2],
				xlabel=tlabel,
				ylabel=L"\theta",
			)
			Label(fig[i,1,TopLeft()], "($(alphabet[2i-1]))")
			Label(fig[i,2,TopLeft()], "($(alphabet[2i]))")
			ax3 = Axis(fig[i,2],
						xlabel=tlabel,
						ylabel=L"v~(\mathrm{m/s})",
						yticklabelcolor = :red,
						yaxisposition = :right
				)
			hidespines!(ax3)
			hidexdecorations!(ax3)

			(;t) = bot.traj
			rp2 = get_trajectory!(bot,1,2)
			rpx = rp2[1,:]
			rpy = rp2[2,:]
			lines!(ax1,rpx,rpy)
			# xlims!(ax1,0,13.0)
			# ylims!(ax1,-0.3,0.45)

			ṙp2 = get_mid_velocity!(bot,1,2)
			ṙpx = ṙp2[1,:]
			ṙpy = ṙp2[2,:]
			t_mids = get_time_mids(bot)
			θ_mids = atan.(ṙpy,ṙpx)
			ṙpxy = map(ṙp2.u) do ṙ
				norm(ṙ[1:2])
			end
			lines!(ax2,t_mids,θ_mids)
			lines!(ax3,t_mids,ṙpxy,color=:red)
			xlims!(ax2,extrema(bot.traj.t)...)
			xlims!(ax3,extrema(bot.traj.t)...)
			ylims!(ax2,-π,π)

			if i !== nbot
				hidex(ax1)
				hidex(ax2)
			end

		end
		colsize!(fig.layout,1,Relative(0.35))
		savefig(fig,figname)
		fig
	end
end

GM.activate!(); plotsave_point_traj_vel([unibot_e5])
CM.activate!(); plotsave_point_traj_vel([unibot_e5],"unibot_contact_point_traj_vel")

#dt
dts = [1e-2,3e-3,1e-3,3e-4,1e-4,1e-5]
# dts = [1e-2]

unibot_z0 = [
	TR.solve!(
		TR.SimProblem(
			uni(0.0; μ = 0.01, e=0.0, z0 = 0.2-0.13468-1e-5, ωz = 50.0, mbar = 1.0),
			uni_dynfuncs
		),
		TR.ZhongCCP();
		tspan=(0.0,1.0),dt,ftol=1e-14,maxiters=50,exception=true
	)
	for dt in dts
]

plot_traj!(unibot_z0[6];
	# zlims = (0,1),
	showground = true
)

GM.activate!(); plotsave_error(unibot_z0,dts;resolution = (0.4tw,0.3tw),)
CM.activate!(); plotsave_error(unibot_z0,dts,"unibot_err";resolution = (0.4tw,0.3tw),)

GM.activate!(); plotsave_velocity_restitution([unibot_z0],true; ymids = [0.0])
plot_traj!(unibot_z0)

GM.activate!(); plotsave_point_traj_vel([unibot_z0])

me = TR.mechanical_energy!(unibot_z0)
me.E |> lines

GM.activate!(); plot_traj!(unibot_z0)



# -------- SUPERBall ---------
ballbot = superball(
	0.0;
	μ = 0.5,
	e = 0.5
)

GM.activate!(); with_theme(theme_pub;
		Scatter = (markersize = fontsize,)	
	) do
	plot_traj!(ballbot;
		AxisType = Axis3,
		figsize = (0.6tw,0.5tw),
		zlims = [1.0,3.0].-2.0,
		doslide = false,
		showinfo = false,
		# showlabels = false,
		showground = false,
		sup! = (ax,_,_) -> begin
			ax.zlabeloffset = 2fontsize
			ax.azimuth = 3.315530633326984
			ax.elevation = 0.2326990816987242
		end,
		# figname = "superball"
	)
end

GM.activate!(); plot_traj!(ballbot;
	xlims = [-1,10],
	ylims = [-1,3],
	zlims = [-1e-3,2.4],
)

function ball_dynfuncs(bot)
    (;tg) = bot
    function F!(F,q,q̇,t)
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_tensiles!(tg)
        TR.apply_gravity!(tg)
        F .= TR.generate_forces!(tg)
    end
    function Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t)
		∂F∂q̌ .= 0
        ∂F∂q̌̇ .= 0
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_tensiles!(tg)
        TR.build_∂Q̌∂q̌!(∂F∂q̌,tg)
        TR.build_∂Q̌∂q̌̇!(∂F∂q̌̇,tg)
    end

	bars = TR.get_rigidbars(tg)
	
	function prepare_contacts!(contacts,q)
		TR.update_rigids!(tg,q)
		foreach(bars) do rb
			rbid = rb.prop.id
			gap1 = rb.state.rps[1][3]
			TR.activate!(contacts[2rbid-1],gap1)
			gap2 = rb.state.rps[2][3]
			TR.activate!(contacts[2rbid  ],gap2)
		end
		active_contacts = filter(contacts) do c
			c.state.active
		end
		na = length(active_contacts)
		inv_μ_vec = ones(eltype(q),3na)
		n = [0,0,1.0]
		for (i,ac) in enumerate(active_contacts)
			(;state) = ac
			state.frame = TR.SpatialFrame(n)
			inv_μ_vec[3(i-1)+1] = 1/ac.μ
		end
        es = [ac.e for ac in active_contacts]
		gaps = [ac.state.gap for ac in active_contacts]
		H = Diagonal(inv_μ_vec)
        active_contacts, gaps, H, es
    end

	function get_directions_and_positions(active_contacts,q)
		na = length(active_contacts)
		TR.update_rigids!(tg,q)
		D = Matrix{eltype(q)}(undef,3na,length(q))
		ŕ = Vector{eltype(q)}(undef,3na)
		for (i,ac) in enumerate(active_contacts)
			(;id,state) = ac
			rbid,pid = divrem(ac.id-1,2) .+ 1
			(;n,t1,t2) = state.frame
			C = bars[rbid].state.cache.Cps[pid]
			CT = C*TR.build_T(tg,rbid)
			dm = hcat(n,t1,t2) |> transpose
			D[3(i-1)+1:3(i-1)+3,:] = dm*CT
			ŕ[3(i-1)+1:3(i-1)+3] = dm*bars[rbid].state.rps[pid]
		end
		D,ŕ
	end

    @eponymtuple(F!,Jac_F!,prepare_contacts!,get_directions_and_positions)
end

tspan = (0.0,2.0)
h = 1e-3

prob = TR.SimProblem(ballbot,ball_dynfuncs)

TR.solve!(prob,TR.ZhongCCP();tspan,dt=h,ftol=1e-14,maxiters=50,exception=false)

plot_traj!(ballbot)

GM.activate!(); 
with_theme(theme_pub;
		figure_padding = (0,1.5fontsize,fontsize,fontsize),
	) do
	plot_traj!(ballbot;
		AxisType = Axis3,
		xlims = [-1,10],
		ylims = [-1,3],
		zlims = [-1e-3,2.4],
		figsize = (1.0tw,0.6tw),
		gridsize = (3,2),
		attimes = [0,0.491,1.029,1.479,2.207,3.112],
		doslide = false,
		showinfo = false,
		showpoints = false,
		showlabels = false,
		sup! = (ax,_,_) -> begin
			ax.azimuth = 4.865530633326983
			ax.elevation = 0.2926990816987241
			ax.zlabeloffset = 2fontsize
		end,
		# figname = "ballbot"
	)
end

me = TR.mechanical_energy!(ballbot)
me.E |> lines


GM.activate!(); plotsave_energy(ballbot)
CM.activate!(); plotsave_energy(ballbot,"ballbot_energy")


function plotsave_contactpoints(bot,figname=nothing)
	contacts_traj_voa = VectorOfArray(bot.contacts_traj)
	(;t) = bot.traj
	with_theme(theme_pub;
			figure_padding = (0,fontsize,0,0),
			resolution = (0.9tw,0.3tw),
		) do 
		fig = Figure()
		ax = Axis(fig[1,1], xlabel = tlabel, ylabel = "Contact No.")
		markersize = fontsize
		for ic in eachindex(contacts_traj_voa[1])
			c_traj = contacts_traj_voa[ic,:]
			idx_imp = findall(isimpact, c_traj)
			idx_per = findall(doespersist,c_traj)
			scatter!(ax,t[idx_per],fill(ic,length(idx_per)); marker=:diamond, markersize)
			scatter!(ax,t[idx_imp],fill(ic,length(idx_imp)); marker=:xcross, markersize)
		end
		xlims!(ax,t[begin],t[end])
		ylims!(ax,0.5,12.5)
		ax.yticks = 1:12
		elem_1 = [MarkerElement(;color = :blue, marker = :xcross, markersize)]
		elem_2 = [MarkerElement(;color = :blue, marker = :diamond, markersize)]
		Legend(fig[1, 1],
			[elem_1, elem_2],
			["Impact", "Persistent"],
			tellheight = false,
			tellwidth = false,
			orientation = :horizontal, 
			halign = :right,
			valign = :top,
			margin = (0,fontsize,fontsize,2fontsize)
		)
		savefig(fig,figname)
		fig
	end
end

GM.activate!(); plotsave_contactpoints(ballbot)
CM.activate!(); plotsave_contactpoints(ballbot,"ballbot_contactpoints")

function plotsave_velocity_restitution(bots,showlegend,
		figname=nothing;
		cps,
		ymids = 0.5ones(length(cps)),
		resolution = (tw,0.5th),
	)
	with_theme(theme_pub;
			resolution,
			figure_padding = (fontsize,fontsize,0,0),
			fontsize = 6.5 |> pt2px
		) do
		fig = Figure()
		grids = [GridLayout(fig[i,j]) for i in 1:length(bots), j = 1:2]
		for (botid,bot) in enumerate(bots)
			ax1 = Axis(grids[botid,1][1,1], xlabel = tlabel, ylabel = L"v~(\mathrm{m/s})",)
			rbid,pid = divrem(cps[botid]-1,2) .+ 1
			ṙpc_mid = get_mid_velocity!(bot,rbid,pid)
			t_mids = get_time_mids(bot)
			lines!(ax1,t_mids,ṙpc_mid[3,:], label="v́ₙ")
			lines!(ax1,t_mids,ṙpc_mid[1,:], label="v́ₜ" )
			xlims!(ax1,0,10.0)
			# ylims!(ax1,-6,ymax)
			Label(grids[botid,1][1,1,TopLeft()], "($(alphabet[2botid-1]))")
			Label(grids[botid,1][1,1,Top()], "Contact No.$(cps[botid])", halign = :center)

			ax2 = Axis(grids[botid,2][1,1], xlabel = "Count", ylabel = L"e_{\mathrm{eff}}")
			vz = get_velocity!(bot,rbid,pid)[3,:]

			contacts_traj_voa = VectorOfArray(bot.contacts_traj)
			c1_traj = contacts_traj_voa[cps[botid],:]
			# steps = 1:length(c1_traj)
			idx_imp = findall(isimpact, c1_traj)
			@show idx_imp |> length
			e_eff = -vz[idx_imp]./vz[idx_imp.-1]
			scatter!(ax2,e_eff)
			@show mean(e_eff)
			ylims!(ax2,ymids[botid]-0.08,ymids[botid]+0.08)
			# ax2.yticks = [0,0.3,0.7,1.0]
			ax2.xticks = collect(1:length(e_eff))
			xlims!(ax2,0.5,length(e_eff)+0.5)
			Label(grids[botid,2][1,1,TopLeft()], "($(alphabet[2botid]))")

			Label(grids[botid,2][1,1,Top()], "Contact No.$(cps[botid])", halign = :center)

			if botid == 3 && showlegend
				axislegend(ax1;position=:cb,orientation=:horizontal)
			end
		end
		rowgap!(fig.layout,0)
		savefig(fig,figname)
		fig
	end
end
GM.activate!(); plotsave_velocity_restitution(
	[ballbot,ballbot,ballbot,ballbot],true;
	cps = [2,3,6,12],
)
CM.activate!(); plotsave_velocity_restitution(
	[ballbot,ballbot,ballbot,ballbot],true,"ballbot_restitution";
	cps = [2,3,6,12],
)


c6= VectorOfArray(ballbot.contacts_traj)[6,:]
c6_1 = c6[4163]
c6_1.state.Λ ⋅ c6_1.state.v
c6_1.state.v
check_Coulomb(1,c6_1)

GM.activate!(); plotsave_friction_direction(
		[ballbot,ballbot,ballbot],L"\mathrm{Cable~No.}",
		[3,6,12]; 
		resolution = (tw,0.5tw),
		mo_α = 8,
		mo = 8,
		vtol=1e-7,
)

CM.activate!(); plotsave_friction_direction(
		[ballbot,ballbot,ballbot],L"\mathrm{Cable~No.}",
		[3,6,12], 
		"ballbot_friction_direction"; 
		resolution = (tw,0.5tw),
		mo_α = 8,
		mo = 8,
		vtol=1e-7,
)

dts = [1e-2,3e-3,1e-3,3e-4,1e-4,1e-5]
superballs_dt = [
	begin
		prob = 
		TR.solve!(TR.SimProblem(superball(0.0;μ = 0.5,e = 0.5),
			ball_dynfuncs),
			TR.ZhongCCP();
			tspan=(0.0,0.12),dt,ftol=1e-14,
			maxiters=50,exception=false
		)
	end
	for dt in dts
]
GM.activate!(); plotsave_error(superballs_dt,dts,bid=5,pid=1)

# ----------- quadruped ---------------
quadbot = quad(10.0)
quadbot.tg.connectivity.numbered

plot_traj!(quadbot;
	zlims = (0,1),
	showground = false
)

function quad_dynfuncs(bot)
    (;tg) = bot
    function F!(F,q,q̇,t)
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_tensiles!(tg)
        TR.apply_gravity!(tg)
        F .= TR.generate_forces!(tg)
    end
    function Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t)
		∂F∂q̌ .= 0
        ∂F∂q̌̇ .= 0
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_tensiles!(tg)
        TR.build_∂Q̌∂q̌!(∂F∂q̌,tg)
        TR.build_∂Q̌∂q̌̇!(∂F∂q̌̇,tg)
    end
	bars = TR.get_rigidbars(tg)
	function prepare_contacts!(contacts, q)
		TR.update_rigids!(tg,q)
		foreach(bars) do rb
			(;id) = rb.prop
			gap = rb.state.rps[1][3]
			TR.activate!(contacts[id],gap)
		end
		active_contacts = filter(contacts) do c
			c.state.active
		end
		na = length(active_contacts)
		D = Matrix{eltype(q)}(undef,3na,length(q))
		inv_μ_vec = ones(eltype(q),3na)
		n = [0,0,1.0]
		for (i,ac) in enumerate(active_contacts)
			(;state) = ac
			state.frame = TR.SpatialFrame(n)
			(;n,t1,t2) = state.frame
			rbid = ac.id
			C = bars[rbid].state.cache.Cps[1]
			CT = C*TR.build_T(tg,rbid)
			Dn = Matrix(transpose(n)*CT)
			Dt1 = Matrix(transpose(t1)*CT)
			Dt2 = Matrix(transpose(t2)*CT)
			D[3(i-1)+1,:] = Dn
			D[3(i-1)+2,:] = Dt1
			D[3(i-1)+3,:] = Dt2
			inv_μ_vec[3(i-1)+1] = 1/ac.μ
		end
		es = [ac.e for ac in active_contacts]
		gaps = [ac.state.gap for ac in active_contacts]
		H = Diagonal(inv_μ_vec)
		active_contacts, na, gaps, D, H, es
	end

    F!,Jac_F!,prepare_contacts!
end

tspan = (0.0,1.91)
tspan = (0.0,1.0)
h = 1e-3

prob = TR.SimProblem(quadbot,quad_dynfuncs)

TR.solve!(prob,TR.ZhongCCP();tspan,dt=h,ftol=1e-14,maxiters=50,exception=true)

TR.solve!(prob,TR.AlphaCCP(0.8);tspan,dt=h,ftol=1e-10,maxiters=50,exception=true)

plot_traj!(quadbot;
	zlims = (0,1),
	showground = false
)

me = TR.mechanical_energy!(quadbot)

me.E |> lines
contacts_traj_voa = VectorOfArray(quadbot.contacts_traj)
