using LinearAlgebra, Statistics
using StaticArrays
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using LaTeXStrings
using RecursiveArrayTools
using StructArrays
using BenchmarkTools
using TypeSortedCollections
using CoordinateTransformations
using Rotations
using GeometryBasics
using Peaks
using OffsetArrays
using Printf
using Unitful
using Match
using Revise
import Meshes
using Meshing
import TensegrityRobots as TR
cd("examples/nonsmooth")
includet("../analysis.jl")
includet("../vis.jl")
const figdir=raw"C:\Users\luo22\OneDrive\Papers\Ph.D.Thesis\ns"

function make_top(ro = [0.0,0.0,0.0],
       Ro = one(RotMatrix{3}),
       ṙo = [0.0,0.0,0.0],
       ωo = [0.0,0.0,5.0];
	   μ = 0.5,
	   e = 0.9
	)
    movable = true
    constrained = false
    m = 1.0
    r̄g = @SVector zeros(3)
    Ī = SMatrix{3,3}(Matrix(1.0I,3,3))
    h = 0.5
    r̄ps = [h.*[x,y,z] for x = [-1,1] for y = [-1,1] for z = [1]]
    push!(r̄ps,[0,0,-h])
	pts = Point3.(r̄ps)
	fcs = TriangleFace.([
		[5,1,2],
		[5,4,3],
		[5,3,1],
		[5,2,4],
		[1,4,2],
		[4,1,3],
		[3,2,1],
		[2,3,4]
	])
	nls = normals(pts,fcs)
	top_mesh = GeometryBasics.Mesh(meta(pts,normals=nls),fcs)
    prop = TR.RigidBodyProperty(1,movable,m,Ī,r̄g,r̄ps)
    ri = ro
    lncs, _ = TR.NaturalCoordinates.NC1P3V(ri,ro,Ro,ṙo,ωo)
    state = TR.RigidBodyState(prop,lncs,ro,Ro,ṙo,ωo)
    rb1 = TR.RigidBody(prop,state,top_mesh)
	rbs = TypeSortedCollection((rb1,))
	numberedpoints = TR.number(rbs)
	matrix_sharing = zeros(Int,0,0)
	indexedcoords = TR.index(rbs,matrix_sharing)
    ss = Vector{Int}()
	tensiles = (cables = ss,)
	connections = TR.connect(rbs,zeros(Int,0,0))
	cnt = TR.Connectivity(numberedpoints,indexedcoords,connections)
	contacts = [TR.Contact(i,μ,e) for i = [5]]
	tg = TR.TensegrityStructure(rbs,tensiles,cnt,contacts)
    bot = TR.TensegrityRobot(tg)
end

function top_contact_dynfuncs(bot)
    (;tg) = bot
    function F!(F,q,q̇,t)
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_cables_apply_forces!(tg)
        TR.apply_gravity!(tg)
        F .= TR.generate_forces!(tg)
    end
    function Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t)
		∂F∂q̌ .= 0
        ∂F∂q̌̇ .= 0
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_cables_apply_forces!(tg)
        TR.build_∂Q̌∂q̌!(∂F∂q̌,tg)
        TR.build_∂Q̌∂q̌̇!(∂F∂q̌̇,tg)
    end

	rbs = TR.get_rigidbodies(tg)

    function prepare_contacts!(contacts,q)
		TR.update_rigids!(tg,q)
		rb = rbs[1]
		for (cid,pid) in enumerate([5])
			gap = rb.state.rps[pid][3]
			TR.activate!(contacts[cid],gap)
		end
		active_contacts = filter(contacts) do c
			c.state.active
		end
		na = length(active_contacts)
		D = Matrix{eltype(q)}(undef,3na,length(q))
		inv_μ_vec = ones(eltype(q),3na)
		n = [0,0,1.0]
		for (i,ac) in enumerate(active_contacts)
			(;id,state) = ac
            state.frame = TR.SpatialFrame(n)
			(;n,t1,t2) = state.frame
			rbid = 1
            C = rbs[rbid].state.cache.Cps[id]
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

# Spinning
ro = [0,0,2.0]
Ro = RotX(π/24)
ṙo = [1.0,0.0,0.0]
ωo = [0.0,0.0,5.0]
top = make_top(ro,Ro,ṙo,ωo;μ = 0.95)

tspan = (0.0,20.0)
h = 1e-3

prob = TR.SimProblem(top,top_contact_dynfuncs)

TR.solve!(prob,TR.ZhongCCP();tspan,dt=h,ftol=1e-14,maxiters=50,exception=true)

plot_traj!(top;showinfo=false,rigidcolor=:white,showwire=true,figsize=(0.6tw,0.6tw))

contacts_traj_voa = VectorOfArray(top.contacts_traj)
for (i,c) in enumerate(contacts_traj_voa[1,1:end])
	check_Coulomb(i,c)
end


function plotsave_energy(bot,figname=nothing)
	(;traj) = bot
	(;t) = traj
	tmin,tmax = t[begin], t[end]
	titles = [
		"机械能",
		"动能",
		"势能"
	]
	with_theme(mv_theme; resolution = (0.9tw,0.6tw)) do
		ME = TR.mechanical_energy!(bot;gravity=true)
		fig = Figure()
		axs = [
			begin
				ax = Axis(
					fig[i,1],
					xlabel=L"t~\mathrm{(s)}",
					ylabel="Energy (J)",
					title=titles[i]
				)
				xlims!(ax,tmin,tmax)
				# ylims!(ax,,)
				if i !== 3
					ax.xlabelvisible = false
					ax.xticklabelsvisible = false
				end
				Label(fig[i,1,TopLeft()], "($(alphabet[i]))")
				ax
			end
			for i = 1:3
		]

		lines!(axs[1], t, ME.E)
		lines!(axs[2], t, ME.T)
		lines!(axs[3], t, ME.V)


		savefig(fig,figname)

		fig
	end
end
GM.activate!(); plotsave_energy(top)
CM.activate!(); plotsave_energy(top,"spinningtop_energy")

# contacts_traj = TR.solve!(prob,TR.AlphaCCP(0.95);tspan,dt=h,ftol=1e-8,maxiters=50,exception=true)

with_theme(my_theme;
		Axis3 = (
			azimuth = 4.2555306333269835,
			elevation = 0.2326990816987238
		)
	) do
	plot_traj!(top;
		doslide=false,
		gridsize=(2,3),
		attimes=[0,2,4,6,8,10],
		xlims=(-0.5,1.5),
		zlims=(-1e-3,3),
		rigidcolor=:white,
		showwire=true,
		showlabels=false,
		fontsize=8 |> pt2px,
		figsize=(0.9tw,0.7tw),
		savefig=true,
		# figname="spinningtop_traj"
	)
end

function plotsave_contact_persistent(bot,figname=nothing)
	with_theme(mv_theme;) do
		contacts_traj_voa = VectorOfArray(bot.contacts_traj)
		c1_traj = contacts_traj_voa[1,:]
		Λtol = 1e-7 # rule out false active
		# steps = 1:length(c1_traj)
		active_nonpersist = findall(c1_traj) do c
			c.state.active && !c.state.persistent && norm(c.state.Λ) > Λtol
		end
		active_persist = findall(c1_traj) do c
			c.state.active && c.state.persistent && norm(c.state.Λ) > Λtol
		end
		fig = Figure()
		ax = Axis(fig[1,1], xlabel = "Step", ylabel = "Contact Type")
		scatter!(ax,active_nonpersist,one.(active_nonpersist))
		scatter!(ax,active_persist,zero.(active_persist))
		savefig(fig,figname)
		fig
	end
end
GM.activate!(); plotsave_contact_persistent(top)
CM.activate!(); plotsave_contact_persistent(top,"spinningtop_contact_persistent")

es = [1.0,0.7,0.3,0.0]
tops = [
	begin
		top = make_top(ro,Ro,[10.0,0.0,0.0],ωo; μ, e)
		TR.solve!(TR.SimProblem(top,top_contact_dynfuncs),
				TR.ZhongCCP();
				tspan=(0.0,4.0),dt=1e-3,ftol=1e-14,maxiters=50,exception=false)
	end
	for μ in [0.9,0.1] for e in es
]

GM.activate!(); plotsave_contact_persistent(tops[3])

function plotsave_velocity_restitution(bots,ymax,showlegend,figname=nothing)
	with_theme(mv_theme;
			resolution = (0.9tw,0.4th)
		) do
		fig = Figure()
		grids = [GridLayout(fig[i,j]) for i in 1:length(bots), j = 1:2]
		ymids = [1.0,0.7,0.3,0.0]
		for (botid,bot) in enumerate(bots)
			ax1 = Axis(grids[botid,1][1,1], xlabel = tlabel, ylabel = L"v~(\mathrm{m/s})")
			ṙp5_mid = get_mid_velocity!(bot,1,5)
			t_mids = get_time_mids(bot)
			lines!(ax1,t_mids,ṙp5_mid[3,:], label="v́ₙ")
			lines!(ax1,t_mids,ṙp5_mid[1,:], label="v́ₜ" )
			if botid == 4 && showlegend
				axislegend(ax1;position=:cb,orientation=:horizontal)
			end
			xlims!(ax1,0,4)
			ylims!(ax1,-6,ymax)
			Label(grids[botid,1][1,1,TopLeft()], "($(alphabet[2botid-1]))")
			ax2 = Axis(grids[botid,2][1,1], xlabel = "Count", ylabel = L"e_{\mathrm{eff}}")
			vz = get_velocity!(bot,1,5)[3,:]

			contacts_traj_voa = VectorOfArray(bot.contacts_traj)
			c1_traj = contacts_traj_voa[1,:]
			# steps = 1:length(c1_traj)
			idx_imp = findall(isimpact, c1_traj)
			e_eff = -vz[idx_imp]./vz[idx_imp.-1]
			scatter!(ax2,e_eff)
			@show mean(e_eff)
			ylims!(ax2,ymids[botid]-0.08,ymids[botid]+0.08)
			# ax2.yticks = [0,0.3,0.7,1.0]
			ax2.xticks = collect(1:length(e_eff))
			Label(grids[botid,2][1,1,TopLeft()], "($(alphabet[2botid]))")
		end
		rowgap!(fig.layout,-fontsize)
		savefig(fig,figname)
		fig
	end
end
GM.activate!()
plotsave_velocity_restitution(tops[1:4],6,true)
plotsave_velocity_restitution(tops[5:8],11,false)
CM.activate!()
plotsave_velocity_restitution(tops[1:4],6,true,"top_restitution_09")
plotsave_velocity_restitution(tops[5:8],11,false,"top_restitution_01")

plot_traj!(tops[8];showinfo=false,rigidcolor=:white,showwire=true,figsize=(0.6tw,0.6tw),
	xlims=(-1.0,100.0),
	ylims=(-1.0,1.0),
	zlims=(-1e-3,1.0),
)

me = TR.mechanical_energy!(tops[1])
me.E |> lines

# friction_direction

# Dare you try tspan = (0.0,15.0)

μs = [0.02,0.01,0.004]
tops_e0 = [
	begin
		top = make_top(ro,Ro,ṙo,ωo; μ, e = 0.0)
		TR.solve!(TR.SimProblem(top,top_contact_dynfuncs),
				TR.ZhongCCP();
				tspan=(0.0,20.0),dt=1e-3,ftol=1e-14,maxiters=50,exception=false)
	end
	for μ in μs
]

function plotsave_friction_direction(bots,x,xs,figname=nothing)
	with_theme(mv_theme;
		resolution = (0.9tw,0.4tw)
	) do
		fig = Figure()

		for (botid,bot) in enumerate(bots)
			ax1 = Axis(fig[botid,1],
						xlabel=tlabel,
						ylabel=L"\delta\alpha~(\mathrm{Rad})",
						title=latexstring("$x=$(xs[botid])")
				)
			ax2 = Axis(fig[botid,2],
						xlabel=tlabel,
						ylabel=L"\delta\theta~(\mathrm{Rad})",
						title=latexstring("$x=$(xs[botid])")
				)
			Label(fig[botid,1,TopLeft()],"($(alphabet[2botid-1]))")
			Label(fig[botid,2,TopLeft()],"($(alphabet[2botid]))")
			(;t) = bot.traj
			contacts_traj_voa = VectorOfArray(bot.contacts_traj)
			c1s = contacts_traj_voa[1,:]
			idx_sli = findall(c1s) do c
						issliding(c;vtol=1e-5)
					end
			idx_imp = findall(isimpact,c1s) ∩ idx_sli
			δα_imp = map(c1s[idx_imp]) do c
						get_contact_angle(c)
					end
			idx_per = findall(doespersist,c1s) ∩ idx_sli
			δα_per = map(c1s[idx_per]) do c
						get_contact_angle(c)
					end
			mo = 8
			scaling = 10.0^(-mo)
			Label(fig[botid,1,Top()],latexstring("\\times 10^{-$(mo)}"))
			scatter!(ax1,t[idx_per],δα_per./scaling)
			scatter!(ax1,t[idx_imp],δα_imp./scaling)
			ylims!(ax1,-1.0,1.0)
			xlims!(ax1,extrema(t)...)

			θ_imp = map(c1s[idx_imp]) do c
					get_friction_direction(c)
				end
			θ_per = map(c1s[idx_per]) do c
					get_friction_direction(c)
				end
			mo = 8
			scaling = 10.0^(-mo)
			Label(fig[botid,2,Top()],latexstring("\\times 10^{-$(mo)}"))
			scatter!(ax2,t[idx_per],(θ_per.-π/4)./scaling, label="Persistent")
			scatter!(ax2,t[idx_imp],(θ_imp.-π/4)./scaling, label="Impact")
			xlims!(ax2,extrema(t)...)
			ylims!(ax2,-1.0,1.0)
			if botid !== length(bots)
				hidex(ax1)
				hidex(ax2)
			else
				Legend(fig[length(bots)+1,:],ax2;
					orientation=:horizontal
				)
			end
		end
		rowgap!(fig.layout,fontsize/2)
		savefig(fig,figname)
		fig
	end
end
δα = map(tops[6].contacts_traj) do c
		get_contact_angle(c[1])
	end |> skipmissing |> collect |> scatter
θ = map(tops[8].contacts_traj) do c
		get_friction_direction(c[1])
	end |> skipmissing |> collect |> scatter

# GM.activate!(); plotsave_friction_direction(tops[1:4])
GM.activate!(); plotsave_friction_direction(tops[5:7],"e",es[1:3])
CM.activate!(); plotsave_friction_direction(tops[5:7],"e",es[1:3],"friction_direction_mu01")

GM.activate!(); plotsave_friction_direction(tops_e0,"\\mu",μs)
CM.activate!(); plotsave_friction_direction(tops_e0,"\\mu",μs,"friction_direction_e00")

# contact point trajectory


function plotsave_point_traj_vel(bots,figname=nothing)
	with_theme(mv_theme;
		resolution=(0.9tw,0.4tw),
	) do
		fig = Figure()
		nbot = length(bots)
		for (i,bot) in enumerate(bots)
			ax1 = Axis(
				fig[i,1],
				xlabel=L"x~(\mathrm{m})",
				ylabel=L"y~(\mathrm{m})",
				title=latexstring("\\mu = $(μs[i])")
			)
			ax2 = Axis(
				fig[i,2],
				xlabel=tlabel,
				ylabel=L"\theta",
				title=latexstring("\\mu = $(μs[i])")
			)
			Label(fig[i,1,TopLeft()], "($(alphabet[2i-1]))")
			Label(fig[i,2,TopLeft()], "($(alphabet[2i]))")

			(;t) = bot.traj
			rp5 = get_trajectory!(bot,1,5)
			rpx = rp5[1,:]
			rpy = rp5[2,:]
			lines!(ax1,rpx,rpy)
			xlims!(ax1,0,12.5)
			ylims!(ax1,-0.25,0.25)

			ṙp5 = get_mid_velocity!(bot,1,5)
			ṙpx = ṙp5[1,:]
			ṙpy = ṙp5[2,:]
			t_mids = get_time_mids(bot)
			θ_mids = atan.(ṙpy,ṙpx)
			idx_v = findall(ṙp5.u) do ṙ
						norm(ṙ[2:3]) > 1e-4
					end
			lines!(ax2,t_mids[idx_v],θ_mids[idx_v])
			# lines!(ax2,t_mids,ṙpx, label="ẋ")
			# lines!(ax2,t_mids,ṙpy, label="ẏ")
			xlims!(ax2,extrema(t_mids)...)
			ylims!(ax2,-π,π)

			if i !== nbot
				hidex(ax1)
				hidex(ax2)
			end

			vlines!(ax1,t[554],)
			vlines!(ax2,t[554],)
		end
		savefig(fig,figname)
		fig
	end
end


GM.activate!(); plotsave_point_traj_vel(tops_e0)
CM.activate!(); plotsave_point_traj_vel(tops_e0,"contact_point_traj_vel")

ωs = [5.0,10.0,20.0]
tops_ω = [
	begin
		top = make_top([0,0,cos(π/24)*0.5-1e-14],Ro,[0,0,0.0],[0,0,ω]; μ=0.9, e = 0.0)
		TR.solve!(TR.SimProblem(top,top_contact_dynfuncs),
				TR.ZhongCCP();
				tspan=(0.0,500.0),dt=1e-2,ftol=1e-14,maxiters=50,exception=false)
	end
	for ω in ωs
]

plot_traj!(tops_ω[1];showinfo=false,rigidcolor=:white,showwire=true,figsize=(0.6tw,0.6tw),
	xlims=(-1.0,2.0),
	ylims=(-1.0,1.0),
	zlims=(-1e-3,1.0),
)

me = TR.mechanical_energy!(tops_ω500s[1])
me.E[100:end]./me.E[end] |> lines

# energy conserving
function plotsave_energy_conserving(bots,figname=nothing)
	with_theme(mv_theme;
			resolution = (0.7tw,0.5tw)
		) do
		fig = Figure()

		for (botid,bot) in enumerate(bots)
			ax1 = Axis(
				fig[botid,1],
				xlabel=tlabel,
				ylabel="Rel. Err.",
				title=latexstring("\\omega_z=$(ωs[botid])")
			)
			Label(fig[botid,1,TopLeft()],"($(alphabet[botid]))")
			if botid !== 3
				hidex(ax1)
			end
			pre = 500
			(;t) = bot.traj
			mo = 3
			scaling = 10.0^(-mo)
			Label(fig[botid,1,Top()],latexstring("\\times 10^{-$mo}"))
			me = TR.mechanical_energy!(bot)
			lines!(ax1,t[pre:end],(me.E[pre:end].-me.E[pre])./me.E[pre]./scaling)
			xlims!(ax1,t[pre],t[end])
			ylims!(ax1,-1,1)
		end
		savefig(fig,figname)
		fig
	end
end
GM.activate!(); plotsave_energy_conserving(tops_ω)
CM.activate!(); plotsave_energy_conserving(tops_ω,"energy_conserving")

# friction_direction
# Dare you try
GM.activate!(); plotsave_friction_direction(tops_e0)

# contact_friction
function plotsave_contact_friction(bots,figname=nothing)
	with_theme(mv_theme;
		resolution = (0.9tw,0.5tw)
	) do
		fig = Figure()

		axs = [
			begin
				k = 3(j-1)+i
				ax = Axis(
					fig[i,j],
					xlabel=tlabel,
					ylabel=L"\Lambda~(\mathrm{N}\cdot\mathrm{s})",
					title=latexstring("\\mu=$(μs[k])")
				)
				if i !== 3
					hidex(ax)
				end
				Label(fig[i,j,TopLeft()],"($(alphabet[k]))")
				# ylims!(ax,16.5,18)
				ax
			end
			for i = 1:3, j = 1:2
		]

		for (botid,bot) in enumerate(bots)
			ax = axs[botid]
			(;t) = bot.traj
			c1s = VectorOfArray(bot.contacts_traj)[1,:]
			idx_per = findall(doespersist, c1s)
			Λt1 = [c.state.Λ[2] for c in c1s]
			Λt2 = [c.state.Λ[3] for c in c1s]
			lines!(ax,t[idx_per],Λt1[idx_per],label=L"\Lambda_{t_1}")
			lines!(ax,t[idx_per],Λt2[idx_per],label=L"\Lambda_{t_2}")
			# pos_normal_idx = findall((x)->x[1]>1e-4,c1_states.Λ)
			# active_pos_idx = intersect(idx_per,pos_normal_idx)
			# t_mids = 0.5 .* (t[active_pos_idx] .+ t[active_pos_idx .+ 1])
			# Λt1_mids = 0.5 .* (Λt1[active_pos_idx] .+ Λt1[active_pos_idx .+ 1])
			# Λt2_mids = 0.5 .* (Λt2[active_pos_idx] .+ Λt2[active_pos_idx .+ 1])
			# lines!(ax,t_mids,Λt1_mids,label=L"\Lambda_{t_1}")
			# lines!(ax,t_mids,Λt2_mids,label=L"\Lambda_{t_2}")
			xlims!(ax,t[556],t[end])
		end
		Legend(fig[:,3],axs[1])
		savefig(fig,figname)
		fig
	end
end
GM.activate!(); plotsave_contact_friction(tops_e0)
CM.activate!(); plotsave_contact_friction(tops_e0,"contact_friction")



# Bar

function make_bar(;
		μ = 0.1,
		e = 0.0
	)
    movable = true
    constrained = false
    m = 0.402
    l = 1.0
    Ixx = 1/12*m*l^2
    Ī = SMatrix{3,3}(Diagonal([Ixx,0.0,0.0]))
    r̄g = SVector{3}(0.0,0.0,0.0)
    r̄ps = SVector{3}.([
        [ -l/2, 0.0, 0.0],
        [  l/2, 0.0, 0.0]
    ])
    θ = π/4
	zoffset = 1e-6
	xoffset = 1e-6
    ri = [     0.0-xoffset,0.0,sin(θ)*l-zoffset]
    rj = [cos(θ)*l-xoffset,0.0,     0.0-zoffset]
    ro = (ri + rj)./2
    ṙo = zero(ro)
    ω = zero(ro)
    R = Matrix(RotY(θ))

	# b = l/2*ω[2]*sin(θ)-9.81
	# p⁺ = 1 + 3*cos(θ)^2-μ*cos(θ)*sin(θ)
	# @show b,p⁺

    prop = TR.RigidBodyProperty(1,movable,m,Ī,r̄g,r̄ps;constrained)

    lncs, _ = TR.NaturalCoordinates.NC3D2P(ri,rj,ro,R,ṙo,ω)
    state = TR.RigidBodyState(prop,lncs,ro,R,ṙo,ω)

	p1 = Meshes.Point(r̄ps[1])
	p2 = Meshes.Point(r̄ps[2])
	s = Meshes.Segment(p1,p2)
	cyl_bar = Meshes.Cylinder(Meshes.length(s)/40,s)
    cylsurf_bar = Meshes.boundary(cyl_bar)
	cyl_bar_simple = Meshes.triangulate(cylsurf_bar)
	cyl_bar_mesh = cyl_bar_simple |> simple2mesh
    rb1 = TR.RigidBody(prop,state,cyl_bar_mesh)
	rbs = TypeSortedCollection((rb1,))
	numberedpoints = TR.number(rbs)
	matrix_sharing = zeros(Int,0,0)
	indexedcoords = TR.index(rbs,matrix_sharing)
    ss = Vector{Int}()
	tensiles = (cables = ss,)
	connections = TR.connect(rbs,zeros(Int,0,0))
	contacts = [TR.Contact(i,μ,e) for i = 1:2]
	cnt = TR.Connectivity(numberedpoints,indexedcoords,connections)
	tg = TR.TensegrityStructure(rbs,tensiles,cnt,contacts)
    bot = TR.TensegrityRobot(tg)
end

bar = make_bar(;)

plot_traj!(bar;)

function bar_contact_dynfuncs(bot)
    (;tg) = bot
    function F!(F,q,q̇,t)
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_cables_apply_forces!(tg)
        TR.apply_gravity!(tg)
        F .= TR.generate_forces!(tg)
    end
    function Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t)
		∂F∂q̌ .= 0
        ∂F∂q̌̇ .= 0
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_cables_apply_forces!(tg)
        TR.build_∂Q̌∂q̌!(∂F∂q̌,tg)
        TR.build_∂Q̌∂q̌̇!(∂F∂q̌̇,tg)
    end

	rbs = TR.get_rigidbodies(tg)
	function prepare_contacts!(contacts, q)
		TR.update_rigids!(tg,q)
		rb = rbs[1]
		dirs = [1,3]
		for pid = 1:2
			gap = rb.state.rps[pid][dirs[pid]]
			TR.activate!(contacts[pid],gap)
		end
		active_contacts = filter(contacts) do c
			c.state.active
		end
		na = length(active_contacts)
		D = Matrix{eltype(q)}(undef,3na,length(q))
		inv_μ_vec = ones(eltype(q),3na)
		for (i,ac) in enumerate(active_contacts)
			(;id,state) = ac
			n = float.(1:3 .== dirs[id])
			state.frame = TR.SpatialFrame(n)
			(;n,t1,t2) = state.frame
			rbid = 1
			C = rbs[rbid].state.cache.Cps[id]
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


tspan = (0.0,0.4)
h = 1e-3


prob = TR.SimProblem(bar,bar_contact_dynfuncs)

TR.solve!(prob,TR.ZhongCCP();tspan,dt=h,ftol=1e-12,maxiters=50,exception=false)

contacts_traj_voa = VectorOfArray(bar.contacts_traj)
for (i,c) in enumerate(contacts_traj_voa[1,1:end])
	check_Coulomb(i,c)
end

plot_traj!(bar;showinfo=false)

me = TR.mechanical_energy!(bar)
lines(me.E)

[c[1].state.active for c in bar.contacts_traj] |> lines
[c[1].state.persistent for c in bar.contacts_traj] |> lines

[c[2].state.active for c in bar.contacts_traj] |> lines
[c[2].state.persistent for c in bar.contacts_traj] |> lines

# pointmass

function new_pointmass(;
		ro = [0.0,0.0,1.0],
		ṙo = zero(ro),
		m = 1.0,
		μ = 0.3,
		e = 0.9
	)
	movable = false
	constrained = true
	Ia = SMatrix{3,3}(Matrix(m*I,3,3))
	r̄g  = SVector{3}([ 0.0, 0.0, 0.0])
	r̄p1 = SVector{3}([ 0.0, 0.0, 0.0])
	r̄ps = [r̄p1]
	prop = TR.RigidBodyProperty(1,movable,m,Ia,
				r̄g,r̄ps;constrained=constrained
				)
	R = RotX(0.0)
	ω = zero(ro)
	rps = Ref(ro) .+ Ref(R).*r̄ps
	lncs, _ = TR.NaturalCoordinates.NCMP(rps,ro,R,ṙo,ω)
	ci = Vector{Int}()
	Φi = Vector{Int}()
	state = TR.RigidBodyState(prop,lncs,ro,R,ṙo,ω,ci,Φi)
	rb1 = TR.RigidBody(prop,state)

	rbs = TypeSortedCollection((rb1,))
	numberedpoints = TR.number(rbs)
	matrix_sharing = zeros(Int,0,0)
	indexedcoords = TR.index(rbs,matrix_sharing)
	ss = Vector{Int}()
	tensiles = (cables = ss,)
	connections = TR.connect(rbs,zeros(Int,0,0))
	cnt = TR.Connectivity(numberedpoints,indexedcoords,connections)
	contacts = [TR.Contact(i,μ,e) for i = [1]]
	tg = TR.TensegrityStructure(rbs,tensiles,cnt,contacts)
	bot = TR.TensegrityRobot(tg)
end

function pm_contact_dynfuncs(bot;θ=0.0)
    (;tg) = bot
    function F!(F,q,q̇,t)
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_cables_apply_forces!(tg)
        TR.apply_gravity!(tg)
        F .= TR.generate_forces!(tg)
    end
    function Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t)
		∂F∂q̌ .= 0
        ∂F∂q̌̇ .= 0
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_cables_apply_forces!(tg)
        TR.build_∂Q̌∂q̌!(∂F∂q̌,tg)
        TR.build_∂Q̌∂q̌̇!(∂F∂q̌̇,tg)
    end

	rbs = TR.get_rigidbodies(tg)

	a = tan(θ)
	n = [-a,0,1]
	inclined_plane = TR.Plane(n,zeros(3))
	function prepare_contacts!(contacts,q)
		T = eltype(q)

        g = TR.signed_distance(q,inclined_plane)

		TR.activate!(contacts[1],g)

		active_contacts = filter(contacts) do c
			c.state.active
		end

		na = length(active_contacts)
		D = Matrix{T}(undef,3na,length(q))
		inv_μ_vec = ones(T,3na)
		for (i,ac) in enumerate(active_contacts)
			(;state) = ac
			state.frame = TR.SpatialFrame(n)
			(;n,t1,t2) = state.frame
			rbid = ac.id
			C = rbs[rbid].state.cache.Cps[1]
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

tspan = (0.0,0.455)
tspan = (0.0,4.0)
h = 1e-3

# horizontal plane
es = [1.0,0.5,0.0]
v0s = [0.0,1.0]
pms_hp = [
	begin
		pm = new_pointmass(;
				e,μ=0.1,ṙo = [v0,0,0]
			)
		TR.solve!(TR.SimProblem(pm,pm_contact_dynfuncs),TR.ZhongCCP();tspan,dt=h,ftol=1e-14,maxiters=50,exception=false)
	end
	for v0 in v0s for e in es
]

CM.activate!()
with_theme(my_theme;
		fontsize = 6 |> pt2px,
		Axis3 = (
			azimuth = 4.575530633326986,
			elevation = 0.13269908169872408
		)
	) do
	plot_traj!(reshape(pms_hp,3,2)[2,2];
		doslide=false,
		gridsize=(1,4),
		attimes=[0,0.3,0.455,0.5],
		showmesh=false,
		xlims=(-1e-3,1.0),
		ylims=(-0.5,0.5),
		zlims=(-1e-3,1.0),
		showlabels=false,
		figsize=(0.9tw,0.25tw),
		# figname="pointmass_e05_v00",
		figname="pointmass_e05_v10"
	)
end

function plotsave_z_energy(bots,figname=nothing)
	with_theme(mv_theme;
			resolution = (0.9tw,0.45tw),
		) do
		fig = Figure()
		grids = [GridLayout(fig[i,j]) for i in 1:length(bots), j = 1:2]
		# xmaxs = [20,4,2,1]
		# ymids = [1.0,0.7,0.3,0.0]
		# xtickmaxs = [20,12,4,1]
		for (botid,bot) in enumerate(bots)
			(;t) = bot.traj
			ax1 = Axis(grids[botid,1][1,1], xlabel = tlabel, ylabel = L"z~(\mathrm{m})")
			rp1 = get_trajectory!(bot,1,1)
			lines!(ax1,t,rp1[3,:])
			xlims!(ax1,t[begin],t[end])
			# ylims!(ax1,-6,6)
			Label(grids[botid,1][1,1,TopLeft()], "($(alphabet[2botid-1]))")
			ax2 = Axis(grids[botid,2][1,1], xlabel = tlabel, ylabel = L"\mathrm{Energy}~(\mathrm{J})")
			me = TR.mechanical_energy!(bot)
			lines!(ax2,t,me.E)
			xlims!(ax2,t[begin],t[end])
			ylims!(ax2,-1,11)
			# ax2.yticks = [0,0.3,0.7,1.0]
			# ax2.xticks = collect(1:xtickmaxs[botid])
			if botid !== 4
				hidex(ax1); hidex(ax2)
			end
			Label(grids[botid,2][1,1,TopLeft()], "($(alphabet[2botid]))")
		end
		savefig(fig,figname)
		fig
	end
end

GM.activate!(); plotsave_z_energy(reshape(pms_hp,3,2)[:,1])
CM.activate!(); plotsave_z_energy(reshape(pms_hp,3,2)[:,1],"pointmass_z_energy")

function plotsave_pointmass_xz_energy(bots,figname=nothing)
	with_theme(mv_theme;
			resolution = (0.9tw,0.45tw),
		) do
		fig = Figure()
		grids = [GridLayout(fig[i,j]) for i in 1:length(bots), j = 1:2]
		for (botid,bot) in enumerate(bots)
			(;t) = bot.traj
			ax1 = Axis(grids[botid,1][1,1], xlabel = tlabel, ylabel = "Coords. (m)")
			rp1 = get_trajectory!(bot,1,1)
			lines!(ax1,t,rp1[3,:],label=L"z")
			lines!(ax1,t,rp1[1,:],label=L"x")
			if botid == 2
				axislegend(ax1,position=:ct,orientation=:horizontal)
			end
			xlims!(ax1,t[begin],t[end])
			Label(grids[botid,1][1,1,TopLeft()], "($(alphabet[2botid-1]))")
			ax2 = Axis(grids[botid,2][1,1], xlabel = tlabel, ylabel = L"\mathrm{Energy}~(\mathrm{J})")
			me = TR.mechanical_energy!(bot)
			lines!(ax2,t,me.E)
			xlims!(ax2,t[begin],t[end])
			ylims!(ax2,-1,11)
			if botid !== 4
				hidex(ax1); hidex(ax2)
			end
			Label(grids[botid,2][1,1,TopLeft()], "($(alphabet[2botid]))")
		end
		savefig(fig,figname)
		fig
	end
end
GM.activate!(); plotsave_pointmass_xz_energy(reshape(pms_hp,3,2)[:,2])
CM.activate!(); plotsave_pointmass_xz_energy(reshape(pms_hp,3,2)[:,2],"pointmass_xz_energy")

function plotsave_pointmass_velocity_restitution(bots,figname=nothing)
	with_theme(mv_theme;
			resolution = (0.9tw,0.6tw),
		) do
		fig = Figure()
		grids = [GridLayout(fig[i,j]) for i in 1:length(bots), j = 1:2]
		ymids = es
		xtickmaxs = [4,5,1]
		for (botid,bot) in enumerate(bots)
			(;t) = bot.traj
			ax1 = Axis(grids[botid,1][1,1], xlabel = tlabel, ylabel = L"v~(\mathrm{m/s})")
			ṙp1 = get_velocity!(bot,1,1)
			lines!(ax1,t,ṙp1[3,:])
			xlims!(ax1,t[begin],t[end])
			ylims!(ax1,-6,6)
			Label(grids[botid,1][1,1,TopLeft()], "($(alphabet[2botid-1]))")
			ax2 = Axis(grids[botid,2][1,1], xlabel = "Count", ylabel = L"e")
			vz = ṙp1[3,:]
			pos_pks_raw, _ = findmaxima(vz)
			pos_pks, _ = peakproms(pos_pks_raw,vz,minprom=0.2)
			neg_pks_raw, _ = findminima(vz)
			neg_pks, _ = peakproms(neg_pks_raw,vz,minprom=0.2)
			if (botid == length(bots))
				pos_pks = neg_pks.+1
			elseif length(neg_pks) > length(pos_pks)
				pop!(neg_pks)
			end
			ratios = -vz[pos_pks]./vz[neg_pks]
			scatter!(ax2,ratios)
			ylims!(ax2,ymids[botid]-0.08,ymids[botid]+0.08)
			# ax2.yticks = [0,0.3,0.7,1.0]
			ax2.xticks = collect(1:xtickmaxs[botid])
			Label(grids[botid,2][1,1,TopLeft()], "($(alphabet[2botid]))")
		end
		rowgap!(fig.layout,0)
		savefig(fig,figname)
		fig
	end
end

GM.activate!(); plotsave_pointmass_velocity_restitution(reshape(pms_hp,3,2)[:,1])
CM.activate!(); plotsave_pointmass_velocity_restitution(reshape(pms_hp,3,2)[:,1],"pointmass_v00_velocity_restitution")
GM.activate!(); plotsave_pointmass_velocity_restitution(reshape(pms_hp,3,2)[:,2])
CM.activate!(); plotsave_pointmass_velocity_restitution(reshape(pms_hp,3,2)[:,2],"pointmass_v10_velocity_restitution")

θ = 15 |> deg2rad
inclined_plane = TR.Plane([-tan(θ),0,1],zeros(3))

pm = new_pointmass(;e=0.0, μ=0.3, ro = [0.0,0,0], ṙo = [2.0cos(θ),0,2.0sin(θ)])

tspan = (0.0,0.6)
prob = TR.SimProblem(pm,(x)->pm_contact_dynfuncs(x;θ))
TR.solve!(prob,TR.ZhongCCP();tspan,dt=1e-4,ftol=1e-14,maxiters=50,exception=false)

GM.activate!(); plotsave_contact_persistent(pm)

# get_trajectory!(pm,1,1).u .|> norm

function plotsave_dis_vel(bot,figname=nothing)
	(;t) = bot.traj
	with_theme(mv_theme;
			resolution = (0.9tw,0.3tw)
		) do
		fig = Figure()
		ax1 = Axis(fig[1,1], xlabel = tlabel, ylabel = "Displacement (m)")
		Label(fig[1,1,TopLeft()], "($(alphabet[1]))")
		ax2 = Axis(fig[1,2], xlabel = tlabel, ylabel = "Velocity (m/s)")
		Label(fig[1,2,TopLeft()], "($(alphabet[2]))")
		rp1 = get_trajectory!(bot,1,1)
		ṙp1 = get_velocity!(bot,1,1)
		dp1 = rp1.u .|> norm
		vl1 = ṙp1.u .|> norm
		steps = 1:6001
		at = 0.3716
		stopstep = time2step(at,t)
		lines!(ax1,t[steps],dp1[steps])
		lines!(ax2,t[steps],vl1[steps])
		xlims!(ax1,0,0.6)
		xlims!(ax2,0,0.6)
		vlines!(ax1,t[stopstep])
		vlines!(ax2,t[stopstep])
		text!(ax1,latexstring("t=$at"), textsize = fontsize,
					position = (at+0.02, 0.19),
					align = (:left, :center)
		)
		text!(ax2,latexstring("t=$at"), textsize = fontsize,
					position = (at+0.02, 1.0),
					align = (:left, :center)
		)
		@show extrema(dp1[steps]), extrema(vl1[steps])
		# lines!(ax1,t[steps],rp1[1,steps])
		# lines!(ax1,t[steps],rp1[3,steps])
		# lines!(ax2,t[steps],ṙp1[1,steps])
		# lines!(ax2,t[steps],ṙp1[3,steps])
		savefig(fig,figname)
		fig
	end
end

GM.activate!(); plotsave_dis_vel(pm)
CM.activate!(); plotsave_dis_vel(pm,"pointmass_dis_vel")

me = TR.mechanical_energy!(pm); me.E[500:end] |> lines

me = TR.mechanical_energy!(pm); @show me.E[end]

contacts_traj_voa = VectorOfArray(pm.contacts_traj)
for (i,c) in enumerate(contacts_traj_voa[1,:])
	check_Coulomb(i,c)
end

GM.activate!(); with_theme(my_theme;
		fontsize = 6 |> pt2px,
		Axis3 = (
			azimuth = 4.575530633326984,
			elevation = 0.16269908169872405,
			zlabelvisible = false,
			yticklabelsvisible = false,
		)
	) do
	plot_traj!(pm;
		doslide=false,
		showmesh=false,
		# showinfo=false,
		gridsize=(1,3),
		attimes=(0.0,0.15,0.372),
		xlims=(-8e-3,0.5),
		ylims=(-0.1,0.1),
		zlims=(-8e-3,0.2),
		showlabels=false,
		ground=inclined_plane,
		figsize=(0.9tw,0.2tw),
		figname="pointmass_sliding"
	)
end

# contact_friction
function plotsave_pointmass_energy_friction(bot,figname=nothing)
	with_theme(mv_theme;
		resolution = (0.9tw,0.3tw)
	) do
		(;t) = bot.traj
		fig = Figure()
		ax1 = Axis(fig[1,1], xlabel = tlabel, ylabel = "Energy (J)")
		Label(fig[1,1,TopLeft()], "($(alphabet[1]))")
		ax2 = Axis(fig[1,2], xlabel = tlabel, ylabel = L"\lambda~(\mathrm{N})")
		Label(fig[1,2,TopLeft()], "($(alphabet[2]))")
		me = TR.mechanical_energy!(bot)
		(;T,V,E) = me
		lines!(ax1,t,me.T, label="Kinetic")
		lines!(ax1,t,me.V, label="Potential")
		lines!(ax1,t,me.E, label="Total")
		xlims!(ax1,extrema(t)...)
		axislegend(ax1)

		c1_traj = VectorOfArray(bot.contacts_traj)[1,5:end]
		Λn = [2c.state.Λ[1]./c.μ for c in c1_traj]
		Λt = [2c.state.Λ[3] for c in c1_traj]
		th = t[5:end]
		lines!(ax2,th,Λn,label=L"\lambda_{n}")
		lines!(ax2,(th[2:2:end]+th[1:2:end-1])./2,(Λt[2:2:end]+Λt[1:2:end-1])./2,label=L"\lambda_{t}")
		# lines!(ax2, th, Λt ,label=L"\lambda_{t}")
		@show Λn[end], Λt[begin], Λt[end]
		axislegend(ax2,position=:rc)
		xlims!(ax2,extrema(t)...)
		savefig(fig,figname)
		fig
	end
end
GM.activate!(); plotsave_pointmass_energy_friction(pm)
CM.activate!(); plotsave_pointmass_energy_friction(pm,"pointmass_energy_friction")
