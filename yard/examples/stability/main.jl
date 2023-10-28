using LinearAlgebra
using SparseArrays
using StaticArrays
using StructArrays
using Rotations
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
import Meshes
using Meshing
using LaTeXStrings
using CubicHermiteSpline
using RecursiveArrayTools
using Cthulhu
using XLSX, DataFrames
using Unitful
using GeometryBasics
using CoordinateTransformations
using TypeSortedCollections
using Interpolations
using EponymTuples
using Match
using HomotopyContinuation
using DynamicPolynomials
using NLsolve
using Random
using Evolutionary
using Printf
using Revise
import Rible as RB
include("../analysis.jl")
include("../vis.jl")
include("gripper_define.jl")
include("../test_prescribed/def3d.jl")
includet("../analysis.jl")
includet("../vis.jl")
includet("gripper_define.jl")
includet("../test_prescribed/def3d.jl")

figdir::String = raw"C:\Users\luo22\OneDrive\Papers\Ph.D.Thesis\st"

function get_μs_σs(bot,B,F̃)
	(;st) = bot
	ℓ = RB.get_cables_len(st)
	x,nb = RB.get_solution_set(B,F̃)
	σmax =  (ℓ - x)./ nb |> maximum
	σmin = -x./nb |> maximum
	σs = LinRange(σmin,σmax,21)
	friction_coefficients = [x + σ*nb for σ in σs]
	x,nb,friction_coefficients,σs
end

function compute_evs(bot,friction_coefficients,Ň)
	(;st) = bot
	evs = [	
		begin
			RB.set_restlen!(st,μ)
			RB.update!(st)
			ini = RB.get_initial(st)			
			# A = RB.make_A(st,ini.q)
			# Ň(ini.q̌)
			# A(ini.q̌)*Ň(ini.q̌)
			rb2 = RB.get_bodies(st)[2]
			ka = [
				transpose(rb2.state.loci_states[i]-rb2.state.ro)*rb2.state.fps[i]
				for i in eachindex(rb2.prop.loci)
			] |> sum
			Ň0 = Ň(ini.q̌,ini.c)
			# Ň0 = Ň(ini.q̌)
			Ǩm_Ǩg = RB.make_Ǩm_Ǩg(st,ini.q)
			Ǩm0, Ǩg0 = Ǩm_Ǩg(ini.q̌,ini.s,ini.μ,ini.k,ini.c)
			Ǩa0 = - RB.∂Aᵀλ∂q̌(st,ini.λ)
			𝒦m0 = transpose(Ň0)*Ǩm0*Ň0
			𝒦g0 = transpose(Ň0)*Ǩg0*Ň0
			𝒦a0 = transpose(Ň0)*Ǩa0*Ň0
			𝒦0 = 𝒦m0 .+ 𝒦g0 .+ 𝒦a0
			er = eigen(𝒦0)
			# v1 = er.vectors[:,1]
			# v1 ./= norm(v1)
			# v1 ./= sign(v1[argmax(abs.(v1))])
			# v1 = [0,1,0]		
			rng = MersenneTwister(0)
			v1 = rand(rng,length(er.vectors[:,1])).-0.5
			@show ka,v1	
			
			δVm = transpose(v1)*𝒦m0*v1
			δVg = transpose(v1)*𝒦g0*v1
			δVa = transpose(v1)*𝒦a0*v1
			δV = δVm+δVg+δVa
			@eponymtuple( values = er.values, vectors = er.vectors, ka, δVm,δVg,δVa,δV)
		end
		for μ in friction_coefficients[begin:end-1]
	] |> StructArray
end

function plotsave_evs(evs,σs,figname=nothing;haska=false)
	values = evs.values |> VectorOfArray
	(;δVm,δVg,δVa,δV,ka) = evs
	with_theme(theme_pub;
		resolution = ifelse(haska,(1.0tw,0.25tw),(0.9tw,0.25tw)),
		figure_padding = (fontsize,fontsize,0,0),
		markersize = 0.7fontsize,
		palette = (
            marker = tenmarkers,
            color = cgrad(:seaborn_bright, 10, categorical=true).colors.colors,
        ),
        Scatter = (
            cycle = Cycle(
                [:color, :marker], 
                covary = true
            ),
        ),
        Lines = (cycle = [:color],),
	) do 
		fig = Figure()
		ax1 = Axis(fig[1,1], xlabel = L"\sigma", ylabel = "Eigen Values")
		Label(fig[1,1,TopLeft()],"($(alphabet[1]))")
		for i in 1:ifelse(length(values[1])>1,3,1)
			lines!(ax1,σs,values[i,:])
			scatter!(ax1,σs,values[i,:],label=latexstring("\\rho_{($i)} "))
		end
		@show values[1,:] |> extrema
		ax2 = Axis(fig[1,2], xlabel = L"\sigma", ylabel = L"\delta V_{(1)}")
		Label(fig[1,2,TopLeft()],"($(alphabet[2]))")		
		bmarker = [:diamond,:hexagon,:star8]
		bcolor = cgrad(:seaborn_bright, 10, categorical=true).colors.colors[5:end]
		scatterlines!(ax2,σs,δVm;marker=bmarker[1],color=bcolor[1],label=L"\delta V_M")
		scatterlines!(ax2,σs,δVg;marker=bmarker[2],color=bcolor[2],label=L"\delta V_G")
		scatterlines!(ax2,σs,δVa;marker=bmarker[3],color=bcolor[3],label=L"\delta V_A")
		
		if haska
			ax3 = Axis(fig[1,3], xlabel = L"\sigma", ylabel = L"𝒦_{\lambda}")
			Label(fig[1,3,TopLeft()],"($(alphabet[3]))")
			ylims!(ax3,-20,0)
			scatterlines!(ax3,σs,ka;marker=:dtriangle)
			axislegend(ax2;
				position=:rc,
				# margin = (fontsize,fontsize,0,2fontsize)
			)
			axislegend(ax1;
				position=:lc,
				margin = (0.5fontsize,0,0,2.5fontsize)
			)
		else
			axislegend(ax1;
				position=:lc
			)
			axislegend(ax2;
				position=:lc,
				margin = (fontsize,fontsize,0,2fontsize)
			)
		end
		savefig(fig,figname)
		fig
	end
end

cRG = cgrad(:RdYlGn, 2, categorical=true).colors.colors

# ULTRA Spine 2D
spine2d2 = make_spine(2)
GM.activate!(); with_theme(theme_pub;
		figure_padding = (0,fontsize,0,fontsize),
		palette = (
			color = cgrad(:RdYlGn, 2, categorical=true).colors.colors[1:1],
		),
		Linesegments = (
			cycle = [:color], 
		),
		Scatter = (
			cycle = [],
		)
	) do 
	plot_traj!(spine2d2;
		figsize=(0.35tw,0.3tw),
		xlims=(-0.13,0.3),
		ylims=(-0.2,0.2),
		doslide=false,
		showinfo=false,
		showinit=true,
		showlabels=false,
		# showpoints=false,
		showground=false,
		rigidcolor=cRG[1],
		rigidlabelcolor=:slategray,
		rowgap=0,
		titleformatfunc = (sgi,tt) -> "",
		sup! = sup_spine2d!,
		figname = "ultra2"
	)
end
B,F̃ = RB.build_inverse_statics_for_restlength(spine2d2.st,spine2d2.st)

μ0, nb, friction_coefficients, σs = get_μs_σs(spine2d2,B,F̃)
@show μ0[[1,3]], nb[[1,3]]
nb
μ = friction_coefficients[end]
@show σs |> extrema

Ň = (q̌,c)-> RB.make_N(spine2d2.st,RB.get_q(spine2d2.st))(q̌)

evs = compute_evs(spine2d2,friction_coefficients,Ň)

evs.ka |> scatterlines
@show evs.δVa |> extrema

@show VectorOfArray(evs.values)[1,:] |> extrema

GM.activate!(); plotsave_evs(evs,σs[begin:end-1];haska=true)
CM.activate!(); plotsave_evs(evs,σs[begin:end-1],haska=true,"ultra2_eigen_values")


RB.check_stability!(spine2d2,Ň;
	scalings=[0.5,0.02,0.02]
)


GM.activate!(); with_theme(theme_pub;
		# fontsize = 6 |> pt2px,
		figure_padding = (0,fontsize,0,0),
		palette = (
			color = cgrad(:RdYlGn, 2, categorical=true).colors.colors,
		),
		Linesegments = (
			cycle = [:color], 
		),
		Scatter = (
			cycle = [],
		)
	) do 
		plot_traj!(spine2d2;
		figsize=(0.9tw,0.3tw),
		xlims=(-0.13,0.3),
		ylims=(-0.2,0.2),
		gridsize=(1,3),
		atsteps=collect(2:4),
		doslide=false,
		showinfo=false,
		showinit=true,
		showlabels=false,
		# showpoints=false,
		showground=false,
		# rigidcolor=cRG[2],
		rigidlabelcolor=:slategray,
		rowgap=0,
		titleformatfunc = (sgi,tt) -> "($(alphabet[sgi])) Mode $sgi",
		sup! = sup_spine2d!,
		figname = "ultra2_modes"
	)
end

# ULTRA Spine 3D

spine3d2 = spine3d(2)
plot_traj!(spine3d2;showground=false)

B,F̃ = RB.build_inverse_statics_for_restlength(spine3d2.st,spine3d2.st)
b = zeros(eltype(B),1,size(B,2))
ℓ = RB.get_cables_len(spine3d2.st)
κ = RB.get_cables_stiffness(spine3d2.st)
b[4] = -κ[4]
b[1] =   κ[1]
B⁺ = vcat(B,b)
display(B)
display(B⁺)
push!(F̃,κ[1]*ℓ[1] - κ[4]*ℓ[4] )
display(F̃)
μ0, nb, friction_coefficients,σs = get_μs_σs(spine3d2,B⁺,F̃)
display(μ0)
display(nb)
@show size(nb)
μ = (μ0 .+ σs[end].*nb[:])
f = κ .*( ℓ .- μ )
f[1]/f[4]

@show σs |> extrema

Ň = (q̌,c) -> RB.make_N(spine3d2.st,RB.get_q(spine3d2.st))(q̌)

evs = compute_evs(spine3d2,friction_coefficients,Ň)
@show evs.δVa |> extrema

GM.activate!(); plotsave_evs(evs,σs[begin:end-1])
CM.activate!(); plotsave_evs(evs,σs[begin:end-1],"ultra_eigen_values")
# μ = RB.inverse_for_restlength(spine3d2,spine3d2;eps_abs=1e-8,eps_rel=1e-8)

RB.set_restlen!(spine3d2.st,friction_coefficients[1])
RB.update!(spine3d2.st)
_,λ = RB.check_static_equilibrium_output_multipliers(spine3d2.st)

RB.check_stability!(spine3d2,Ň;scalings=[0.5,0.5,0.5,0.01,0.01,0.01])

spine3d2.traj.t
GM.activate!(); with_theme(theme_pub;
		fontsize = 6 |> pt2px,
		figure_padding = (0,0,0,fontsize)
	) do 
	plot_traj!(spine3d2;
		AxisType=Axis3,
		figsize=(0.9tw,0.5tw),
		xlims=(-0.03,0.03),
		ylims=(-0.03,0.03),
		zlims=(-0.03,0.07),
		gridsize=(2,3),
		atsteps=collect(2:7),
		doslide=false,
		showinfo=false,
		showinit=true,
		showlabels=false,
		showground=false,
		rigidcolor=cRG[2],
		rowgap=0,
		titleformatfunc = (sgi,tt) -> "($(alphabet[sgi])) Mode $sgi",
		sup! = (ax,_,_) -> begin
			ax.azimuth = 3.595530633326982
			ax.elevation = 0.26979632679489624
			ax.zlabeloffset = 2fontsize
			hidex(ax)
			hidey(ax)
		end,
		figname = "ultra_modes"
	)
end

# Manipulator
function make_Ň(st,q)
	(;ndof,connectivity) = st
	(;indexed,numbered) = connectivity
	(;nfree,mem2sysfree) = indexed
	(;mem2num,num2sys) = numbered
	N = RB.make_N(st,q)
	function inner_Ň(q̌,c)
		# ret = zeros(eltype(q̌),nfree,ndof)
		bodyid = 2
		pid = 1
		q_rb = @view q̌[mem2sysfree[bodyid]]
		u_rb = SVector{2}(q_rb[3:4])
		v_rb = SVector{2}(q_rb[5:6])
		ic = mem2num[bodyid][pid]
		c1,c2 = c[num2sys[ic]]
		ret = N(q̌)*[
			c1.*skew(u_rb) .+ c2.*skew(v_rb);
			1
		]
		ret
	end
end

man1 = dualtri(1;θ=deg2rad(0))
man1.st.ndof
# plot_traj!(man1;)  
q0 = RB.get_q(man1.st)

N = RB.make_N(man1.st,q0)
N(RB.get_q̌(man1.st))


Ň = make_Ň(man1.st,q0)
Ň(RB.get_q̌(man1.st),RB.get_c(man1.st))
A = RB.make_A(man1.st,q0)
A(RB.get_q̌(man1.st))*Ň(RB.get_q̌(man1.st),RB.get_c(man1.st)) |> norm

B,F̃ = RB.build_inverse_statics_for_restlength(man1.st,man1.st)

μ0, nb, friction_coefficients,σs = get_μs_σs(man1,B,F̃)
μ0
@show Δμ = nb[:,1]
@show σs |> extrema

evs = compute_evs(man1,friction_coefficients,Ň)
values = evs.values |> VectorOfArray

GM.activate!(); plotsave_evs(evs,σs[begin:end-1];)
CM.activate!(); plotsave_evs(evs,σs[begin:end-1],"man1_eigen_values")

function pinpoint_equilibriums(μ0,Δμ,
		θs = [-50, 0, 50] .|> deg2rad,
		σ = 0.19
	)
	bots = [
		begin
			bot = dualtri(1;θ)
			RB.set_restlen!(bot.st,μ0.+σ.*Δμ)
			# RB.update!(bot.st)
			polyP, poly𝒦, ini, pv = RB.pinpoint(bot;Ň)
			@assert ini.isconverged
			inirc = RB.recover(ini,bot.st)
			RB.set_new_initial!(bot,inirc.q)
			RB.update!(bot.st)
			_, _, er = RB.check_stability(bot.st)
			@show er.values
			# plot_traj!(dualtri_up;)
			bot
		end
		for θ in θs
	]
	push!(bots[1].traj,bots[2].st.state.system)
	push!(bots[1].traj,bots[3].st.state.system)
	bots[1]
end

man1_bi = pinpoint_equilibriums(μ0,Δμ)


function forward_bi(bot;
		change = :c,
	)
	botfor = deepcopy(bot)
	[
		begin
			RB.goto_step!(bot,i)
			ini = RB.get_initial(bot.st)
			startsol = (ini.q̌,ini.s,ini.λ)
			start_parameters = (d=ini.d,c=ini.c,k=ini.k,u=ini.μ,g=[0.0])
			(;mem2num,num2sys) = bot.st.connectivity.numbered
			target_parameters = deepcopy(start_parameters)
			if change == :c
				target_parameters.c[num2sys[mem2num[1][1]]] .+= [0.00,0.02]
				target_parameters.c[num2sys[mem2num[1][3]]] .+= [0.02,0.0]
				target_parameters.c[num2sys[mem2num[2][3]]] .+= [0.02,0.0]
				target_parameters.c[num2sys[mem2num[2][2]]] .+= [0.02,0.02]
			else
				target_parameters.u .-= 0.01
			end

			@show i
			seq = RB.forward_sequence(
				bot.st,
				startsol,
				start_parameters,
				target_parameters,
				RB.AllMode();n=10
			)
			[
				begin
					stptrc = RB.recover(stpt,botfor.st)
					(;q,q̌,s,λ,d,c,k,u,g) = stptrc
					Q̌ = RB.make_Q̌(botfor.st,q)(q̌,s,u,k,c)
					A = RB.make_A(botfor.st,q)(q̌,c)
					# @show Q̌+transpose(A)*λ |> norm
					# @show stptrc.c[num2sys[cidx]]
					gěneralized_forces = Q̌
					# @show gěneralized_forces
					RB.set_restlen!(botfor.st,u)
					botfor.st.state.system.c .= c
					botfor.st.state.system.q̌ .= q̌
					RB.update!(botfor.st)
					# RB.check_static_equilibrium_output_multipliers!(botfor.st,q;stpt=stptrc)
					# # RB.update_orientations!(botfor.st)
					# _,_,er = RB.check_stability(botfor.st,stptrc.λ)
					V = RB.potential_energy(botfor.st)

					_,_,er = RB.check_stability(botfor.st)
					k = er.values[1]
					θ = get_angles(botfor.st)[1]
					@eponymtuple(V,k,θ,state=deepcopy(botfor.st.state.system))
				end
				for stpt in seq
			] |> StructArray
		end
		for i = 1:3
	] |> VectorOfArray
end

change::Symbol = :c
change::Symbol = :u
Vkθs = forward_bi(man1_bi;change)

man1_bi_seqs = [
	begin
		bot = dualtri(1;)
		bot.traj[1] = Vkθ.state[1]
		append!(bot.traj,Vkθ.state[2:end])
		bot
	end
	for Vkθ in Vkθs
]

GM.activate!(); with_theme(theme_pub;
		# fontsize = 6 |> pt2px,
		Axis = (
			titlefont = "CMU Serif",
			titlesize = fontsize,
			titlegap = 0
		),
		figure_padding = (0,fontsize,0,0),
		Scatter = (
			cycle = [],
		)
	) do 
	fig = Figure(resolution=(1.0tw,0.28tw))
	colormap = cgrad(:RdYlGn, 3, categorical=true)
	for i  = 1:3
		gridi = fig[1,i] = GridLayout()
		title = "($(alphabet[i])) $(["稳定平衡","不稳定平衡","稳定平衡"][i])"
		plot_traj!(man1_bi_seqs[i];
			fig = gridi,
			xlims=(-0.05,0.3),
			ylims=(-0.15,0.15),
			# gridsize=(1,3),
			atsteps=collect(1:5:11),
			doslide=false,
			showinfo=false,
			# showinit=true,
			showlabels=false,
			# showpoints=false,
			showground=false,
			# rigidcolor=cRG[2],
			rigidlabelcolor=:slategray,
			rowgap=0,
			titleformatfunc = (sgi,tt) -> title,
			sup! = (ax,tgob,sgi) -> sup_dualtri!(ax,tgob,;
				color=colormap.colors[sgi],
			),
		)
	end
	len = 3
	if change == :c
		xmin,xmax = 0.00,0.02
		xs = range(xmin,xmax,length=len)
		ticks = (xs,string.(xs))
	else
		xmin,xmax = 0.18,0.19
		xs = range(xmin,xmax,length=len)
		ticks = (xs,string.(xs) |> reverse)
	end
	dx = (xmax-xmin)/(len-1)
	Colorbar(fig[1, 4];
		colormap=cgrad(:RdYlGn, 3, categorical=true),
		label=L"\sigma",
		limits=(xmin-dx/2, xmax+dx/2),
		ticks
	)
	savefig(fig,"dualtri_$change")
	fig
end


function plotsave_Vkθ(Vkθs,figname=nothing)
	colormap = cgrad(:brg, 11, categorical=true)
	if change == :c
		xs = range(0.00,0.02,length=11)
		dx = xs[2]-xs[1]
	else
		xs = range(0.19,0.18,length=11)
		dx = xs[1]-xs[2]
	end
	with_theme(theme_pub;
		figure_padding = (0,fontsize,0,fontsize),
		resolution = (0.5tw,0.3tw),	
		palette = (
			color = colormap,
		),
		Lines = (
			cycle = [:color], 
		),
	) do
		fig = Figure()
		ax = Axis(fig[1,1],xlabel = L"\theta~(\degree)", ylabel = L"V")
		for i = 1:3
			lines!(ax,Vkθs[i].θ,Vkθs[i].V;color=:black)
		end
		for j = 1:length(Vkθs[1])
			Vkθsj = Vkθs[j,:] |> StructArray
			x = Vkθsj.θ
			y = Vkθsj.V
			gradient = zero.(x)			
			@show x,y,gradient
			spl = CubicHermiteSplineInterpolation(x, y, gradient)
			xi = range(minimum(x), maximum(x), length=21)
			lines!(ax,xi,spl(xi))
		end
		Colorbar(
			fig[1,2];
			colormap,
			limits = (minimum(xs)-dx/2,maximum(xs)+dx/2),
			ticks = ifelse(change==:c,(xs,string.(xs)),(xs,string.(xs)|> reverse)),
			label = ifelse(change==:c,L"\Delta c",L"\sigma")
		)
		savefig(fig,figname)
		fig
	end
end

GM.activate!(); plotsave_Vkθ(Vkθs)
CM.activate!(); plotsave_Vkθ(Vkθs,"man1_Vkθ_$change")


# gripper
tggriper = new_gripper(); tggriper_plot = deepcopy(tggriper)
GM.activate!(); plot_traj!(tggriper_plot; sup! = sup_ggriper!)

# _,λ = RB.check_static_equilibrium_output_multipliers(tggriper.st)

Ň = make_halftggriper_nullspace(tggriper,RB.get_q̌(tggriper.st))
A = RB.make_A(tggriper)

A(RB.get_q(tggriper.st))*Ň(RB.get_q̌(tggriper.st),RB.get_c(tggriper.st)) |> norm

plot_traj!(tggriper)


function pinpoint_critical(bot_input;
		Ň,		
		F̌
	)
	bot = deepcopy(bot_input)
	(;st) = bot
	polyP, poly𝒦, gue, pv_equ = RB.get_poly(bot;Ň)

	ň = length(pv_equ.q̌)
	ns = length(pv_equ.s)
	nλ = length(pv_equ.λ)
	ndof = ň - nλ
	@polyvar ξ[1:ndof]
	@polyvar ζ
	pv = merge(pv_equ,@eponymtuple(ξ,ζ))
    pnξ = 1.0ξ .+ 0.0
	pnζ = 1.0ζ + 0.0

	polyP[1:ň] .-= F̌*[pnζ]
	append!(polyP,poly𝒦*pnξ)
	push!(polyP,transpose(pnξ)*pnξ-1)
	# initial

	Ǩ0 = RB.build_Ǩ(bot.st,gue.λ)
	# @show Ǩ0
	Ň0 = Ň(gue.q̌,gue.c)
	𝒦0 = transpose(Ň0)*Ǩ0*Ň0
	eg𝒦0 = eigen(𝒦0)
	# 𝐞 = [1:ndof .== i for i in 1:ndof]
	# i = 1; j = 2
	# 𝒦0ij = (I-𝐞[i]*transpose(𝐞[i]))*𝒦0+𝐞[i]*transpose(𝐞[j])
	# ξ0 = 𝒦0ij\𝐞[i]
	@show eg𝒦0.values
	ξ0 = eg𝒦0.vectors[:,1]
	ξ0 /= norm(ξ0)
	@show transpose(ξ0)*𝒦0*ξ0
	ζ0 = 0.0
	@show ξ0,ζ0
    function make_bf()
        function inner_pp!(f,x)
            q̌x = @view x[             1:ň]
            sx = @view x[           ň+1:ň+ns]
            λx = @view x[        ň+ns+1:ň+ns+nλ]
            ξx = @view x[     ň+ns+nλ+1:ň+ns+nλ+ndof]
            ζx =       x[end]
			Px = map(polyP) do z
                z(
                    pv.q̌=>q̌x,
                    pv.s=>sx,
                    pv.λ=>λx,
					pv.ξ=>ξx,
					pv.ζ=>ζx,
                    pv.μ=>gue.μ,
                    pv.k=>gue.k,
                    pv.d=>gue.d,
                    pv.c=>gue.c,
                )
			end

			f .= Px
        end
    end
    f_holder = zeros(ň+ns+nλ+ndof+1)
    x_initial = vcat(gue.q̌,gue.s,gue.λ,ξ0,ζ0)
    pp! = make_bf()
    pp = nlsolve(pp!,x_initial,ftol=1e-10,iterations=100,method=:newton)
    # @show
    pp!(f_holder,pp.zero)
	# @show f_holder |> norm
    # @show f_holder[                 1:ň+ns+nλ] |> norm
    # @show f_holder[ň+ns+nλ+1:ň+ns+nλ+ndof] |> norm
    # @show f_holder[end]
    # @show  pp.zero[ň+ns+nλ+1:ň+ns+nλ+ndof]
    # @show  pp.zero[end]
    q̌ = pp.zero[        1:ň]
    s = pp.zero[      ň+1:ň+ns]
    λ = pp.zero[   ň+ns+1:ň+ns+nλ]
    ξ = pp.zero[ň+ns+nλ+1:ň+ns+nλ+ndof]
    ζ = pp.zero[end]
	ini = @eponymtuple(
		q̌,s,λ,ξ,ζ,isconverged=converged(pp),
		d=gue.d, c=gue.c, μ=gue.μ, k=gue.k
	)
    polyP,ini,pv
end

function pinpoint_gripper(Ň)
	guesses = [
		(x1 = 110.0, ϕ1 = (2.5)*π/3, ϕ2 = -π/10), # opened
		(x1 = 70.0, ϕ1 = π/2, ϕ2 = -π/10), # mid
		(x1 = 40.0, ϕ1 = π/3, ϕ2 = -π/5), # closed
	]
	bots = [	
		begin
			bot = new_gripper(;guess...)
			# c = RB.get_c(bot.st)
			# # c[7] -= 0.1
			# RB.update_points!(bot.st,c)
			polyP, poly𝒦, ini, pv  = RB.pinpoint(bot;Ň)
			@assert ini.isconverged
			inirc = RB.recover(ini,bot.st)
			RB.set_new_initial!(bot,inirc.q)
			RB.update!(bot.st)
			_, _, er = RB.check_stability(bot.st)
			@show er.values
			# plot_traj!(dualtri_up;)
			bot
		end
		for guess in guesses
	]
	botcri = deepcopy(bots[1])
	F̌ = RB.build_F̌(botcri.st,2,1,[-1.0,0.0])
	polyPcri, inicri, pvcri  = pinpoint_critical(botcri;Ň,F̌)
	@assert inicri.isconverged
	inicrirc = RB.recover(inicri,botcri.st)
	RB.set_new_initial!(botcri,inicrirc.q)
	RB.update!(botcri.st)
	@show inicri.ζ
	_, _, er = RB.check_stability(botcri.st;F̌=inicri.ζ.*F̌)
	@show er.values
	push!(bots[1].traj,botcri.st.state.system)
	push!(bots[1].traj,bots[2].st.state.system)
	push!(bots[1].traj,bots[3].st.state.system)
	bots[1],polyPcri, inicri, pvcri
end

tggripers,polyP, ini, pv = pinpoint_gripper(Ň)

plot_traj!(tggripers;
	xlims=[-50,250],ylims=[-100,100],
	# atsteps=[1,2],
	sup! = sup_ggriper!
)

tggriper.st.connectivity.numbered.mem2num[2][2]
tggriper.st.connectivity.numbered.num2sys[4]

GM.activate!(); with_theme(theme_pub;
		fontsize = 6 |> pt2px,
		Axis = (
			titlefont = "CMU Serif",
			titlesize = fontsize,
			titlegap = 0
		),
		figure_padding = (0,fontsize,0,0),
		Scatter = (
			cycle = [],
		)
	) do 
	fig = Figure(resolution=(1.0tw,0.20tw))
	plot_traj!(tggripers;
		fig,
		xlims=(-50,200),
		ylims=(-75,75),
		gridsize=(1,4),
		atsteps=collect(1:4),
		doslide=false,
		showinfo=false,
		# showinit=true,
		showlabels=false,
		# showpoints=false,
		showground=false,
		# rigidcolor=cRG[2],
		rigidlabelcolor=:slategray,
		rowgap=0,
		titleformatfunc = (sgi,tt) -> "($(alphabet[sgi])) $(["稳定（张开）","临界稳定（最大推力）","不稳定","稳定（闭合）"][sgi])",
		sup! = (ax,tgob,sgi) -> sup_ggriper!(ax,tgob,sgi;
			# color=colormap.colors[sgi],
		),
	)
	savefig(fig,"gripper")
	fig
end

	
Δc7=-0.1
Δu6=-30.0


function forward_bi_gripper(bot;
		Δc7,
		Δu6
	)
	botfor = deepcopy(bot)
	[
		begin
			RB.goto_step!(bot,i)
			ini = RB.get_initial(bot.st)
			startsol = (ini.q̌,ini.s,ini.λ)
			start_parameters = (d=ini.d,c=ini.c,k=ini.k,u=ini.μ,g=[0.0])
			(;mem2num,num2sys) = bot.st.connectivity.numbered
			target_parameters = deepcopy(start_parameters)
			target_parameters.c[7] += Δc7
			target_parameters.u[6] += Δu6

			@show i
			seq = RB.forward_sequence(
				bot.st,
				startsol,
				start_parameters,
				target_parameters,
				RB.AllMode();n=10
			)
			[
				begin
					stptrc = RB.recover(stpt,botfor.st)
					(;q,q̌,s,λ,d,c,k,u,g) = stptrc
					Q̌ = RB.make_Q̌(botfor.st,q)(q̌,s,u,k,c)
					A = RB.make_A(botfor.st,q)(q̌,c)
					# @show Q̌+transpose(A)*λ |> norm
					# @show stptrc.c[num2sys[cidx]]
					gěneralized_forces = Q̌
					# @show gěneralized_forces
					RB.set_restlen!(botfor.st,u)
					botfor.st.state.system.c .= c
					@show c[7]
					botfor.st.state.system.q̌ .= q̌
					RB.update!(botfor.st)
					# RB.check_static_equilibrium_output_multipliers!(botfor.st,q;stpt=stptrc)
					# # RB.update_orientations!(botfor.st)
					# _,_,er = RB.check_stability(botfor.st,stptrc.λ)
					V = RB.potential_energy(botfor.st)

					_,_,er = RB.check_stability(botfor.st)
					k = er.values[1]
					θ = botfor.st.state.system.q̌[1]
					@eponymtuple(V,k,θ,state=deepcopy(botfor.st.state.system))
				end
				for stpt in seq
			] |> StructArray
		end
		for i = [1,3,4]
	] |> VectorOfArray
end

Vkθs = forward_bi_gripper(tggripers;Δc7,Δu6)


gripper_bi_seqs = [
	begin
		bot = new_gripper()
		bot.traj[1] = Vkθ.state[1]
		append!(bot.traj,Vkθ.state[2:end])
		bot
	end
	for Vkθ in Vkθs
]

plot_traj!(gripper_bi_seqs[1];
	xlims=[-50,250],ylims=[-100,100],
	# atsteps=[1,5,11],
	sup! = sup_ggriper!
)

GM.activate!(); with_theme(theme_pub;
		# fontsize = 6 |> pt2px,
		markersize = 10,
		Axis = (
			titlefont = "CMU Serif",
			titlesize = fontsize,
			titlegap = 0
		),
		figure_padding = (0,0,0,0),
		Scatter = (
			cycle = [],
		)
	) do 
	fig = Figure(resolution=(1.0tw,0.28tw))
	colormap = cgrad(:RdYlGn, 2, categorical=true)
	for i  = 1:3
		gridi = fig[1,i] = GridLayout()
		title = "($(alphabet[i])) $(["稳定（张开）","不稳定","稳定（闭合）"][i])"
		plot_traj!(gripper_bi_seqs[i];
			fig = gridi,
			xlims=[-10,180],
			ylims=[-75,75],
			# gridsize=(1,3),
			atsteps=[1,11],
			doslide=false,
			showinfo=false,
			# showinit=true,
			showlabels=false,
			# showpoints=false,
			showground=false,
			# rigidcolor=cRG[2],
			rigidlabelcolor=:slategray,
			rowgap=0,
			titleformatfunc = (sgi,tt) -> title,
			sup! = (ax,tgob,sgi) -> sup_ggriper!(ax,tgob,0;
				color=colormap.colors[sgi],
				linewidth = 6
			),
		)
	end
	len = 2
	xmin,xmax = 0.0,1.0
	xs = range(xmin,xmax,length=len)
	ticks = (xs,string.(xs))
	dx = (xmax-xmin)/(len-1)
	Colorbar(fig[1, 4];
		colormap=cgrad(:RdYlGn, len, categorical=true),
		label=L"\sigma",
		limits=(xmin-dx/2, xmax+dx/2),
		ticks
	)
	savefig(fig,"gripper_forward")
	fig
end

function path_follow_critical(polyP, ini, pv;
		n = 10,
		Δc7,
		Δu6
	)
	
	bot0 = new_gripper()
	bot = deepcopy(bot0)

	variable_groups = [pv.q̌,pv.s,pv.λ,pv.ξ,[pv.ζ]]
	parameters = [pv.d;pv.c;pv.k;pv.μ]
	startsols = [[ini.q̌;ini.s;ini.λ;ini.ξ;ini.ζ]]
	@show ini.ζ
	start_parameters = [ini.d;ini.c;ini.k;ini.μ]
	# target_parameters = [ini.d;ini.c;ini.k;ini.μ.-20.0]

	Psys = System(polyP;parameters)
	
	thepath = [
		begin 
			tar = deepcopy(ini)
			@show tar.c[7]
			tar.c[7] += Δc7*i/n
			@show tar.μ[6]
			tar.μ[6] += Δu6*i/n
			target_parameters = [tar.d;tar.c;tar.k;tar.μ]
			result = HomotopyContinuation.solve(
				Psys,
				startsols;
				start_parameters,
				target_parameters,
				threading = false
			)
			path_results = results(result)
			if length(path_results) != 1
				@show failed(result)
				@warn "Tracking failed."
			end
			path_result1 = path_results[1]
			sol = real(solution(path_result1))
			q̌,s,λ,ξ,ζ = RB.split_by_lengths(sol,length.(variable_groups))
			pp = @eponymtuple(q̌,s,λ,ξ,ζ=ζ[1])
			pprc = RB.recover(pp,bot0.st)
			RB.set_restlen!(bot.st,tar.μ)
			RB.update_points!(bot.st,tar.c)		
			RB.set_new_initial!(bot,pprc.q)
			if i == 0
				bot0.traj[1] = deepcopy(bot.st.state.system)
			else
				append!(bot0.traj,deepcopy(bot.traj))
			end
			@show pprc.ζ
			merge(pprc,(V=RB.mechanical_energy(bot.st).V,))
		end
		for i = 0:n
	] |> StructArray
	bot0, thepath
end


tggriper, cripath = path_follow_critical(polyP, ini, pv; Δc7, Δu6)
ζmin,ζmax = cripath.ζ |> extrema .|> abs 
@show ζmax,ζmin
@show (ζmax-ζmin)/ζmax

VectorOfArray(cripath.q̌)[1,:] 
cripath.V
CM.activate!(); with_theme(theme_pub;
		resolution = (1.0tw,0.42tw),
		figure_padding = (0,fontsize,0,0),
		palette = (
			color = cgrad(:RdYlGn, 11, categorical=true),
		),
		Lines = (
			cycle = [:color], 
		),
		Scatter = (
			cycle = [],
		)
	) do
	fig = Figure()
	colormap = cgrad(:RdYlGn, 3, categorical=true)
	grid1 = fig[1,1] = GridLayout()
	Label(grid1[1,1,TopLeft()],"($(alphabet[1]))")
	plot_traj!(tggriper;
		fig=grid1,
		xlims=[-10,180],
		ylims=[-75,75],
		atsteps=collect(1:5:11),
		doslide=false,
		showlabels=false,
		titleformatfunc = (sgi,tt) -> "",
		sup! = (ax,tgob,sgi) -> sup_ggriper!(ax,tgob,0;
			color=colormap.colors[sgi],
		),
	)
	grid2 = fig[1,2] = GridLayout()
	xmin, xmax = 0.0,1.0
	len = 10
	dx = 0.1
	ax1 = Axis(grid2[1,1],xlabel=L"\sigma",ylabel=L"\zeta")
	Label(grid2[1,1,TopLeft()],"($(alphabet[2]))")
	scatterlines!(ax1,xmin:dx:xmax,cripath.ζ;
		color=collect(xmin:dx:xmax),
		colormap=:RdYlGn,
		markercolor = :blue,
		linewidth = 4,
		markersize = 8
	)
	Colorbar(grid2[2, 2];
		colormap=cgrad(:RdYlGn, 11, categorical=true),
		label=L"\sigma",
		limits=(xmin, xmax),
		# ticks
	)
	ax2 = Axis(grid2[2,1],xlabel = L"x_{2,1}", ylabel = L"V")
	Label(grid2[2,1,TopLeft()],"($(alphabet[3]))")
	gθ = [
		begin 
			Vkθsj = Vkθs[j,:] |> StructArray
			x = Vkθsj.θ |> reverse
			y = Vkθsj.V |> reverse
			gradient = zero.(x)			
			# @show x,y,gradient
			spl = CubicHermiteSplineInterpolation(x, y, gradient)
			θ = range(minimum(x), maximum(x), length=21)
			lines!(ax2,θ,spl(θ))
			xi = range(minimum(x), maximum(x), length=1001)
			idx = spl(xi; grad=true) |> argmin
			(Vmin=spl(xi)[idx],θmin=xi[idx])
		end
		for j = 1:length(Vkθs[1])
	] |> StructArray
	for i = 1:3
		lines!(ax2,Vkθs[i].θ   ,Vkθs[i].V  ;color=:black)
	end
	lines!(ax2, gθ.θmin ,gθ.Vmin  ;color=:blue)
	rowsize!(grid2,1,Relative(0.3))
	rowgap!(grid2,0)
	savefig(fig,"gripper_cri")
	fig
end
