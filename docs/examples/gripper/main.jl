using LinearAlgebra
using SparseArrays
using StaticArrays
using StructArrays
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using LaTeXStrings
using RecursiveArrayTools
using Cthulhu
using XLSX, DataFrames
using Unitful
using GeometryBasics
using TypeSortedCollections
using EponymTuples
using Match
using HomotopyContinuation
using DynamicPolynomials
using NLsolve
using Random
using Evolutionary
using Printf
using Revise
import TensegrityRobots as TR
includet("../analysis.jl")
includet("../vis.jl")
includet("gripper_define.jl")


tggriper = new_gripper(); tggriper_plot = deepcopy(tggriper)
plot_traj!(tggriper_plot;sup!)

# _,λ = TR.check_static_equilibrium_output_multipliers(tggriper.tg)

N = build_halftggriper_nullspace(tggriper,TR.get_q̌(tggriper.tg))
A = TR.make_A(tggriper)

A(TR.get_q(tggriper.tg))*N

plot_traj!(tggriper)

function get_poly(bot_input)
	bot = deepcopy(bot_input)
    (;tg) = bot
    # (;ndof,nconstraints,connectivity) = bot.tg
    # (;cables) = tg.tensiles
    # (;nfull,nfree) = connectivity.indexed
    # ncables = length(cables)
    # nλ = nconstraints
	gue = TR.get_initial(tg)
    Φ = TR.make_Φ(tg,gue.q)
    A = TR.make_A(tg,gue.q)
    Q̌ = TR.make_Q̌(tg,gue.q)
    S = TR.make_S(tg,gue.q)
    Ǩm_Ǩg = TR.make_Ǩm_Ǩg(tg,gue.q)

	pv = TR.get_polyvar(tg)

    pnq̌ = 1.0pv.q̌ .+ 0.0
    pns = 1.0pv.s .+ 0.0
    pnλ = 1.0pv.λ .+ 0.0
    pnd = 1.0pv.d .+ 0.0
    pnc = 1.0pv.c .+ 0.0
    pnk = 1.0pv.k .+ 0.0
    pnμ = 1.0pv.μ .+ 0.0
    polyΦ = Φ(pnq̌,pnd,pnc)
    polyA = A(pnq̌,pnc)
    polyQ̌ = Q̌(pnq̌,pns,pnμ,pnk,pnc)
    polyS = S(pnq̌,pns,pnc)
    polyQ̌a = transpose(polyA)*pnλ
    polyǨa = reduce(hcat,differentiate.(-polyQ̌a,Ref(pv.q̌))) |> transpose
    polyǨm, polyǨg = Ǩm_Ǩg(pnq̌,pns,pnμ,pnk,pnc)
    polyǨ = polyǨm .+ polyǨg .+ polyǨa
	polyŇ = build_halftggriper_nullspace(bot,pnq̌)
	poly𝒦 = transpose(polyŇ)*polyǨ*polyŇ

	polyP = [
		- polyQ̌ .- transpose(polyA)*pnλ ;
		polyS;
		polyΦ;
		# poly𝒦*pnξ.-pnζ.*pnξ;
		# transpose(pnξ)*pnξ-1;
	]

	# Ǩ0 = TR.build_Ǩ(bot.tg,gue.λ)
	# Ǩx = map(polyǨ) do z
	# 		z(
	# 			pv.q̌=>gue.q̌,
	# 			pv.s=>gue.s,
	# 			pv.λ=>gue.λ,
	# 			pv.μ=>gue.μ,
	# 			pv.k=>gue.k,
	# 			pv.d=>gue.d,
	# 			pv.c=>gue.c
	# 		)
	# 	end
	# # @show Ǩ0
	# @show Ǩ0.- Ǩx |> norm

	# P0 = map(polyP) do z
	# 	z(
	# 		pvq̌=>q̌0,
	# 		pvs=>s0,
	# 		pvλ=>λ0,
	# 		# pvξ=>ξ0,
	# 		pvμ=>μ0,
	# 		pvk=>k0,
	# 	    pvd=>d0,
	# 		pvc=>c0,
	# 		# pv.ζ=>ζ0
	# 	)
	# end
	# @show P0[                 1:nfree] |> norm
	# @show P0[           nfree+1:nfree+ncables] |> norm
	# @show P0[   nfree+ncables+1:nfree+ncables+nλ] |> norm
	# @show P0[nfree+ncables+nλ+1:nfree+ncables+nλ+ndof]
	# @show P0[end]
	polyP,poly𝒦,gue,pv
end

function pinpoint(bot_input)
	polyP, poly𝒦, gue, pv = get_poly(bot_input)
	ň = length(pv.q̌)
	ns = length(pv.s)
	nλ = length(pv.λ)
    function make_bf()
        function inner_pp!(f,x)
            q̌x = @view x[        1:ň]
            sx = @view x[      ň+1:ň+ns]
            λx = @view x[   ň+ns+1:ň+ns+nλ]
			Px = map(polyP) do z
                z(
                    pv.q̌=>q̌x,
                    pv.s=>sx,
                    pv.λ=>λx,
                    pv.μ=>gue.μ,
                    pv.k=>gue.k,
                    pv.d=>gue.d,
                    pv.c=>gue.c,
                )
			end

			f .= Px
        end
    end
    f_holder = zeros(ň+ns+nλ)
    x_initial = vcat(gue.q̌,gue.s,gue.λ)
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
    ini = @eponymtuple(
			q̌,s,λ,
			isconverged=converged(pp),
			d=gue.d, c=gue.c, μ=gue.μ, k=gue.k
	)
	polyP, ini, pv
end

function path_follow(bot_input)
	polyP, ini, pv = pinpoint(bot_input)
	variable_groups = [pv.q̌,pv.s,pv.λ]
	parameters = [pv.d;pv.c;pv.k;pv.μ]
	startsols = [[ini.q̌;ini.s;ini.λ]]
	start_parameters = [ini.d;ini.c;ini.k;ini.μ]
	target_parameters = [ini.d;ini.c;ini.k;ini.μ.+1.0]
	Psys = System(polyP;parameters)
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
		error("Tracking failed.")
	end
	path_result1 = path_results[1]
	sol = real(solution(path_result1))
	q̌,s,λ = TR.split_by_lengths(sol,length.(variable_groups))
	@eponymtuple(q̌,s,λ)
end

function pinpoint_critical(bot_input)
	bot = deepcopy(bot_input)
	(;tg) = bot
	polyP, poly𝒦, gue, pv_equ = get_poly(bot)

	ň = length(pv_equ.q̌)
	ns = length(pv_equ.s)
	nλ = length(pv_equ.λ)
	ndof = ň - nλ
	@polyvar ξ[1:ndof]
	@polyvar ζ
	pv = merge(pv_equ,@eponymtuple(ξ,ζ))
    pnξ = 1.0ξ .+ 0.0
	pnζ = 1.0ζ + 0.0

	F̌ = TR.build_F̌(tg,2,1,[-1.0,0.0])
	polyP[1:ň] .-= F̌*[pnζ]
	append!(polyP,poly𝒦*pnξ)
	push!(polyP,transpose(pnξ)*pnξ-1)
	# initial

	Ǩ0 = TR.build_Ǩ(bot.tg,gue.λ)
	# @show Ǩ0
	Ň0 = build_halftggriper_nullspace(bot,gue.q̌)
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


function path_follow_critical(bot_input)
	polyP, ini, pv = pinpoint_critical(bot_input)
	variable_groups = [pv.q̌,pv.s,pv.λ,pv.ξ,[pv.ζ]]
	parameters = [pv.d;pv.c;pv.k;pv.μ]
	startsols = [[ini.q̌;ini.s;ini.λ;ini.ξ;ini.ζ]]
	start_parameters = [ini.d;ini.c;ini.k;ini.μ]
	target_parameters = [ini.d;ini.c;ini.k;ini.μ.+1.0]
	Psys = System(polyP;parameters)
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
		error("Tracking failed.")
	end
	path_result1 = path_results[1]
	sol = real(solution(path_result1))
	q̌,s,λ,ξ,ζ = TR.split_by_lengths(sol,length.(variable_groups))
	@eponymtuple(q̌,s,λ,ξ,ζ)
end


x1 = 40.0; ϕ1 = π/3; ϕ2 = -π/5 # closed
x1 = 70.0; ϕ1 = π/2; ϕ2 = -π/10 # mid
x1 = 110.0; ϕ1 = (2.5)*π/3; ϕ2 = -π/10 # opened
tggriper = new_gripper(;x1,ϕ1,ϕ2); tggriper_plot = deepcopy(tggriper)
plot_traj!(tggriper_plot;xlims=[-50,250],ylims=[-100,100],sup!)

# get_poly(tggriper)

_,pp,_  = pinpoint(tggriper)
pp.isconverged

pp = path_follow(tggriper)


_,pp,_ = pinpoint_critical(tggriper)
pp.isconverged
pp.ζ

pp = path_follow_critical(tggriper)
pp.ζ

pprc = TR.recover(pp,tggriper.tg)
TR.set_new_initial!(tggriper_plot,pprc.q)
plot_traj!(tggriper_plot;xlims=[-50,250],ylims=[-100,100],sup!)

_,λ = TR.check_static_equilibrium_output_multipliers(tggriper_plot.tg)
Ǩ = TR.build_Ǩ(tggriper_plot.tg,pprc.λ)
Ň = nullspace(TR.make_A(tggriper_plot)(pprc.q))
𝒦 = transpose(Ň)*Ǩ*Ň
eg𝒦 = eigen(𝒦)

TR.get_q̌(tggriper.tg) .- pp.q̌
