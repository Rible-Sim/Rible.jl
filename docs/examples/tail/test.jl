using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using Makie
#plot(rand(10))
# using BenchmarkTools
# import PyPlot; const plt = PyPlot
# using LaTeXStrings
# using NLsolve
using Revise
using TensegrityRobots; const TR = TensegrityRobots
#cd("examples/tail")
includet("tail_define.jl")
includet("make_new_tail.jl")
includet("plotting.jl")
# includet("../analysis.jl")

n = 4
tail = make_new_tail(n)
#plotstructure(tail)
#@code_warntype make_tail(n)
q0,q̇0,λ0 = TR.get_initial(tail.tg)
q0=q =Vec{28,Float64}(0,0,2,0,0,-2,2,-2,0,-4,2,-4,0,-6,2,-6,0,-8,2,-8,0,-10,0,-12,0,-14,0,-16)
q̇0[end-3:end] .= 5*[0.1,0.0,0.1,0.0]
TR.set_new_initial!(tail,q0,q̇0)

function dynfuncs(tr)
    @unpack tg = tr
    M = TR.build_massmatrix(tg)
    Φ = TR.build_Φ(tg)
    A = TR.build_A(tg)

    #Q̃=TR.build_Q̃(tg)

    function F!(F,q,q̇,t)
        TR.reset_forces!(tg)
        TR.distribute_q_to_rbs!(tg,q,q̇)
        TR.update_cables_apply_forces!(tg)
        # F .= Q̃*TR.fvector(tgstruct)
        F .= 0.0
        TR.assemble_forces!(F,tg)
    end

    M,Φ,A,F!,nothing
end

dt = 0.01
prob = TR.DyProblem(dynfuncs(tail),tail,(0.0,20.0))
TR.solve!(tail,prob,TR.Zhong06(),dt=dt,ftol=1e-14,verbose=true)
# sol = TR.solve(prob,TR.ConstSPARK(1),dt=dt,ftol=1e-12,verbose=true)
#@code_warntype TR.solve(prob,TR.Zhong06(),dt=dt,ftol=1e-14,verbose=true)
plot_traj!(tail)
