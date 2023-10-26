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
using Rible; const TR = Rible
#cd("examples/tail")
includet("tail_define.jl")
includet("make_new_tail.jl")
includet("plotting.jl")
includet("../analysis.jl")

n = 4
tail = make_new_tail(n)
#plotstructure(tail)
#@code_warntype make_tail(n)
q0,q̇0,λ0 = RB.get_initial(tail.st)
q0=q =Vec{28,Float64}(0,0,2,0,0,-2,2,-2,0,-4,2,-4,0,-6,2,-6,0,-8,2,-8,0,-10,0,-12,0,-14,0,-16)
q̇0[end-3:end] .= 5*[0.1,0.0,0.1,0.0]
RB.set_new_initial!(tail,q0,q̇0)

function dynfuncs(tr)
    @unpack st = tr
    M = RB.build_massmatrix(st)
    Φ = RB.build_Φ(st)
    A = RB.build_A(st)

    #Q̃=RB.build_Q̃(st)

    function F!(F,q,q̇,t)
        RB.reset_forces!(st)
        RB.distribute_q_to_rbs!(st,q,q̇)
        RB.update_cables_apply_forces!(st)
        # F .= Q̃*RB.fvector(tgstruct)
        F .= 0.0
        RB.assemble_forces!(F,st)
    end

    M,Φ,A,F!,nothing
end

dt = 0.01
prob = RB.DyProblem(dynfuncs(tail),tail,(0.0,20.0))
RB.solve!(tail,prob,RB.Zhong06(),dt=dt,ftol=1e-14,verbose=true)
# sol = RB.solve(prob,RB.ConstSPARK(1),dt=dt,ftol=1e-12,verbose=true)
#@code_warntype RB.solve(prob,RB.Zhong06(),dt=dt,ftol=1e-14,verbose=true)
plotstructure(tail)
