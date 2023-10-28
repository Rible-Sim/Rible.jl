using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
# using BenchmarkTools
# import PyPlot; const plt = PyPlot
# using LaTeXStrings
using Makie
AbstractPlotting.__init__()

using NLsolve
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using Robot
const TR = Robot
include("tail_define.jl")
include("plotting.jl")

include("../analysis.jl")

n = 13
tail = make_curve_tail(n)
# Vy = make_curve_tail(n)
_,_,tailscene = plotstructure(tail)
tailscene
# @code_warntype make_tail(n)
q0,q̇0,λ0 = RB.get_initial(tail)
# q̇0[end-3:end] .= 5*[0.1,0.0,0.1,0.0]

# RB.get_nbodyconstraint(tail)
# RB.get_nbodydof(tail) # RB.get_nbodycoords(tail)


function dynfuncs(tgstruct,q0)

    M = RB.build_massmatrix(tgstruct)
    Φ = RB.build_Φ(tgstruct,q0)
    A = RB.build_A(tgstruct)

    #Q̃=RB.build_Q̃(tgstruct)

    function F!(F,q,q̇,t)
        # du = 0.01*sin(t)
        # RB.actuate!(tgstruct.actuators[1],du)
        RB.reset_forces!(tgstruct)
        RB.distribute_q_to_rbs!(tgstruct,q,q̇)
        RB.update_cables_apply_forces!(tgstruct)
        # RB.apply_gravity!(tgstruct,factor=0.001)
        # F .= Q̃*RB.fvector(tgstruct)
        F .= 0.0
        RB.assemble_forces!(F,tgstruct)
    end

    M,Φ,A,F!,nothing
end

M,Φ,A,F!,Jacs = dynfuncs(tail,q0)
# Φ(q0)
# @code_warntype Φ(q0)
# A(q0)
# @code_warntype A(q0)
F = similar(q0)
# @code_warntype F!(F,q0,q̇0,0.0)
dt = 0.01
prob = TS.DyProblem(dynfuncs(tail,q0),q0,q̇0,λ0,(0.0,10.0))
sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-14,verbose=true)

plotstructure(tail,sol,sliderplot)
plotstructure(tail,sol,recordplot)
# sol = TS.solve(prob,TS.ConstSPARK(1),dt=dt,ftol=1e-12,verbose=true)
@code_warntype TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-14,verbose=true)
kes,epes,gpes,restitution_coefficients,es_err = analyse_energy(tail,sol;gravity=false,elasticity=true)
restitution_coefficients
plot(restitution_coefficients)
