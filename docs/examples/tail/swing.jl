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
using TensegrityRobot
const TR = TensegrityRobot
include("tail_define.jl")
include("plotting.jl")

include("../analysis.jl")

n = 13
tail = make_curve_tail(n)
# Vy = make_curve_tail(n)
_,_,tailscene = plotstructure(tail)
tailscene
# @code_warntype make_tail(n)
q0,q̇0,λ0 = TR.get_initial(tail)
# q̇0[end-3:end] .= 5*[0.1,0.0,0.1,0.0]

# TR.get_nbodyconstraint(tail)
# TR.get_nbodydof(tail) # TR.get_nbodycoords(tail)


function dynfuncs(tgstruct,q0)

    M = TR.build_massmatrix(tgstruct)
    Φ = TR.build_Φ(tgstruct,q0)
    A = TR.build_A(tgstruct)

    #Q̃=TR.build_Q̃(tgstruct)

    function F!(F,q,q̇,t)
        # du = 0.01*sin(t)
        # TR.actuate!(tgstruct.actuators[1],du)
        TR.reset_forces!(tgstruct)
        TR.distribute_q_to_rbs!(tgstruct,q,q̇)
        TR.update_cables_apply_forces!(tgstruct)
        # TR.apply_gravity!(tgstruct,factor=0.001)
        # F .= Q̃*TR.fvector(tgstruct)
        F .= 0.0
        TR.assemble_forces!(F,tgstruct)
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
kes,epes,gpes,es,es_err = analyse_energy(tail,sol;gravity=false,elasticity=true)
es
plot(es)
