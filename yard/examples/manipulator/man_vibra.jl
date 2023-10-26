using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using BenchmarkTools
import PyPlot; const plt = PyPlot
using LaTeXStrings
# using NLsolve
using Revise
# using TS
# using Robot2D
# const TR = Robot2D
using TensegritySolvers; const TS = TensegritySolvers
using Robot
const TR = Robot

include("man_define.jl")

# ------------------Create Struture --------------------------
ndof = 6
refman = man_ndof(ndof,-π/12) # reference
manipulator = man_ndof(ndof,0.0)
# ------------------Create Struture\\-------------------------

q0,q̇0,λ0 = RB.get_initial(manipulator) #backup

dt = 0.01 # Same dt used for PID AND Dynamics solver

# ----------------------------Dynamics-----------------------------------
function dynfuncs(tgstruct,q0)

    M = RB.build_massmatrix(tgstruct)
    #M!, ∂T∂q̇! = TS.const_massmatrix_functions(M)
    Φ = RB.build_Φ(tgstruct,q0)
    A = RB.build_A(tgstruct)

    #Q̃=RB.build_Q̃(tgstruct)

    function F!(F,q,q̇,t)
        RB.reset_forces!(tgstruct)
        # RB.distribute_q_to_rbs!(tgstruct,q,q̇)
        # RB.update_cables_apply_forces!(tgstruct)
        RB.apply_gravity!(tgstruct)
        #F .= Q̃*RB.fvector(tgstruct)
        F .= 0.0
        RB.assemble_forces!(F,tgstruct)
        #@show isapprox(F,)
    end

    M,Φ,A,F!,nothing

    #A,Φ,∂T∂q̇!,F!,M!,nothing
end
M,Φ,A,F!,Jacs = dynfuncs(manipulator,q0)
Φ(q0)
A(q0)
prob = TS.DyProblem(dynfuncs(manipulator,q0),q0,q̇0,λ0,(0.0,100.0))
# RB.actuate!(manipulator,u)
# sol = TS.solve(prob,dt=dt,ftol=1e-13,verbose=true)

sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-12,verbose=true)

kes = [RB.kinetic_energy_coords(manipulator,q,q̇) for (q,q̇) in zip(sol.qs,sol.q̇s)]
pes = [RB.gravity_potential_energy(manipulator,q) for q in sol.qs]
es = kes .+ pes
plt.plot(sol.ts,es)
plt.ylim(-0.44,-0.46)
