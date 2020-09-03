using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using Parameters
import PyPlot; const plt = PyPlot
plt.pygui(true)
#using DifferentialEquations
#using ForwardDiff
# using DynamicPolynomials
# using HomotopyContinuation
using Printf
using CoordinateTransformations
using Makie
AbstractPlotting.__init__()
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

cd("examples\\spine3d")
include("define.jl")
include("plotting.jl")
include("../analysis.jl")

spine = spine3d(4)
_,_,_,spine_scene = plotstructure(spine)

q0,q̇0,λ0 = TR.get_initial(spine)

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
        TR.update_strings_apply_forces!(tgstruct)
        # TR.apply_gravity!(tgstruct,factor=0.001)
        # F .= Q̃*TR.fvector(tgstruct)
        F .= 0.0
        TR.assemble_forces!(F,tgstruct)
    end

    M,Φ,A,F!,nothing
end

M,Φ,A,F!,Jacs = dynfuncs(spine,q0)

dt = 0.01
prob = TS.DyProblem(dynfuncs(spine,q0),q0,q̇0,λ0,(0.0,5.0))
sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-14,verbose=true)
plotstructure(spine,sol,recordplot)
