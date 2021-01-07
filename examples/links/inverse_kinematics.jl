using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using Parameters
using Makie
AbstractPlotting.__init__()
plot(rand(10))
import PyPlot; const plt = PyPlot
plt.pygui(true)
using LaTeXStrings
#using DifferentialEquations
#using ForwardDiff
# using DynamicPolynomials
# using Printf
using CoordinateTransformations
# using HomotopyContinuation
using NLsolve
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot
cd("examples/links")
include("links_define.jl")
include("link_u_plotting.jl")
include("link_inverse.jl")


function gen_angles(α1 = 0.0, α2 = π/12, n = 6)
    α = α2 - α1
    if n == 1
        Δα = zero(α)
    else
        Δα = α/(n-1)
    end
    [α1 + (i-1)*Δα for i = 1:n]
end
θs1 = gen_angles(0.0, π/12, 1)
θs2 = gen_angles(0.0, π/12, 2)
θs32 = gen_angles(0.0, π/12, 32)


inverse_analysis(θs)


rl_indexes = [
[[1,2,3],[4,5,6]],
[[1],[2,3],[4],[5,6]]
]
linkn = inverse2energy_plot(θs2,rl_indexes)

rl_index = [[1],[2,3],[4],[5,6]]
inverse_with_energy_plot(θs32,rl_index)

rl_indexes = [rl_index,[g.+6 for g in rl_index],
            [g.+12 for g in rl_index]]
inverse_with_energy_3_plot(θs32,rl_indexes)

fig = position_plot(θs2;row=1,elev=11,azim=-76)
fig.savefig("link2_inv_pose.png",dpi=300,bbox_inches="tight")

fig = position_plot(θs2;nstage=4,row=1,zup=0.32,elev=11,azim=-76,y=-0.12)
fig.savefig("link4_inv_pose.png",dpi=300,bbox_inches="tight")
