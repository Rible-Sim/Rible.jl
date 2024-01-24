using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using Parameters
using Makie
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
using Rible; const TR = Rible
cd("examples/links")
include("links_define.jl")
include("link_u_plotting.jl")
include("link_inverse.jl")
include("../plot_helpers.jl")
set_pyplot()
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


_,rls = inverse_analysis(θs2)
rls

rl_indexes = [
[[1,2,3],[4,5,6]],
[[1],[2,3],[4],[5,6]]
]
fig,_ = inverse2energy_plot(θs2,rl_indexes)
fig.savefig("link2_inv_pose_energy.png",dpi=300,bbox_inches="tight")

rl_index = [[1],[2,3],[4],[5,6]]
fig = inverse_with_energy_plot(θs32,rl_index)
fig.savefig("link2_inv_by_energy.png",dpi=300,bbox_inches="tight")

rl_indexes = [rl_index,[g.+6 for g in rl_index],
            [g.+12 for g in rl_index]]
fig = inverse_with_energy_3_plot(θs32,rl_indexes;y=-0.4)
fig.savefig("link2_inv_by_energy_3.png",dpi=300,bbox_inches="tight")

fig = position_plot(θs2;row=1,elev=11,azim=-76,y=-0.1)
fig.savefig("link2_inv_pose.png",dpi=300,bbox_inches="tight")

fig = position_plot(θs2;nstage=4,row=1,zup=0.32,elev=11,azim=-76,y=-0.2)
fig.savefig("link4_inv_pose.png",dpi=300,bbox_inches="tight")
