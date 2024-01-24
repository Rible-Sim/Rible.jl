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
