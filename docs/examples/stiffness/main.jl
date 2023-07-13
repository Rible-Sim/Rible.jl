#note -- preamble
using OhMyREPL
using LinearAlgebra
using Statistics
using SparseArrays
using StaticArrays
using ElasticArrays
using TypeSortedCollections
using TypedTables
using Rotations
# import GeometricalPredicates as GP
using CoordinateTransformations
using OffsetArrays
import DifferentialEquations as DE
using BenchmarkTools
using RecursiveArrayTools
using CircularArrays
const CA = CircularArray
using Interpolations
# using CubicSplines
# import FLOWMath
using Makie
import GLMakie as GM
import CairoMakie as CM
import WGLMakie as WM
GM.activate!()
Makie.inline!(false)
using FileIO, MeshIO
using LaTeXStrings
using TexTables
using DataStructures
using Latexify
using PrettyTables
auto_display(false)
using Unitful
using EzXML
using CSV, Tables
using Printf
using StructArrays
using EponymTuples
import GeometryBasics as GB
using Meshing
# using FloatingTableView
import Meshes
using Match
using Cthulhu
using COSMO
using AbbreviatedStackTraces
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true
ENV["JULIA_STACKTRACE_MINIMAL"] = true
using Revise
import TensegrityRobots as TR
cd(@__DIR__)
include("../vis.jl"); includet("../vis.jl")
include("../analysis.jl"); includet("../analysis.jl")
include("../dyn.jl"); includet("../dyn.jl")
include("../nonsmooth/def.jl"); includet("../nonsmooth/def.jl")
figdir::String = ""
if Sys.iswindows()
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\TensegrityStability\ES"
    # figdir::String =raw"C:\Users\luo22\OneDrive\Papers\Ph.D.Thesis\dyn"
elseif Sys.isapple()
    figdir::String = raw"."
end
#-- preamble end
unibot = uni(0.0;
            μ = 0.9,
            e = 0.0,
			z0 = 0.2
)
plot_traj!(unibot)
TR.check_static_equilibrium_output_multipliers(unibot.tg)

function static_kinematic_determine(
    B,
    rtol::Real = min(size(B)...)*eps(real(float(one(eltype(B)))))
    )
    ndof,ncables = size(B)
    D = svd(B;)
    (;U,V,S) = D
    tol =  rtol*S[1]
    rank_B = count(x -> x > tol, S)
    S_nonzero = S[begin:rank_B]
    dsi = ncables - rank_B
    dki = ndof - rank_B
    @show dsi,dki
    self_stress_states = V[:,rank_B+1:rank_B+dsi]
    stiffness_directions = U[:,rank_B+1:rank_B+dki]
    self_stress_states,stiffness_directions
end

unibot.tg

q = TR.get_q(unibot.tg)
q̌ = TR.get_q̌(unibot.tg)
Ǎ = TR.make_A(unibot.tg)(q)
Ň = TR.nullspace(Ǎ)
Q̃ = TR.build_Q̃(unibot.tg)
L̂ = TR.build_L̂(unibot.tg)

# Left hand side
Q̃L̂ = Q̃*L̂

B = -Q̃L̂
ℬ = transpose(Ň)*B

s,d = static_kinematic_determine(ℬ)
f = s[:,1]*100
λ = inv(Ǎ*transpose(Ǎ))*Ǎ*B*f
Ǩa = - TR.∂Aᵀλ∂q̌(unibot.tg,λ)
k = TR.get_cables_stiffness(unibot.tg)
Ǩm, Ǩg = TR.build_Ǩm_Ǩg!(unibot.tg,q,f,k)
𝒦m = transpose(Ň)*Ǩm*Ň
𝒦g = transpose(Ň)*Ǩg*Ň
𝒦a = transpose(Ň)*Ǩa*Ň
#note zero material stiffness at the kinematic indeterminate direction(s)
eigen(𝒦m)
#note geometric stiffness scaled with prestress level
d'*𝒦m*d
eigen(𝒦g)
d'*𝒦g*d
eigen(𝒦a)
d'*𝒦a*d

𝒦 = 𝒦m .+ 𝒦g .+ 𝒦a
d'*𝒦*d

eigen_result = eigen(𝒦)
nn = count(x -> x < 0, eigen_result.values)
if nn > 1
    @warn "Instability detected! Number of negative eigenvalues: $nn"
    isstable = false
else
    isstable = true
end
isstable, Ň0, eigen_result



(;U,Vt,S) = D
tol =  rtol*S[1]
count(x -> x > tol, S)

rtol::Real = min(size(B)...)*eps(real(float(one(eltype(B)))))

rank(B)
nullspace(B)
Vend = svdB.Vt[end,:]
B*Vend

unibot.tg.ndof