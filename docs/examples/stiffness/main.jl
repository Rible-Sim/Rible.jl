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
import Polyhedra
import CDDLib 
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
include("../ES/def.jl"); includet("../ES/def.jl")
include("../ES/def3d.jl"); includet("../ES/def3d.jl")


include("examples.jl"); includet("examples.jl")
figdir::String = ""
if Sys.iswindows()
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\TensegrityStability\ES"
elseif Sys.isapple()
    figdir::String = raw"."
end
#-- preamble end
unibot = uni(0.0;
            μ = 0.9,
            e = 0.0,
			z0 = 0.2
)
bot = unibot

spine3dbot = spine3d(2;)
bot = spine3dbot

two = two_tri()
bot = two

ballbot = superball(;constrained=true)
bot = ballbot

tbbot = Tbars()
bot = tbbot

TR.check_static_equilibrium_output_multipliers(bot.tg)

TR.get_cables_tension(bot)
cb1 = bot.tg.tensiles.cables[1]
cb1.state.tension

set_theme!(theme_pub;fontsize = 16 |> pt2px)
plot_traj!(bot;showground=false)
bot.tg.ndof

Makie.inline!(false)
GM.activate!();

function static_kinematic_determine(
    B,
    rtol::Real = min(size(B)...)*eps(real(float(one(eltype(B)))))
    )
    ndof,ncables = size(B)
    D = svd(B; full = true)
    (;U,V,S) = D
    tol =  rtol*S[1]
    rank_B = count(x -> x > tol, S)
    S_nonzero = S[begin:rank_B]
    dsi = ncables - rank_B
    dki = ndof - rank_B
    # @show ndof,ncables,rank_B,dsi,dki
    A = vcat(
        B,
        Matrix(-1I,ncables,ncables)
    )

    hr = Polyhedra.hrep(A, zeros(ndof+ncables),  BitSet(1:ndof))
    ph = Polyhedra.polyhedron(hr, CDDLib.Library(:float))
    vr = Polyhedra.vrep(ph)
    @assert Polyhedra.npoints(vr) == 1
    @show Polyhedra.nrays(vr),dsi
    @assert Polyhedra.nrays(vr) == dsi
    self_stress_states = reduce(hcat,[ray.a for ray in Polyhedra.rays(vr)])
    # display(V[:,rank_B+1:rank_B+dsi])
    stiffness_directions = U[:,rank_B+1:rank_B+dki]
    self_stress_states,stiffness_directions
end

q = TR.get_q(bot.tg)
q̌ = TR.get_q̌(bot.tg)
Ǎ = TR.make_A(bot.tg)(q)
Ň = TR.nullspace(Ǎ)
Q̃ = TR.build_Q̃(bot.tg)
L̂ = TR.build_L̂(bot.tg)

# Left hand side
Q̃L̂ = Q̃*L̂

B = -Q̃L̂
ℬ = transpose(Ň)*B

s,d = static_kinematic_determine(ℬ)
s 
d
k = TR.get_cables_stiffness(bot.tg)
l = TR.get_cables_len(bot.tg)
f = s[:,2]
# equivalent μ
μ = l .- (f./k)

λ = inv(Ǎ*transpose(Ǎ))*Ǎ*B*f
Ǩa = - TR.∂Aᵀλ∂q̌(bot.tg,λ)

Ǩm, Ǩg = TR.build_Ǩm_Ǩg!(bot.tg,q,f,k)
𝒦m = transpose(Ň)*Ǩm*Ň
𝒦g = transpose(Ň)*Ǩg*Ň
𝒦a = transpose(Ň)*Ǩa*Ň
#note zero material stiffness at the kinematic indeterminate direction(s)
eigen(𝒦m)
#note geometric stiffness scaled with prestress level
d'*𝒦m*d
d'*𝒦g*d
d'*𝒦a*d

vals_𝒦g, vecs_𝒦g = eigen(𝒦g)

rank(𝒦g .+ 𝒦a)

issymmetric(𝒦a)

𝒦a .- 𝒦a' .|> norm |> maximum

vals_𝒦ga, vecs_𝒦ga = eigen(𝒦g .+ 𝒦a)

vals_𝒦a, vecs_𝒦a = eigen(𝒦a)

vals_𝒦ga

eigen(𝒦a)

𝒦 = 𝒦m .+ 𝒦g .+ 𝒦a

𝒦 = 𝒦m .+ α.*(𝒦g .+ 𝒦a)

d = eigvecs(𝒦g .+ 𝒦a)[:,1]

d'*(𝒦g .+ 𝒦a)*d

α = - 𝒦m ./(𝒦g .+ 𝒦a)

d'*𝒦*d

eigen_result = eigen(𝒦)





(;U,Vt,S) = D
tol =  rtol*S[1]
count(x -> x > tol, S)

rtol::Real = min(size(B)...)*eps(real(float(one(eltype(B)))))

rank(B)
nullspace(B)
Vend = svdB.Vt[end,:]
B*Vend

unibot.tg.ndof