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
using NLsolve
using SymmetricFormats
using EponymTuples
import GeometryBasics as GB
using Meshing
# using FloatingTableView
import Meshes
using Match
using Cthulhu
using COSMO
using Polyhedra
import CDDLib 
lib = CDDLib.Library()
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
include("../stability/gripper_define.jl"); includet("../stability/gripper_define.jl")
include("examples.jl"); includet("examples.jl")
include("../LC/mydef.jl"); includet("../LC/mydef.jl")
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

ballbot = superball(;θ = 0.0,constrained=true)
bot = ballbot

tbbot = Tbars()
bot = tbbot

# ULTRA Spine 2D
spine2d = make_spine(2)
bot = spine2d

nd1 = man_nd1(;ratio=0.85)
bot = nd1

bot.tg.connectivity.numbered.mem2sys
plot_traj!(bot;showground=false)
bot.tg.ndof

Makie.inline!(false)
GM.activate!();

TR.check_static_equilibrium_output_multipliers(bot.tg)

function static_kinematic_determine(
    ℬᵀ,
    rtol::Real = min(size(ℬᵀ)...)*eps(real(float(one(eltype(ℬᵀ)))))
    )
    ndof,ncables = size(ℬᵀ)
    U,S,V = svd(ℬᵀ; full = true)
    @show size(U),size(V),size(S)
    tol =  rtol*S[1]
    rank_ℬᵀ = count(x -> x > tol, S)
    S_nonzero = S[begin:rank_ℬᵀ]
    dsi = ncables - rank_ℬᵀ
    dki = ndof - rank_ℬᵀ
    @show ndof,ncables,rank_ℬᵀ,dsi,dki
    A = vcat(
        ℬᵀ,
        Matrix(-1I,ncables,ncables)
    )

    hr = hrep(A, zeros(ndof+ncables),  BitSet(1:ndof))
    ph = polyhedron(hr, lib)
    vr = vrep(ph)
    @assert npoints(vr) == 1
    @show nrays(vr)
    self_stress_states = reduce(hcat,[ray.a for ray in rays(vr)])
    # self_stress_states = V[:,rank_ℬᵀ+1:rank_ℬᵀ+dsi]
    stiffness_directions = U[:,rank_ℬᵀ+1:rank_ℬᵀ+dki]
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

Bᵀ = -Q̃L̂
ℬᵀ = transpose(Ň)*Bᵀ

S,D = static_kinematic_determine(ℬᵀ)
S 
D
k = TR.get_cables_stiffness(bot.tg)
l = TR.get_cables_len(bot.tg)

f = S[:,1]
# equivalent μ
μ = l .- (f./k)
l.*k

λ = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*f
# @show f,λ
Ǩa = - TR.∂Aᵀλ∂q̌(bot.tg,λ)
Ǩm, Ǩg = TR.build_Ǩm_Ǩg!(bot.tg,q,f,k)

𝒦m = transpose(Ň)*Ǩm*Ň
𝒦g = transpose(Ň)*Ǩg*Ň
𝒦a = transpose(Ň)*Ǩa*Ň
#note zero material stiffness at the kinematic indeterminate direction(S)
eigen(𝒦m)
#note geometric stiffness scaled with prestress level
D'*𝒦m*D
D'*𝒦g*D
D'*𝒦a*D

𝒦a |> issymmetric

vals_𝒦a, vecs_𝒦a = eigen(Symmetric(𝒦a))
vals_𝒦a
ρs = Float64[]
vs = Vector{Float64}[]
for ratio in range(0.95,1.0,20)
    f = k.*(l.-ratio*l)
    # s = S\f
    # @show s
    λ = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*f
    # @show f,λ
    Ǩa = - TR.∂Aᵀλ∂q̌(bot.tg,λ)
    Ǩm, Ǩg = TR.build_Ǩm_Ǩg!(bot.tg,q,f,k)
    𝒦m = transpose(Ň)*(Ǩm)*Ň
    𝒦G = transpose(Ň)*(Ǩg.+Ǩa)*Ň
    𝒦 = 𝒦m .+ 𝒦G
    genvals_𝒦,_ = eigen(-𝒦G, 𝒦m)
    @show f.*first(genvals_𝒦)
    vals_𝒦, vecs_𝒦 = eigen(𝒦)
    push!(ρs, first(vals_𝒦))
    push!(vs,vecs_𝒦[:,begin])
end


ρs
vs
Makie.scatterlines(range(0.95,1.0,20),ρs)

s1 = S[:,1]
s2 = S[:,2]
# s = S\f
# @show s
λ1 = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*s1
λ2 = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*s2
# @show f,λ
Ǩa1 = - TR.∂Aᵀλ∂q̌(bot.tg,λ1)
Ǩa2 = - TR.∂Aᵀλ∂q̌(bot.tg,λ2)

Ǩm1, Ǩg1 = TR.build_Ǩm_Ǩg!(bot.tg,q,s1,k)
Ǩm2, Ǩg2 = TR.build_Ǩm_Ǩg!(bot.tg,q,s2,k)

𝒦m = transpose(Ň)*(Ǩm1)*Ň

𝒦1 = transpose(Ň)*(Ǩg1.+Ǩa1)*Ň
𝒦2 = transpose(Ň)*(Ǩg2.+Ǩa2)*Ň

vec𝒦1 = SymmetricPacked(𝒦1).tri
vec𝒦2 = SymmetricPacked(𝒦2).tri
vec𝒦m = SymmetricPacked(𝒦m).tri
vecI = SymmetricPacked(Matrix(1.0I,size(𝒦m))).tri
kl = k.*l

model = COSMO.Model()

constraint1 = COSMO.Constraint(
    hcat(
        vec𝒦1,vec𝒦2,zero(vecI),-vecI
    ), 
    vec𝒦m, 
    COSMO.PsdConeTriangle
)

constraint2 = COSMO.Constraint(
    hcat(
        S,kl,zero(k)
    ), 
    -kl, 
    COSMO.ZeroSet
)

constraints = [constraint1,constraint2]

nx = size(S,2)+2

custom_settings = COSMO.Settings(
    verbose = true, 
    eps_abs = 1e-7, 
    eps_rel = 1e-7
)

COSMO.assemble!(
    model, 
    zeros(nx,nx), 
    -Matrix(1.0I,nx,nx)[:,end], 
    constraints,
    # settings = custom_settings
)

results = COSMO.optimize!(model)

results.x[end-1]

𝒦 = [𝒦1,𝒦2]

# zero stiffness

model = COSMO.Model()

constraint1 = COSMO.Constraint(
    hcat(
        vec𝒦1,
        vec𝒦2,
        zero(vecI),
    ), 
    vec𝒦m, 
    COSMO.ZeroSet
)

constraint2 = COSMO.Constraint(
    hcat(
        S,
        kl,
    ), 
    -kl, 
    COSMO.ZeroSet
)

constraints = [constraint1,constraint2]

nx = size(S,2)+2

custom_settings = COSMO.Settings(
    verbose = true, 
    eps_abs = 1e-7, 
    eps_rel = 1e-7
)

COSMO.assemble!(
    model, 
    zeros(nx,nx), 
    -Matrix(1.0I,nx,nx)[:,end], 
    constraints,
    # settings = custom_settings
)

function make_zerofunc(𝒦m,𝒦,kl,S)    
    ns = size(S,2)
    ndof = size(𝒦m,2)
    pinvS = pinv(S)
    function inner_f(x)
        ξ = x[1:ndof]
        s = x[ndof+1:ndof+ns]
        σ = x[end]
        vcat(
            (sum([𝒦[i]*s[i] for i = 1:ns]) .+ 𝒦m)*ξ,
            s .- pinvS*kl .+ σ.*pinvS*kl,
            transpose(ξ)*ξ - 1
        )
    end
end

f = make_zerofunc(𝒦m,𝒦,kl,S)
x0 = vcat(
    ones(size(S,2)),
    zeros(size(𝒦m,2)),
    0.0
)
f(x0)

sol = nlsolve(f, x0, method=:newton)
sol.zero[end]


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

D = eigvecs(𝒦g .+ 𝒦a)[:,1]

D'*(𝒦g .+ 𝒦a)*D

α = - 𝒦m ./(𝒦g .+ 𝒦a)

D'*𝒦*D

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