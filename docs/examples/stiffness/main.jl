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
# using LDLFactorizations
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
using Arpack
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
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\TensegrityStability"
elseif Sys.isapple()
    figdir::String = raw"."
end

fontsize = 8 |> pt2px
tw = 468 |> pt2px
th = 622 |> pt2px

# for use with Class-1 and the 1st rigid fixed
function build_Ň(tg)
    (;sysfree,mem2sysndof) = tg.connectivity.indexed
    q = TR.get_q(bot.tg)
    Nin = TR.make_intrinsic_nullspace(tg,q)
    Nin[
        sysfree,
        reduce(vcat,mem2sysndof[2:end])
    ]
end
#-- preamble end

bot = deepcopy(botinput)
q = TR.get_q(bot.tg)
q̌ = TR.get_q̌(bot.tg)
Ǎ = TR.make_A(bot.tg)(q)
Ň_ = TR.nullspace(Ǎ)
Ň = modified_gram_schmidt(Ň_)
Ň = Ň_
@myshow Ǎ*Ň |> norm

Q̃ = TR.build_Q̃(bot.tg)
L̂ = TR.build_L̂(bot.tg)

# Left hand side
Q̃L̂ = Q̃*L̂

Bᵀ = -Q̃L̂
ℬᵀ = transpose(Ň)*Bᵀ
ℬ = transpose(ℬᵀ)
S,D = static_kinematic_determine(ℬᵀ;atol=1e-14)
S 
D
nk = size(D,2)
δq̌ = [Ň*D[:,i] for i in axes(D,2)]
δq̌ = [Ň*vecs_𝒦[:,i] for i in axes(D,2)]
scaling=0.1
botvis = deepcopy(bot)
for i = 1:nk
    push!(botvis.traj,deepcopy(botvis.traj[end]))
    botvis.traj.t[end] = i
    δq̌i = δq̌[i]
    ratio = norm(δq̌i)/norm(q̌)
    botvis.traj.q̌[end] .= q̌ .+ scaling.*δq̌i/ratio
end



#-- uni bot
unibot = uni(0.0;
            μ = 0.9,
            e = 0.0,
			z0 = 0.2
)
bot = unibot
plot_traj!(bot;showground=false)
plot_rigid(TR.get_rigidbodies(bot)[2])


plot_traj!(bot;showground=false)
bot.tg.ndof

TR.check_static_equilibrium_output_multipliers(bot.tg)

q = TR.get_q(bot.tg)
q̌ = TR.get_q̌(bot.tg)
Ǎ = TR.make_A(bot.tg)(q)
Ň_ = TR.nullspace(Ǎ)
Ň = modified_gram_schmidt(Ň_)
Ň = Ň_
# N = TR.make_intrinsic_nullspace(bot.tg,q)

Ň = build_Ň(bot.tg)

rank(Ň)

Ǎ*Ň |> norm

Q̃ = TR.build_Q̃(bot.tg)
L̂ = TR.build_L̂(bot.tg)

# Left hand side
Q̃L̂ = Q̃*L̂

Bᵀ = -Q̃L̂
ℬᵀ = transpose(Ň)*Bᵀ
ℬ = transpose(ℬᵀ)
S,D = static_kinematic_determine(ℬᵀ;atol=1e-14)
S 
plot_self_stress_states(bot,S)
D
nk = size(D,2)
δq̌ = [Ň*D[:,i] for i in axes(D,2)]
δq̌ = [Ň*vecs_𝒦[:,i] for i in axes(D,2)]
scaling=0.1
botvis = deepcopy(bot)
for i = 1:nk
    push!(botvis.traj,deepcopy(botvis.traj[end]))
    botvis.traj.t[end] = i
    δq̌i = δq̌[i]
    ratio = norm(δq̌i)/norm(q̌)
    botvis.traj.q̌[end] .= q̌ .+ scaling.*δq̌i/ratio
end

plot_traj!(botvis;)

k = TR.get_cables_stiffness(bot.tg)
l = TR.get_cables_len(bot.tg)

f =  sum(S,dims=2)

# equivalent μ
# μ = l .- (f./k)

λ = -inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*f
# @show f,λ
using Symbolics
@variables λ[1:6]
rb2 = TR.get_rigidbodies(bot.tg)[2]
rb2.state.cache.funcs.∂Aᵀλ∂q(λ)#[:,free_idx]

Ǎ*Ǎ'

Ǩa = TR.∂Aᵀλ∂q̌(bot.tg,Symbolics.scalarize(λ))
𝒦a = transpose(Ň)*Ǩa*Ň 
vals_𝒦a,vecs_𝒦a = eigen(𝒦a)
sort(vals_𝒦a)
memincst = [1, 2, 3, 4, 5, 6]
Symbolics.unwrap(λ)[memincst...]
[memincst]

# @show count((x)->x<0,D_𝒦a)
# @show count((x)->x==0,D_𝒦a)

Ǩm, Ǩg = TR.build_Ǩm_Ǩg!(bot.tg,q,f,k)

𝒦m = transpose(Ň)*Ǩm*Ň |> Symmetric 
𝒦g = transpose(Ň)*Ǩg*Ň |> Symmetric 
𝒦p = 𝒦g.+ 𝒦a |> Symmetric 
𝒦 = 𝒦m.+ 𝒦p |> Symmetric

#note zero material stiffness at the kinematic indeterminate direction(S)
D'*𝒦m*D
D'*𝒦p*D
D'*𝒦*D

vals_𝒦m,vecs_𝒦m = eigen(𝒦m)
sort(vals_𝒦m)
vm = vecs_𝒦m[:,1]

rank(ℬ)

vals_𝒦g,vecs_𝒦g = eigen(𝒦g)
sort(vals_𝒦g)

vals_𝒦p,vecs_𝒦p = eigen(𝒦p)
sort(vals_𝒦p)

vals_𝒦,vecs_𝒦 = eigen(𝒦)
sort(vals_𝒦)

v = vecs_𝒦[:,1]

vm'*𝒦m*vm
vm'*𝒦g*vm
vm'*𝒦a*vm
vm'*𝒦p*vm
rank(vm'*𝒦p*vm)

v'*𝒦*v


spine3dbot = spine3d(2;)
bot = spine3dbot

#-- end uni bot

#-- two triangles
two = two_tri()
bot = two
plot_traj!(bot;showground=false)
bot.tg.ndof

TR.check_static_equilibrium_output_multipliers(bot.tg)

q = TR.get_q(bot.tg)
q̌ = TR.get_q̌(bot.tg)
Ǎ = TR.make_A(bot.tg)(q)
Ň_ = TR.nullspace(Ǎ)
Ň = modified_gram_schmidt(Ň_)
Ň = Ň_
# N = TR.make_intrinsic_nullspace(bot.tg,q)

function make_Ň(tg,q)
	N = TR.make_intrinsic_nullspace(tg,q)
    I2 = TR.NCF.I2_Int
    O2 = zero(I2)
    o2 = O2[:,1]
    ret = N*[
        I2    o2    o2;
        o2'    1    0;
        I2    o2    o2;
        o2'    0    1;
    ]
    # ret = N*[
    #     I2    o2    o2;
    #     o2'    1    0;
    #     I2    o2    o2;
    #     o2'    0    1;
    # ]
    ret
end

Ň = make_Ň(bot.tg,q)

rank(Ň)

Ǎ*Ň |> norm

Q̃ = TR.build_Q̃(bot.tg)
L̂ = TR.build_L̂(bot.tg)

# Left hand side
Q̃L̂ = Q̃*L̂

Bᵀ = -Q̃L̂
ℬᵀ = transpose(Ň)*Bᵀ
ℬ = transpose(ℬᵀ)
S,D = static_kinematic_determine(ℬᵀ;atol=1e-14)
S 
D
nk = size(D,2)
δq̌ = [Ň*D[:,i] for i in axes(D,2)]
δq̌ = [Ň*vecs_𝒦[:,i] for i in axes(D,2)]
scaling=0.1
botvis = deepcopy(bot)
for i = 1:nk
    push!(botvis.traj,deepcopy(botvis.traj[end]))
    botvis.traj.t[end] = i
    δq̌i = δq̌[i]
    ratio = norm(δq̌i)/norm(q̌)
    botvis.traj.q̌[end] .= q̌ .+ scaling.*δq̌i/ratio
end

plot_traj!(botvis;)

k = TR.get_cables_stiffness(bot.tg)
l = TR.get_cables_len(bot.tg)

f =  S[:,4]

# equivalent μ
# μ = l .- (f./k)

λ = -inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*f
# @show f,λ
Ǩa = TR.∂Aᵀλ∂q̌(bot.tg,λ)
𝒦a = transpose(Ň)*Ǩa*Ň |> Symmetric 
vals_𝒦a,vecs_𝒦a = eigen(𝒦a)
sort(vals_𝒦a)

# @show count((x)->x<0,D_𝒦a)
# @show count((x)->x==0,D_𝒦a)

Ǩm, Ǩg = TR.build_Ǩm_Ǩg!(bot.tg,q,f,k)

𝒦m = transpose(Ň)*Ǩm*Ň |> Symmetric 
𝒦g = transpose(Ň)*Ǩg*Ň |> Symmetric 
𝒦p = 𝒦g.+ 𝒦a |> Symmetric 
𝒦 = 𝒦m.+ 𝒦p |> Symmetric

#note zero material stiffness at the kinematic indeterminate direction(S)
D'*𝒦m*D
D'*𝒦p*D
D'*𝒦*D

vals_𝒦m,vecs_𝒦m = eigen(𝒦m)
sort(vals_𝒦m)
vm = vecs_𝒦m[:,1]

rank(ℬ)

vals_𝒦g,vecs_𝒦g = eigen(𝒦g)
sort(vals_𝒦g)

vals_𝒦p,vecs_𝒦p = eigen(𝒦p)
sort(vals_𝒦p)

vals_𝒦,vecs_𝒦 = eigen(𝒦)
sort(vals_𝒦)

v = vecs_𝒦[:,1]

vm'*𝒦m*vm
vm'*𝒦g*vm
vm'*𝒦a*vm
vm'*𝒦p*vm
rank(vm'*𝒦p*vm)

v'*𝒦*v

#-- end two triangles

#--- superball
ballbot = superball(;
    θ = 0.0,
    l = 2.0/2,
    d = 2.0/2/2,
    z0 = 2.0/2,
    constrained=true
)
bot = ballbot
rb1 = TR.get_rigidbodies(bot)[1]

# ballbot = superball(;
#     θ = 0.0,
#     l = 2.0/2,
#     d = 2.0/2/2,
#     z0 = 2.0/2,
#     constrained=true,
#     addconst = reshape(δq̌,1,:)
# )
# bot = ballbot

TR.check_static_equilibrium_output_multipliers(bot.tg)

TR.update!(bot.tg)
f = TR.get_cables_tension(bot)

function verify_lambda(tg)
    T = TR.get_numbertype(tg)
    λs = zeros(T,tg.nbodies)
    foreach(tg.bodies) do rb
        (;prop,state) = rb
        (;rps,ro,fps) = state
        @myshow prop.id
        @myshow rps
        @myshow fps
        for (rp,fp) in zip(rps,fps)
            λs[prop.id] += 1/2*(rp-ro)'*fp
        end
    end
    λs
end
verify_lambda(bot.tg)
q = TR.get_q(bot.tg)
q̌ = TR.get_q̌(bot.tg)
Ǎ = TR.make_A(bot.tg)(q)
# Ň_ = TR.nullspace(Ǎ)
# Ň = modified_gram_schmidt(Ň_)
Ň = build_Ň(bot.tg)
Q̃ = TR.build_Q̃(bot.tg)
L̂ = TR.build_L̂(bot.tg)

# Left hand side
Q̃L̂ = Q̃*L̂

Bᵀ = -Q̃L̂
ℬᵀ = transpose(Ň)*Bᵀ

S,D = TR.static_kinematic_determine(ℬᵀ)
S 
D

ns = size(S,2)
nk = size(D,2)
k = TR.get_cables_stiffness(bot.tg)
l = TR.get_cables_len(bot.tg)

f = S[:,1]# + S[:,2] + S[:,3] + S[:,4]

# equivalent μ
# μ = l .- (f./k)

λ = -inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*f
# @show f,λ
Ǩa = TR.∂Aᵀλ∂q̌(bot.tg,λ)
𝒦a = transpose(Ň)*Ǩa*Ň |> Symmetric 
vals_𝒦a,vecs_𝒦a = eigen(𝒦a)
@myshow sort(vals_𝒦a)
@myshow 𝒦a[1:5,1:5]
# @show count((x)->x<0,D_𝒦a)
# @show count((x)->x==0,D_𝒦a)

Ǩm = TR.build_Ǩm!(bot.tg,q,k)
Ǩg = TR.build_Ǩg!(bot.tg,q,f)

vec𝒦ps = [
    begin
        si = S[:,i]
        # s = S\f
        # @show s
        λi = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*si
        # @show f,λ
        Ǩai = - TR.∂Aᵀλ∂q̌(bot.tg,λi)

        Ǩgi = TR.build_Ǩg!(bot.tg,q,si)

        𝒦pi = transpose(Ň)*(Ǩgi.+Ǩai)*Ň |> Symmetric 
        # vec𝒦pi = SymmetricPacked(𝒦pi).tri
        vec𝒦pi = vec(𝒦pi)
    end
    for i = 1:ns
]

mat𝒦ps = reduce(hcat,vec𝒦ps)

𝒦m = transpose(Ň)*Ǩm*Ň |> Symmetric 
𝒦g = transpose(Ň)*Ǩg*Ň |> Symmetric 
𝒦p = 𝒦g.+ 𝒦a |> Symmetric 
𝒦 = 𝒦m.+ 𝒦p |> Symmetric

vals_𝒦m,vecs_𝒦m = eigen(𝒦m)
sort(vals_𝒦m)
vm = vecs_𝒦m[:,1:nk]

vals_𝒦,vecs_𝒦 = eigen(𝒦)
sort(vals_𝒦)

v = vecs_𝒦[:,1]
v'*𝒦*v

vm[:,1] = v
orthovm = TR.modified_gram_schmidt(vm)


with_theme(theme_pub;
    resolution = (0.3tw,0.3tw),
    fontsize = 6.5 |> pt2px,
    figure_padding = (2fontsize,fontsize,0,0),
    Axis3 = (        
        azimuth = 3.7555306333269844,
        elevation = 0.3726990816987242,
    )
    ) do 
    botvis = deepcopy(bot)
    δq̌ = [Ň*orthovm[:,i] for i in axes(orthovm,2)]
    scaling=0.3
    for i = 1:nk
        push!(botvis.traj,deepcopy(botvis.traj[end]))
        botvis.traj.t[end] = i
        δq̌i = δq̌[i]
        ratio = norm(δq̌i)/norm(q̌)
        botvis.traj.q̌[end] .= q̌ .+ scaling.*δq̌i/ratio
    end
    plot_traj!(
        botvis;
        figsize = (0.8tw,0.26tw),
        AxisType=Axis3,
        gridsize=(1,nk+1),        
        atsteps=1:nk+1,
        doslide=false,
        showlabels=false,
        showpoints=false,
        # showcables = false,
        showground = false,
        xlims = (-1e0,1e0),
        ylims = (-1e0,1e0),
        zlims = (-1e-5,2e0),
        slack_linestyle = :solid,
        showinit = true,titleformatfunc = (sgi,tt)-> begin
            rich(
                    rich("($(alphabet[sgi])) ", font=:bold),
                    [
                        "Initial",
                        "Mechanism Mode 1",
                        "Mechanism Mode 2"
                    ][sgi]
                )
        end,
        sup! = (ax,tgob,sgi)->begin
            if sgi != 1
                hidedecorations!(ax)
                xlims!(ax,-1.0e0,1.2e0)
                ylims!(ax,-1.2e0,1.0e0)
            end
        end,
        figname="superball"
    )
end


Ňv = Ň*nullspace(v')

Ǩm = TR.build_Ǩm!(bot.tg,q,k) 
r𝒦m = transpose(Ňv)*(Ǩm)*Ňv |> Symmetric 
# vecr𝒦m = SymmetricPacked(r𝒦m).tri
vecr𝒦m = vec(r𝒦m)

# vecI = SymmetricPacked(Matrix(1.0I,size(r𝒦m))).tri
vecI = vec(Matrix(1.0I,size(r𝒦m)))
r𝒦m |> issymmetric

vecr𝒦ps = [
    begin
        si = S[:,i]
        # s = S\f
        # @show s
        λi = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*si
        # @show f,λ
        Ǩai = - TR.∂Aᵀλ∂q̌(bot.tg,λi)

        Ǩgi = TR.build_Ǩg!(bot.tg,q,si)

        r𝒦pi = transpose(Ňv)*(Ǩgi.+Ǩai)*Ňv |> Symmetric 
        # vecr𝒦pi = SymmetricPacked(r𝒦pi).tri
        vecr𝒦pi = vec(r𝒦pi)
    end
    for i = 1:ns
]

matr𝒦ps = reduce(hcat,vecr𝒦ps)

ᾱ = [1.0]
A = hcat(
    -Matrix(1.0I,ns,ns),
    ᾱ,
    zero(ᾱ)
)
b = [0.0]
nx = ns+2
result_max = TR.optimize_maximum_stiffness(matr𝒦ps,vecr𝒦m,vecI,A,b,nx)
σ_max = result_max.x[end-1]
ρ_max = result_max.x[end]

r𝒦_max = r𝒦m + σ_max*reshape(matr𝒦ps*ᾱ,size(r𝒦m))
vals_r𝒦_max, vecs_r𝒦_max = eigen(r𝒦_max)

vals, vecs = eigen(r𝒦_max - ρ_max*I)
@myshow vals

result_min = TR.optimize_minimum_stiffness(matr𝒦ps,vecr𝒦m,vecI,
    hcat(
        -Matrix(1.0I,ns,ns),
        ᾱ,
    ),
    [0.0],
    ns+1,
    # result.x
)
σ_min = result_min.x[end]

r𝒦_min = r𝒦m + σ_min*reshape(matr𝒦ps*ᾱ,size(r𝒦m))
vals_r𝒦_min, vecs_r𝒦_min = eigen(r𝒦_min)
ρ_min = vals_r𝒦_min[1]
maxminmodes = hcat(
    vecs_r𝒦_max[:,1],
    vecs_r𝒦_min[:,1:3],
)

with_theme(theme_pub;
    fontsize = 6.5 |> pt2px,
    figure_padding = (0,0,-fontsize,0),
    Axis3 = (        
        azimuth = 3.7555306333269844,
        elevation = 0.3726990816987242,
    )
    ) do 
    botvis = deepcopy(bot)
    δq̌ = [Ňv*maxminmodes[:,i] for i in axes(maxminmodes,2)]
    scaling=0.3
    for i = 1:4
        push!(botvis.traj,deepcopy(botvis.traj[end]))
        botvis.traj.t[end] = i
        δq̌i = δq̌[i]
        ratio = norm(δq̌i)/norm(q̌)
        botvis.traj.q̌[end] .= q̌ .+ scaling.*δq̌i/ratio
    end
    plot_traj!(
        botvis;
        figsize = (0.8tw,0.18tw),
        AxisType=Axis3,
        gridsize=(1,4),        
        atsteps=1+1:4+1,
        doslide=false,
        showlabels=false,
        showpoints=false,
        # showcables = false,
        showground = false,
        xlims = (-1e0,1e0),
        ylims = (-1e0,1e0),
        zlims = (-1e-5,2e0),
        slack_linestyle = :solid,
        showinit = true,titleformatfunc = (sgi,tt)-> begin
            rich(
                    rich("($(alphabet[sgi])) ", font=:bold),
                    [
                        rich("Mode 1 at σ", subscript("max")),
                        rich("Mode 1 at σ", subscript("min")),
                        rich("Mode 2 at σ", subscript("min")),
                        rich("Mode 3 at σ", subscript("min")),
                    ][sgi]
                )
        end,
        sup! = (ax,tgob,sgi)->begin
            hidedecorations!(ax)
            xlims!(ax,-1.0e0,1.2e0)
            ylims!(ax,-1.2e0,1.0e0)
        end,
        figname="superball_maxmin"
    )
end


σs = 0:100:5200
rρs =  [
    begin
        r𝒦 = r𝒦m + σ*reshape(matr𝒦ps*ᾱ,size(r𝒦m))
        vals_r𝒦, vecs_r𝒦 = eigen(r𝒦)
        vals_r𝒦[begin]
    end
    for σ in σs
]

ρs =  [
    begin
        𝒦 = 𝒦m + σ*reshape(mat𝒦ps*ᾱ,size(𝒦m))
        vals_𝒦, vecs_𝒦 = eigen(𝒦)
        vals_𝒦[begin+1]
    end
    for σ in σs
]

with_theme(theme_pub;
        resolution = (0.3tw,0.2tw),
        figure_padding = (0,fontsize,0,fontsize),
    ) do 
    fig = Figure()
    ax = Axis(fig[1,1],
        xlabel = L"\sigma",
        ylabel = L"\rho_{\mathrm{1}}"
    )
    lines!(ax,σs,rρs,)
    xlims!(ax,0,5500)
    ylims!(ax,-400,600)
    scatter!(
        ax,
        [σ_max,σ_min],
        [ρ_max,ρ_min]
    )
    text!([σ_max], [ρ_max], 
        text = [L"\sigma_{\mathrm{max}}"],
        align = (:center,:bottom),
        offset = (0, fontsize/4)
    )
    text!([σ_min], [ρ_min], 
        text = [L"\sigma_{\mathrm{min}}"],
        align = (:right,:center),
        offset = (-fontsize/2, 0)
    )
    # text!(x, y, text = string.(aligns), align = aligns)
    savefig(fig,"superball_curve")
    fig
end

#-- superball end

#-- prism begin
m = 3
α = 2π/m
θ = 1.25α
n = 4
b = 0.14
r = 0.04*sqrt(2)
prism1 = prisms(;
    r1= 0.03*sqrt(2),
    r,b,m,α,θ,n = 1,
)
bot = prism1
rb1 = TR.get_rigidbodies(bot)[1]
plot_rigid(rb1) 
plot_traj!(bot;showground=false)
TR.check_static_equilibrium_output_multipliers(bot.tg)
@myshow bot.tg.ndof
TR.update!(bot.tg)
f = TR.get_cables_tension(bot)

# for use with Class-1 and the 1st rigid fixed
function build_Ň(tg)
    (;sysfree,mem2sysfull,mem2sysndof) = tg.connectivity.indexed
    q = TR.get_q(bot.tg)
    Nin = TR.make_intrinsic_nullspace(tg,q)[
        sysfree,
        reduce(vcat,mem2sysndof[begin:end-1])
    ]
    Nex = zeros(eltype(q),30,12)
    for i = 1:6
        is = (i-1)*5
        js = (i-1)*2
        Nex[is+4:is+5,js+1:js+2] .= Matrix(1I,2,2)
    end
    cm = CircularArray(collect(1:3))
    for i = 1:3
        is = (3+cm[i+1]-1)*5
        js = (i-1)*2
        q_I = q[mem2sysfull[i]]
        ri = @view q_I[1:3]
        u = @view q_I[4:6]
        v,w = TR.NCF.HouseholderOrthogonalization(u)
        @myshow i,3+cm[i+1],u,v,w
        # R = [u v w;]
        Nex[is+1:is+3,js+1:js+2] = -TR.NCF.skew(0.14u)*[w -v;]
    end
    Nin*Nex
end

function verify_lambda(tg)
    T = TR.get_numbertype(tg)
    λs = zeros(T,tg.nbodies)
    foreach(tg.bodies) do rb
        (;prop,state) = rb
        (;rps,ro,fps) = state
        @myshow prop.id
        @myshow rps
        @myshow fps
        for (rp,fp) in zip(rps,fps)
            λs[prop.id] += 1/2*(rp-ro)'*fp
        end
    end
    λs
end
verify_lambda(bot.tg)
q = TR.get_q(bot.tg)
q̌ = TR.get_q̌(bot.tg)
Ǎ = TR.make_A(bot.tg)(q)
# Ň_ = TR.nullspace(Ǎ)
# Ň = TR.modified_gram_schmidt(Ň_)
Ň = build_Ň(bot.tg)
Q̃ = TR.build_Q̃(bot.tg)
L̂ = TR.build_L̂(bot.tg)
rank(Ň)
Ǎ*Ň |> norm
# Left hand side
Q̃L̂ = Q̃*L̂

Bᵀ = -Q̃L̂
ℬᵀ = transpose(Ň)*Bᵀ

S,D = TR.static_kinematic_determine(ℬᵀ)
S 
D
plot_self_stress_states(
    bot,
    S,
    rtol = 1e-10
    # Ň = build_Ň(bot.tg)
)

ns = size(S,2)
nk = size(D,2)
k = TR.get_cables_stiffness(bot.tg)
l = TR.get_cables_len(bot.tg)

f = S[:,1]# + S[:,2] + S[:,3] + S[:,4]

# equivalent μ
# μ = l .- (f./k)

λ = -inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*f
# @show f,λ
Ǩa = TR.∂Aᵀλ∂q̌(bot.tg,λ)
𝒦a = transpose(Ň)*Ǩa*Ň |> Symmetric 
vals_𝒦a,vecs_𝒦a = eigen(𝒦a)
@myshow sort(vals_𝒦a)
@myshow 𝒦a[1:5,1:5]
# @show count((x)->x<0,D_𝒦a)
# @show count((x)->x==0,D_𝒦a)

Ǩm = TR.build_Ǩm!(bot.tg,q,k)
Ǩg = TR.build_Ǩg!(bot.tg,q,f)

𝒦m = transpose(Ň)*Ǩm*Ň |> Symmetric
vec𝒦m  = vec(𝒦m)
vecI = vec(Matrix(1.0I,size(𝒦m)))
𝒦g = transpose(Ň)*Ǩg*Ň |> Symmetric

𝒦p = 𝒦g.+ 𝒦a |> Symmetric

vals_𝒦p,vecs_𝒦p = eigen(𝒦p)
@myshow sort(vals_𝒦p)

𝒦 = 𝒦m.+ 𝒦p |> Symmetric

vec𝒦ps = [
    begin
        si = S[:,i]
        # s = S\f
        # @show s
        λi = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*si
        # @show f,λ
        Ǩai = - TR.∂Aᵀλ∂q̌(bot.tg,λi)

        Ǩgi = TR.build_Ǩg!(bot.tg,q,si)

        𝒦pi = transpose(Ň)*(Ǩgi.+Ǩai)*Ň |> Symmetric 
        # vec𝒦pi = SymmetricPacked(𝒦pi).tri
        vec𝒦pi = vec(𝒦pi)
    end
    for i = 1:ns
]

mat𝒦ps = reduce(hcat,vec𝒦ps)

ᾱ = [1.0,1.0]
A = hcat(
    -Matrix(1.0I,ns,ns),
    ᾱ,
    zero(ᾱ)
)
b = [0.0,0.0]
nx = ns+2
result_max = TR.optimize_maximum_stiffness(mat𝒦ps,vec𝒦m,vecI,A,b,nx)
σ_max = result_max.x[end-1]
ρ_max = result_max.x[end]

𝒦_max = 𝒦m + σ_max*reshape(mat𝒦ps*ᾱ,size(𝒦m))
vals_𝒦_max, vecs_𝒦_max = eigen(𝒦_max)

vals, vecs = eigen(𝒦_max - ρ_max*I)
@myshow vals

result_min = TR.optimize_minimum_stiffness(mat𝒦ps,vec𝒦m,vecI,
    hcat(
        -Matrix(1.0I,ns,ns),
        ᾱ,
    ),
    [0.0,0.0],
    ns+1,
    # result.x
)
σ_min = result_min.x[end]

𝒦_min = 𝒦m + σ_min*reshape(mat𝒦ps*ᾱ,size(𝒦m))
vals_𝒦_min, vecs_𝒦_min = eigen(𝒦_min)
ρ_min = vals_𝒦_min[1]
maxminmodes = hcat(
    vecs_𝒦_max[:,1],
    vecs_𝒦_min[:,1:3],
)

σs = 0:10:1600
ρs =  [
    begin
        𝒦 = 𝒦m + σ*reshape(mat𝒦ps*ᾱ,size(𝒦m))
        vals_𝒦, vecs_𝒦 = eigen(𝒦)
        vals_𝒦[begin]
    end
    for σ in σs
]

with_theme(theme_pub;
        resolution = (0.3tw,0.2tw),
        figure_padding = (0,fontsize,0,fontsize),
    ) do 
    fig = Figure()
    ax = Axis(fig[1,1],
        xlabel = L"\sigma",
        ylabel = L"\rho_{\mathrm{1}}"
    )
    lines!(ax,σs,ρs,)
    # xlims!(ax,0,5500)
    # ylims!(ax,-400,600)
    scatter!(
        ax,
        [σ_max,σ_min],
        [ρ_max,ρ_min]
    )
    text!([σ_max], [ρ_max], 
        text = [L"\sigma_{\mathrm{max}}"],
        align = (:center,:bottom),
        offset = (0, fontsize/4)
    )
    text!([σ_min], [ρ_min], 
        text = [L"\sigma_{\mathrm{min}}"],
        align = (:right,:center),
        offset = (-fontsize/2, 0)
    )
    # text!(x, y, text = string.(aligns), align = aligns)
    # savefig(fig,"superball_curve")
    fig
end


#-- prism end
tbbot = Tbars()
bot = tbbot

# ULTRA Spine 2D
spine2d = make_spine(2)
bot = spine2d


pp = planar_parallel()
bot = pp



bot.tg.connectivity.numbered.mem2sys
plot_traj!(bot;showground=false)
bot.tg.ndof

Makie.inline!(false)
GM.activate!();

TR.check_static_equilibrium_output_multipliers(bot.tg)


super_stability(bot)

q = TR.get_q(bot.tg)
Ǎ = TR.make_A(bot.tg)(q)
# Ň_ = TR.nullspace(Ǎ)
# Ň = modified_gram_schmidt(Ň_)
Ň = build_Ň(bot.tg)
Q̃ = TR.build_Q̃(bot.tg)
L̂ = TR.build_L̂(bot.tg)

# Left hand side
Q̃L̂ = Q̃*L̂

Bᵀ = -Q̃L̂
ℬᵀ = transpose(Ň)*Bᵀ
rank
S,D = static_kinematic_determine(ℬᵀ)
S 
D
k = TR.get_cables_stiffness(bot.tg)
l = TR.get_cables_len(bot.tg)

f = S[:,1]# + S[:,2] + S[:,3] + S[:,4]

# equivalent μ
# μ = l .- (f./k)

λ = -inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*f
# @show f,λ
Ǩa = TR.∂Aᵀλ∂q̌(bot.tg,λ)
𝒦a = transpose(Ň)*Ǩa*Ň |> Symmetric 
vals_𝒦a,vecs_𝒦a = eigen(𝒦a)
sort(vals_𝒦a)

# @show count((x)->x<0,D_𝒦a)
# @show count((x)->x==0,D_𝒦a)

Ǩm = TR.build_Ǩm!(bot.tg,q,k)
Ǩg = TR.build_Ǩg!(bot.tg,q,f)

𝒦m = transpose(Ň)*Ǩm*Ň |> Symmetric 
𝒦g = transpose(Ň)*Ǩg*Ň |> Symmetric 

vals_𝒦m,vecs_𝒦m = eigen(𝒦m)

sort(vals_𝒦m)

#note zero material stiffness at the kinematic indeterminate direction(S)
D'*𝒦m*D

vals_𝒦g,vecs_𝒦g = eigen(𝒦g)
sort(vals_𝒦g)

𝒦p = 𝒦g.+ 𝒦a |> Symmetric 
vals_𝒦p,vecs_𝒦p = eigen(𝒦p)
sort(vals_𝒦p)

# @show count((x)->x<0,D_𝒦p)
# @show count((x)->x==0,D_𝒦p)


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

