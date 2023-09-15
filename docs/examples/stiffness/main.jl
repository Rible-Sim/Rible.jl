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
    Nin = TR.
    
    .(tg,q)
    Nin[
        sysfree,
        reduce(vcat,mem2sysndof[2:end])
    ]
end
#-- preamble end

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
        Nex[is+1:is+3,js+1:js+2] = -TR.NCF.skew(0.14u)*[v w;]
    end
    Nin,Nex,Nin*Nex
end

q = TR.get_q(bot.tg)
q̌ = TR.get_q̌(bot.tg)
Ǎ = TR.make_A(bot.tg)(q)
# Ň_ = TR.nullspace(Ǎ)
# Ň = TR.modified_gram_schmidt(Ň_)
Nin,Nex,Ň = build_Ň(bot.tg)
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
ns = size(S,2)
nk = size(D,2)

GM.activate!();with_theme(theme_pub;
        resolution = (0.95tw,0.24tw),
        figure_padding = (2fontsize,0,0,0),
        fontsize = 6.5 |> pt2px,
        Axis3 = (
            azimuth = 3.8255306333269843,
            elevation = 0.2026990816987241
        )
    ) do 
    maxS = maximum(abs.(S))
    rtol = 1e-10
    Sbool = S.> maxS*rtol
    S[.!Sbool] .= 0.0
    fig = Figure()
    gd1 = fig[1,1] = GridLayout()
    gd2 = fig[1,2:3] = GridLayout()
    gd3 = fig[1,4:5] = GridLayout()
    botmm = deepcopy(bot)
    plot_traj!(
        bot;
        fig = gd1,
        AxisType=Axis3,
        doslide=false,
        showlabels=false,
        showpoints=false,
        # showcables = false,
        xlims = (-6e-2,6e-2),
        ylims = (-6e-2,6e-2),
        zlims = (-1e-5,2.2e-1),
        showground = false,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                    rich("($(alphabet[sgi])) ", font=:bold),
                    "Initial"
                )
        end,
        sup! = (ax,tgob,sgi)->begin
            hidex(ax)
            hidey(ax)
            # xlims!(ax,-1.0e0,1.2e0)
            # ylims!(ax,-1.2e0,1.0e0)
        end,
    )
    plot_traj!(
        bot;
        fig = gd2,
        AxisType=Axis3,
        gridsize=(1,ns), 
        doslide=false,
        showlabels=false,
        showpoints=false,
        showcables = false,
        xlims = (-6e-2,6e-2),
        ylims = (-6e-2,6e-2),
        zlims = (-1e-5,2.1e-1),
        showground = false,
        showinit = true,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                    rich("($(alphabet[sgi+1])) ", font=:bold),
                    [
                        "Self-stress State 1",
                        "Self-stress State 2"
                    ][sgi]
                )
        end,
        sup! = (ax,tgob,sgi)-> begin
            # cables
            ax.azimuth = 4.73553063332698
            ax.elevation = 0.18269908169872395
            # azimuth = 4.665530633326984
            # elevation = 0.16269908169872424
            hidexyz(ax)
            @myshow Sbool[:,sgi]
            linesegs_cables = @lift begin
                get_linesegs_cables($tgob;)[Sbool[:,sgi]]
            end
            linesegments!(ax, 
                linesegs_cables, 
                color = :red, 
                # linewidth = cablewidth
                )
            rcs_by_cables = @lift begin
                (;tensioned) = $tgob.connectivity
                ndim = TR.get_ndim($tgob)
                T = TR.get_numbertype($tgob)
                ret = Vector{MVector{ndim,T}}()
                mapreduce(
                    (scnt)->
                    [(
                        scnt.end1.rbsig.state.rps[scnt.end1.pid].+
                        scnt.end2.rbsig.state.rps[scnt.end2.pid]
                    )./2],
                    vcat,
                    tensioned.connected
                    ;init=ret
                )
            end
            # @show rcs_by_cables
            Stext = [
                    @sprintf "%4.2f"  S[i,sgi] 
                    for i in axes(S,1)
                    if Sbool[i,sgi]
                ]
            @myshow Stext
            # scatter!(
            #     ax,
            #     rcs_by_cables[][Sbool[:,sgi]],
            #     marker = :rect, 
            #     markersize = 12 |> pt2px, 
            #     color = :white
            # )
            text!(
                ax,
                Stext,
                position = rcs_by_cables[][Sbool[:,sgi]],
                fontsize = 5 |> pt2px,
                color = :red,
                align = (:center, :center),
                # offset = (-fontsize/2, 0)
            )
        end
    )
    δq̌ = [Ň*D[:,i] for i in axes(D,2)]
    scaling=0.3
    for i = 1:nk
        push!(botmm.traj,deepcopy(botmm.traj[end]))
        botmm.traj.t[end] = i
        δq̌i = δq̌[i]
        ratio = norm(δq̌i)/norm(q̌)
        botmm.traj.q̌[end] .= q̌ .+ scaling.*δq̌i/ratio
    end
    plot_traj!(
        botmm;
        fig = gd3,
        AxisType=Axis3,
        gridsize=(1,nk),        
        atsteps=1+1:nk+1,
        doslide=false,
        showlabels=false,
        showpoints=false,
        # showcables = false,
        showground = false,
        xlims = (-8e-2,8e-2),
        ylims = (-8e-2,8e-2),
        zlims = (-1e-5,2.5e-1),
        slack_linestyle = :solid,
        showinit = true,titleformatfunc = (sgi,tt)-> begin
            rich(
                    rich("($(alphabet[sgi+3])) ", font=:bold),
                    [
                        "Mechanism Mode 1",
                        "Mechanism Mode 2"
                    ][sgi]
                )
        end,
        sup! = (ax,tgob,sgi)->begin
            hidexyz(ax)
            # xlims!(ax,-1.0e0,1.2e0)
            # ylims!(ax,-1.2e0,1.0e0)
        end,
    )
    savefig(fig,"prism")
    fig
end

k = TR.get_cables_stiffness(bot.tg)
l = TR.get_cables_len(bot.tg)

Ǩm = TR.build_Ǩm!(bot.tg,q,k)
𝒦m = transpose(Ň)*Ǩm*Ň |> Symmetric
vec𝒦m  = vec(𝒦m)
vecI = vec(Matrix(1.0I,size(𝒦m)))

ᾱ = [1.0,1.0]
f = S*ᾱ 

# equivalent μ
# μ = l .- (f./k)

λ = -inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*f
# @show f,λ
Ǩa = TR.∂Aᵀλ∂q̌(bot.tg,λ)
𝒦ain = transpose(Nin)*Ǩa*Nin
𝒦a = transpose(Ň)*Ǩa*Ň |> Symmetric 
vals_𝒦a,vecs_𝒦a = eigen(𝒦a)
@myshow sort(vals_𝒦a)
@myshow 𝒦a[1:5,1:5]
# @show count((x)->x<0,D_𝒦a)
# @show count((x)->x==0,D_𝒦a)

Ǩg = TR.build_Ǩg!(bot.tg,q,f)

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

σs = 0:1:1600
Vals =  [
    begin
        𝒦 = 𝒦m + σ*reshape(mat𝒦ps*ᾱ,size(𝒦m))
        vals_𝒦, vecs_𝒦 = eigen(𝒦)
        vals_𝒦
    end
    for σ in σs
] |> VectorOfArray

with_theme(theme_pub;
        resolution = (0.6tw,0.2tw),
        figure_padding = (0,fontsize,0,fontsize),
    ) do 
    fig = Figure()
    ax1 = Axis(fig[1,1],
        xlabel = L"\sigma",
        ylabel = L"\rho_{\mathrm{1}}"
    )
    lines!(ax1,σs,Vals[1,:],)
    xlims!(ax1,0,1700)
    ylims!(ax1,-20,60)
    scatter!(
        ax1,
        [σ_max,σ_min],
        [ρ_max,ρ_min]
    )
    text!(ax1,
        [σ_max], [ρ_max], 
        text = [L"\sigma_{\mathrm{max}}"],
        align = (:center,:bottom),
        offset = (0, fontsize/4)
    )
    text!(ax1,
        [σ_min], [ρ_min], 
        text = [L"\sigma_{\mathrm{min}}"],
        align = (:right,:center),
        offset = (-fontsize/2, 0)
    )
    
    ax2 = Axis(fig[1,2],
        xlabel = L"\sigma",
        ylabel = L"\rho"
    )
    for i in 1:6
        lines!(ax2,σs,Vals[i,:],label=latexstring("\\rho_$i"))
    end
    Legend(
        fig[1,3],
        ax2
    )
    xlims!(ax2,0,1700)
    ylims!(ax2,-20,400)
    savefig(fig,"prism_curve")
    fig
end


#-- prism end


#-- two triangles
two = two_tri()
bot = two
plot_traj!(bot;showground=false)
bot.tg.ndof

TR.check_static_equilibrium_output_multipliers(bot.tg)

function make_Ň(tg)    
    (;sysfree,mem2sysndof) = tg.connectivity.indexed
    q = TR.get_q(bot.tg)
    Nin = TR.make_intrinsic_nullspace(tg,q)
    Nin[
        sysfree,
        reduce(vcat,mem2sysndof[2:end])
    ][:,end]
    # I2 = TR.NCF.I2_Int
    # O2 = zero(I2)
    # o2 = O2[:,1]
    # ret = N*[
    #     I2    o2    o2;
    #     o2'    1    0;
    #     I2    o2    o2;
    #     o2'    0    1;
    # ]
    # ret = N*[
    #     I2    o2    o2;
    #     o2'    1    0;
    #     I2    o2    o2;
    #     o2'    0    1;
    # ]
    # ret
end

q = TR.get_q(bot.tg)
q̌ = TR.get_q̌(bot.tg)
Ǎ = TR.make_A(bot.tg)(q)
Ň_ = TR.nullspace(Ǎ)
Ň = TR.modified_gram_schmidt(Ň_)
# Ň = Ň_
# N = TR.make_intrinsic_nullspace(bot.tg,q)

Ň = make_Ň(bot.tg)

rank(Ň)

Ǎ*Ň |> norm

Q̃ = TR.build_Q̃(bot.tg)
L̂ = TR.build_L̂(bot.tg)

# Left hand side
Q̃L̂ = Q̃*L̂

Bᵀ = -Q̃L̂
ℬᵀ = transpose(Ň)*Bᵀ
ℬ = transpose(ℬᵀ)
S,D = TR.static_kinematic_determine(ℬᵀ;atol=1e-14)
S 
D
nk = size(D,2)
ns = size(S,2)

GM.activate!();with_theme(theme_pub;
        resolution = (0.9tw,0.22tw),
        figure_padding = (fontsize,0,0,0),
        fontsize = 6.5 |> pt2px,
        Axis3 = (
            azimuth = -π/2-1e-10,
            elevation = π/2,
        )
    ) do 
    maxS = maximum(abs.(S))
    rtol = 1e-10
    Sbool = S.> maxS*rtol
    S[.!Sbool] .= 0.0
    fig = Figure()
    gd1 = fig[1,1] = GridLayout(;tellheight=false)
    gd2 = fig[1,2:3] = GridLayout(;tellheight=false)
    plot_traj!(
        bot,
        fig = gd1,
        AxisType=Axis3,
        showpoints = false,
        showlabels = false,
        showground = false,
        doslide = false,
        xlims = (-0.5,0.5),
        ylims = (-0.15,0.15),
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[sgi])) ", font=:bold),
                "Inital"
            )
        end,
        sup! = (ax,tgob,sgi)-> begin
            # cables
            hidez(ax)
        end
    )
    rowsize!(gd1,1,Fixed(0.1tw))
    plot_traj!(
        bot,
        fig = gd2,
        AxisType=Axis3,
        gridsize = (2,2),
        showpoints = false,
        showlabels = false,
        showground = false,
        doslide = false,
        showcables = false,
        xlims = (-0.5,0.5),
        ylims = (-0.15,0.15),
        rowgap=0,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                    rich("($(alphabet[sgi+1])) ", font=:bold),
                    [
                        "Self-stress State 1",
                        "Self-stress State 2",
                        "Self-stress State 3",
                        "Self-stress State 4"
                    ][sgi]
                )
        end,
        sup! = (ax,tgob,sgi)-> begin
            # cables
            hidexyz(ax)
            @myshow Sbool[:,sgi]
            linesegs_cables = @lift begin
                get_linesegs_cables($tgob;)[Sbool[:,sgi]]
            end
            linesegments!(ax, 
                linesegs_cables, 
                color = :red, 
                # linewidth = cablewidth
                )
            rcs_by_cables = @lift begin
                (;tensioned) = $tgob.connectivity
                ndim = TR.get_ndim($tgob)
                T = TR.get_numbertype($tgob)
                ret = Vector{MVector{ndim,T}}()
                mapreduce(
                    (scnt)->
                    [(
                        scnt.end1.rbsig.state.rps[scnt.end1.pid].+
                        scnt.end2.rbsig.state.rps[scnt.end2.pid]
                    )./2],
                    vcat,
                    tensioned.connected
                    ;init=ret
                )
            end
            # @show rcs_by_cables
            Stext = [
                    @sprintf "%4.2f"  S[i,sgi] 
                    for i in axes(S,1)
                    if Sbool[i,sgi]
                ]
            @myshow Stext
            # scatter!(
            #     ax,
            #     rcs_by_cables[][Sbool[:,sgi]],
            #     marker = :rect, 
            #     markersize = 12 |> pt2px, 
            #     color = :white
            # )
            text!(
                ax,
                Stext,
                position = rcs_by_cables[][Sbool[:,sgi]],
                fontsize = 5 |> pt2px,
                color = :red,
                align = (:center, :center),
                # offset = (-fontsize/2, 0)
            )
        end
    )
    savefig(fig,"two_tri")
    DataInspector(fig)
    fig
end


k = TR.get_cables_stiffness(bot.tg)

l = TR.get_cables_len(bot.tg)


struct𝒦 = [
    begin
        s = S[:,i]        
        Ǩm = TR.build_Ǩm!(bot.tg,q,100*s)
        𝒦m = transpose(Ň)*Ǩm*Ň 
        # s = S\f
        # @show s
        λ = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*s
        # @show f,λ
        Ǩa = - TR.∂Aᵀλ∂q̌(bot.tg,λ)

        Ǩg = TR.build_Ǩg!(bot.tg,q,s)

        𝒦p = transpose(Ň)*(Ǩg.+Ǩa)*Ň
        @eponymtuple(𝒦m, 𝒦p,)
    end
    for i = 1:ns
] |> StructArray

mat𝒦ps = reduce(hcat,struct𝒦.𝒦p)

mat𝒦ms = reduce(hcat,struct𝒦.𝒦m)

ᾱs = [
   [1.0,0,0,0],
   [0,1.0,0,0],
# #    [0,0,1.0,0],
   [0,0,0,1.0],
   10 .*[0, 0.114892, 0, 0.0748331]
]

@myshow mat𝒦ps
σs = 0:0.1:10
Vals =  [
    begin
        [   
            begin
                𝒦 = mat𝒦ms*ᾱ + σ*mat𝒦ps*ᾱ
                𝒦[1]
            end
            for ᾱ in ᾱs
        ]
    end
    for σ in σs
] |> VectorOfArray

GM.activate!();with_theme(theme_pub;
        resolution = (0.6tw,0.2tw),
        figure_padding = (0,fontsize,0,fontsize),
    ) do 
    fig = Figure()
    ax1 = Axis(fig[1,1],
        xlabel = L"\sigma",
        ylabel = L"\rho_{\mathrm{1}}"
    )
    for i in 1:3
        lines!(ax1,σs,Vals[i,:],label=("Self-stress State $i"))
    end
    lines!(ax1,σs,Vals[4,:],label=("Combined Self-stress State 2 & 3"))

    xlims!(ax1,0,10)
    ylims!(ax1,-0,6)
    # scatter!(
    #     ax1,
    #     [σ_max,σ_min],
    #     [ρ_max,ρ_min]
    # )
    
    # ax2 = Axis(fig[1,2],
    #     xlabel = L"\sigma",
    #     ylabel = L"\rho"
    # )
    # Legend(
    #     fig[1,3],
    #     ax2
    # )
    # xlims!(ax2,0,1700)
    # ylims!(ax2,-20,400)
    Legend(fig[1,2],
        ax1
    )
    savefig(fig,"two_tri_curve")
    fig
end

#-- end two triangles

#-- begin lander
landerbot = lander()
bot = landerbot 
plot_traj!(bot;showground=false)

q = TR.get_q(bot.tg)
q̌ = TR.get_q̌(bot.tg)
Ǎ = TR.make_A(bot.tg)(q)
Ň_ = TR.nullspace(Ǎ)
Ň = TR.modified_gram_schmidt(Ň_)
# Ň = build_Ň(bot.tg)
Q̃ = TR.build_Q̃(bot.tg)
L̂ = TR.build_L̂(bot.tg)

Ǎ*Ň |> norm

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

f = sum(S,dims=2)

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
vec𝒦m = vec(𝒦m)
vecI = vec(Matrix(1.0I,size(𝒦m)))

𝒦g = transpose(Ň)*Ǩg*Ň |> Symmetric 
𝒦p = 𝒦g.+ 𝒦a |> Symmetric 
𝒦 = 𝒦m.+ 𝒦p |> Symmetric

vals_𝒦m,vecs_𝒦m = eigen(𝒦m)
sort(vals_𝒦m)

vals_𝒦,vecs_𝒦 = eigen(𝒦)
sort(vals_𝒦)

ᾱ = ones(60)
A = hcat(
    -Matrix(1.0I,ns,ns),
    ᾱ,
    zero(ᾱ)
)
b = zeros(60)
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
    zeros(60),
    ns+1,
    # result.x
)
σ_min = result_min.x[end]

𝒦_min = 𝒦m + σ_min*reshape(mat𝒦ps*ᾱ,size(𝒦m))
vals_𝒦_min, vecs_𝒦_min = eigen(𝒦_min)
ρ_min = vals_𝒦_min[1]

σs = LinRange(0,0.3,100)
Vals =  [
    begin
        𝒦 = 𝒦m + σ*reshape(mat𝒦ps*ᾱ,size(𝒦m))
        vals_𝒦, vecs_𝒦 = eigen(𝒦)
        vals_𝒦
    end
    for σ in σs
] |> VectorOfArray

with_theme(theme_pub;
        resolution = (0.6tw,0.2tw),
        figure_padding = (0,fontsize,0,fontsize),
    ) do 
    fig = Figure()
    ax1 = Axis(fig[1,1],
        xlabel = L"\sigma",
        ylabel = L"\rho_{\mathrm{1}}"
    )
    lines!(ax1,σs,Vals[1,:],)
    # xlims!(ax1,0,1700)
    # ylims!(ax1,-20,60)
    scatter!(
        ax1,
        [σ_max,σ_min],
        [ρ_max,ρ_min]
    )
    text!(ax1,
        [σ_max], [ρ_max], 
        text = [L"\sigma_{\mathrm{max}}"],
        align = (:center,:bottom),
        offset = (0, fontsize/4)
    )
    text!(ax1,
        [σ_min], [ρ_min], 
        text = [L"\sigma_{\mathrm{min}}"],
        align = (:right,:center),
        offset = (-fontsize/2, 0)
    )
    
    # ax2 = Axis(fig[1,2],
    #     xlabel = L"\sigma",
    #     ylabel = L"\rho"
    # )
    # for i in 1:6
    #     lines!(ax2,σs,Vals[i,:],label=latexstring("\\rho_$i"))
    # end
    # Legend(
    #     fig[1,3],
    #     ax2
    # )
    # xlims!(ax2,0,1700)
    # ylims!(ax2,-20,400)
    # savefig(fig,"lander_curve")
    fig
end


#-- end lander

#-- begin tower

towerbot = tower()
bot = towerbot

@myshow bot.tg.ndof

plot_traj!(bot;showground=false)

function build_Ň(tg)
    (;sysfree,mem2sysfull,mem2sysndof) = tg.connectivity.indexed
    q = TR.get_q(bot.tg)
    Nin = TR.make_intrinsic_nullspace(tg,q)[
        sysfree,
        reduce(vcat,mem2sysndof[begin+1:end])
    ]    
    Nex = zeros(eltype(q),11,5)
    is = 0; js = 0
    Nex[is+4:is+5,js+1:js+2] .= Matrix(1I,2,2)
    is = 5; js = 2
    Nex[is+4:is+6,js+1:js+3] .= Matrix(1I,3,3)

    q_I = q[mem2sysfull[2]]
    ri = @view q_I[1:3]
    u = @view q_I[4:6]
    v,w = TR.NCF.HouseholderOrthogonalization(u)
    @myshow u,v,w
    # R = [u v w;]
    is = 5; js = 0
    Nex[is+1:is+3,js+1:js+2] = -TR.NCF.skew(u)*[v w;]
    Nin*Nex
end

q = TR.get_q(bot.tg)
q̌ = TR.get_q̌(bot.tg)
Ǎ = TR.make_A(bot.tg)(q)
# Ň_ = TR.nullspace(Ǎ)
# Ň = TR.modified_gram_schmidt(Ň_)
Ň = build_Ň(bot.tg)
Q̃ = TR.build_Q̃(bot.tg)
L̂ = TR.build_L̂(bot.tg)

Ǎ*Ň |> norm

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

isis = [8,14,24]
GM.activate!();with_theme(theme_pub;
        resolution = (0.95tw,0.24tw),
        figure_padding = (2fontsize,0,0,0),
        fontsize = 6.5 |> pt2px,
        Axis3 = (
            azimuth = 3.8255306333269843,
            elevation = 0.2026990816987241
        )
    ) do 
    maxS = maximum(abs.(S))
    rtol = 1e-10
    Sbool = S.> maxS*rtol
    S[.!Sbool] .= 0.0
    fig = Figure()
    gd1 = fig[1,1] = GridLayout()
    gd2 = fig[1,2:3+1] = GridLayout()
    botmm = deepcopy(bot)
    plot_traj!(
        bot;
        fig = gd1,
        AxisType=Axis3,
        doslide=false,
        showlabels=false,
        showpoints=false,
        # showcables = false,
        xlims = (-1.5e0,1.5e0),
        ylims = (-1.5e0,1.5e0),
        zlims = (-0.6e0,1.5e0),
        showground = false,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                    rich("($(alphabet[sgi])) ", font=:bold),
                    "Initial"
                )
        end,
        sup! = (ax,tgob,sgi)->begin
            hidex(ax)
            hidey(ax)
            # xlims!(ax,-1.0e0,1.2e0)
            # ylims!(ax,-1.2e0,1.0e0)
        end,
    )
    plot_traj!(
        bot;
        fig = gd2,
        AxisType=Axis3,
        gridsize=(1,3), 
        doslide=false,
        showlabels=false,
        showpoints=false,
        showcables = false,
        xlims = (-1.5e0,1.5e0),
        ylims = (-1.5e0,1.5e0),
        zlims = (-0.6e0,1.5e0),
        showground = false,
        showinit = true,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                    # rich("($(alphabet[sgi+1])) ", font=:bold),
                    "Self-stress State $sgi"
                )
        end,
        sup! = (ax,tgob,sgi_input)-> begin
            sgi = isis[sgi_input]
            # cables
            ax.azimuth = 4.73553063332698
            ax.elevation = 0.18269908169872395
            # azimuth = 4.665530633326984
            # elevation = 0.16269908169872424
            hidexyz(ax)
            @myshow Sbool[:,sgi]
            linesegs_cables = @lift begin
                get_linesegs_cables($tgob;)[Sbool[:,sgi]]
            end
            linesegments!(ax, 
                linesegs_cables, 
                color = :red, 
                # linewidth = cablewidth
                )
            rcs_by_cables = @lift begin
                (;tensioned) = $tgob.connectivity
                ndim = TR.get_ndim($tgob)
                T = TR.get_numbertype($tgob)
                ret = Vector{MVector{ndim,T}}()
                mapreduce(
                    (scnt)->
                    [(
                        scnt.end1.rbsig.state.rps[scnt.end1.pid].+
                        scnt.end2.rbsig.state.rps[scnt.end2.pid]
                    )./2],
                    vcat,
                    tensioned.connected
                    ;init=ret
                )
            end
            # @show rcs_by_cables
            Stext = [
                    @sprintf "%4.2f"  S[i,sgi] 
                    for i in axes(S,1)
                    if Sbool[i,sgi]
                ]
            @myshow Stext
            # scatter!(
            #     ax,
            #     rcs_by_cables[][Sbool[:,sgi]],
            #     marker = :rect, 
            #     markersize = 12 |> pt2px, 
            #     color = :white
            # )
            text!(
                ax,
                Stext,
                position = rcs_by_cables[][Sbool[:,sgi]],
                fontsize = 5 |> pt2px,
                color = :red,
                align = (:center, :center),
                # offset = (-fontsize/2, 0)
            )
        end
    )
    savefig(fig,"tower")
    fig
end

ᾱ = zeros(ns)
ᾱ[isis] .= 1.0
f = S*ᾱ

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
𝒦m = transpose(Ň)*Ǩm*Ň |> Symmetric
vec𝒦m = vec(𝒦m)
vecI = vec(Matrix(1.0I,size(𝒦m)))

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
        vec𝒦pi = vec(𝒦pi)
    end
    for i = 1:ns
]

mat𝒦ps = reduce(hcat,vec𝒦ps)

A = hcat(
    -Matrix(1.0I,ns,ns),
    ᾱ,
    zero(ᾱ)
)
b = zeros(ns)
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
    zeros(ns),
    ns+1,
    result_max.x[1:end-1]
)
σ_min = result_min.x[end]

𝒦_min = 𝒦m + σ_min*reshape(mat𝒦ps*ᾱ,size(𝒦m))
vals_𝒦_min, vecs_𝒦_min = eigen(𝒦_min)
ρ_min = vals_𝒦_min[1]

σs = LinRange(-50,100,100)
Vals =  [
    begin
        𝒦 = 𝒦m + σ*reshape(mat𝒦ps*ᾱ,size(𝒦m))
        vals_𝒦, vecs_𝒦 = eigen(𝒦)
        vals_𝒦
    end
    for σ in σs
] |> VectorOfArray

with_theme(theme_pub;
        resolution = (0.35tw,0.2tw),
        figure_padding = (0,fontsize,0,fontsize),
    ) do 
    fig = Figure()
    ax1 = Axis(fig[1,1],
        xlabel = L"\sigma",
        ylabel = L"\rho_{\mathrm{1}}"
    )
    lines!(ax1,σs,Vals[1,:],)
    ylims!(ax1,-30,140)
    xlims!(ax1,-50,100)
    scatter!(
        ax1,
        [σ_max,σ_min],
        [ρ_max,ρ_min]
    )
    text!(ax1,
        [σ_max], [ρ_max], 
        text = [L"\sigma_{\mathrm{max}}"],
        align = (:center,:bottom),
        offset = (0, fontsize/4)
    )
    text!(ax1,
        [σ_min], [ρ_min], 
        text = [L"\sigma_{\mathrm{min}}"],
        align = (:right,:center),
        offset = (-fontsize/2, 0)
    )
    
    # ax2 = Axis(fig[1,2],
    #     xlabel = L"\sigma",
    #     ylabel = L"\rho"
    # )
    # for i in 1:5
    #     lines!(ax2,σs,Vals[i,:],label=latexstring("\\rho_$i"))
    # end
    # Legend(
    #     fig[1,3],
    #     ax2
    # )
    # xlims!(ax2,0,1700)
    # ylims!(ax2,-20,400)
    savefig(fig,"tower_curve")
    fig
end


#-- end tower

#-- T bars
tbbot = Tbars()
bot = tbbot
@myshow bot.tg.ndof

plot_traj!(bot;showground=false)

TR.check_static_equilibrium_output_multipliers(bot.tg)

function make_Ň(tg)    
    (;sysfree,mem2sysfull,mem2sysndof) = tg.connectivity.indexed
    q = TR.get_q(bot.tg)
    Nin = TR.make_intrinsic_nullspace(tg,q)[
        sysfree,
        reduce(vcat,mem2sysndof[2:end])
    ]
    I3 = TR.NCF.I3_Int
    O3 = zero(I3)
    o3 = O3[:,1]
    qbar = q[mem2sysfull[4]]
    ri = @view qbar[1:3]
    u = @view qbar[4:6]
    v,w = TR.NCF.HouseholderOrthogonalization(u)
    @myshow v,w
    x,_,_ = ri
    Nslider1 = [
        o3;
        o3;;
    ] |> sparse
    Nslider2 = [
         0;
        -x;
         0;
        o3;;
    ] |> sparse
    Nbar = [
        o3;
        0;
        1;;
    ] |> sparse
    Nex = vcat(Nslider1,Nslider2,Nbar)
    Nin*Nex
end

q = TR.get_q(bot.tg)
q̌ = TR.get_q̌(bot.tg)
Ǎ = TR.make_A(bot.tg)(q)
Ň_ = TR.nullspace(Ǎ)
Ň = TR.modified_gram_schmidt(Ň_)

Ň = make_Ň(bot.tg)

# todo construct null space
rank(Ň)

Ǎ*Ň |> norm

Q̃ = TR.build_Q̃(bot.tg)
L̂ = TR.build_L̂(bot.tg)

# Left hand side
Q̃L̂ = Q̃*L̂

Bᵀ = -Q̃L̂
ℬᵀ = transpose(Ň)*Bᵀ
ℬ = transpose(ℬᵀ)
S,D = TR.static_kinematic_determine(ℬᵀ;atol=1e-14)
S 
D
nk = size(D,2)
ns = size(S,2)

k = TR.get_cables_stiffness(bot.tg)

l = TR.get_cables_len(bot.tg)

struct𝒦 = [
    begin
        s = S[:,i]        
        Ǩm = TR.build_Ǩm!(bot.tg,q,100*s)
        𝒦m = transpose(Ň)*Ǩm*Ň 
        # s = S\f
        # @show s
        λ = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*s
        # @show f,λ
        Ǩa = - TR.∂Aᵀλ∂q̌(bot.tg,λ)

        Ǩg = TR.build_Ǩg!(bot.tg,q,s)

        𝒦g = transpose(Ň)*(Ǩg)*Ň

        𝒦a = transpose(Ň)*(Ǩa)*Ň

        𝒦p = 𝒦g .+ 𝒦a
        @eponymtuple(𝒦m, 𝒦g, 𝒦a, 𝒦p,)
    end
    for i = 1:ns
] |> StructArray


mat𝒦ms = reduce(hcat,struct𝒦.𝒦m)
mat𝒦gs = reduce(hcat,struct𝒦.𝒦g)
mat𝒦as = reduce(hcat,struct𝒦.𝒦a)
mat𝒦ps = reduce(hcat,struct𝒦.𝒦p)

GM.activate!();with_theme(theme_pub;
        resolution = (1tw,0.32tw),
        figure_padding = (fontsize,fontsize,0,0),
        Axis3 = (
            azimuth = -π/2-1e-10,
            elevation = π/2,
        )
    ) do 
    maxS = maximum(abs.(S))
    rtol = 1e-10
    Sbool = S.> maxS*rtol
    S[.!Sbool] .= 0.0
    fig = Figure()
    gd0 = fig[1,1] = GridLayout(;tellheight=false)
    gd1 = fig[1,2] = GridLayout(;tellheight=false)
    gd2 = fig[1,3:4] = GridLayout(;tellheight=false)
    rowsize!(gd0,1,Fixed(0.15tw))
    rowsize!(gd1,1,Fixed(0.15tw))
    plot_traj!(
        bot,
        fig = gd0,
        AxisType = Axis3,
        showpoints = false,
        showlabels = false,
        showground = false,
        doslide = false,
        xlims = (-2.5,1.5),
        ylims = (-1.5,1.5),
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[sgi])) ", font=:bold),
                "Inital"
            )
        end,
        sup! = (ax,tgob,sgi)-> begin
            # cables
            hidez(ax)
        end
    )
    botmm = deepcopy(bot)
    δq̌ = Ň*[1.0]
    scaling=0.2
    push!(botmm.traj,deepcopy(botmm.traj[end]))
    botmm.traj.t[end] = 1
    ratio = norm(δq̌)/norm(q̌)
    botmm.traj.q̌[end] .= q̌ .+ scaling.*δq̌/ratio
    plot_traj!(
        botmm,
        fig = gd1,
        AxisType = Axis3,
        atsteps=[2],
        showpoints = false,
        showlabels = false,
        showground = false,
        doslide = false,
        xlims = (-2.5,1.5),
        ylims = (-1.5,1.5),
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[sgi+1])) ", font=:bold),
                "Mode 1"
            )
        end,
        sup! = (ax,tgob,sgi)-> begin
            # cables
            hidez(ax)
        end
    )
    plot_traj!(
        bot,
        fig = gd2,
        AxisType=Axis3,
        gridsize = (2,2),
        showpoints = false,
        showlabels = false,
        showground = false,
        doslide = false,
        showcables = false,
        xlims = (-2.5,1.5),
        ylims = (-1.5,1.5),
        rowgap=0,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                    rich("($(alphabet[sgi+2])) ", font=:bold),
                    [
                        "Self-stress State 1",
                        "Self-stress State 2",
                        "Self-stress State 3",
                        "Self-stress State 4"
                    ][sgi]
                )
        end,
        sup! = (ax,tgob,sgi)-> begin
            # cables
            hidexyz(ax)
            @myshow Sbool[:,sgi]
            linesegs_cables = @lift begin
                get_linesegs_cables($tgob;)[Sbool[:,sgi]]
            end
            linesegments!(ax, 
                linesegs_cables, 
                color = :red, 
                # linewidth = cablewidth
                )
            rcs_by_cables = @lift begin
                (;tensioned) = $tgob.connectivity
                ndim = TR.get_ndim($tgob)
                T = TR.get_numbertype($tgob)
                ret = Vector{MVector{ndim,T}}()
                mapreduce(
                    (scnt)->
                    [(
                        scnt.end1.rbsig.state.rps[scnt.end1.pid].+
                        scnt.end2.rbsig.state.rps[scnt.end2.pid]
                    )./2],
                    vcat,
                    tensioned.connected
                    ;init=ret
                )
            end
            # @show rcs_by_cables
            Stext = [
                    @sprintf "%4.2f"  S[i,sgi] 
                    for i in axes(S,1)
                    if Sbool[i,sgi]
                ]
            @myshow Stext
            # scatter!(
            #     ax,
            #     rcs_by_cables[][Sbool[:,sgi]],
            #     marker = :rect, 
            #     markersize = 12 |> pt2px, 
            #     color = :white
            # )
            text!(
                ax,
                Stext,
                position = rcs_by_cables[][Sbool[:,sgi]],
                fontsize = 5 |> pt2px,
                color = :red,
                align = (:center, :center),
                # offset = (-fontsize/2, 0)
            )
        end
    )
    savefig(fig,"Tbars")
    DataInspector(fig)
    fig
end

ᾱs = [
   [1.0,1.0,1.0,0],
   [1.0,1.0,0,1.0],
   [1.0,1.0,1.0,1.0],
]

@myshow mat𝒦ps
σs = 0:0.1:10
Vals =  [
    begin
        [   
            begin
                𝒦 = mat𝒦ms*ᾱ + σ*mat𝒦ps*ᾱ
                𝒦[1]
            end
            for ᾱ in ᾱs
        ]
    end
    for σ in σs
] |> VectorOfArray

GM.activate!();with_theme(theme_pub;
        resolution = (0.55tw,0.2tw),
        figure_padding = (0,fontsize,0,fontsize),
    ) do 
    fig = Figure()
    ax1 = Axis(fig[1,1],
        xlabel = L"\sigma",
        ylabel = L"\rho_{\mathrm{1}}"
    )
    for i in 1:3
        lines!(ax1,σs,Vals[i,:],label=("Combined Self-stress State $i"))
    end
    xlims!(ax1,0,10)

    # ylims!(ax1,-0,6)
    # scatter!(
    #     ax1,
    #     [σ_max,σ_min],
    #     [ρ_max,ρ_min]
    # )
    
    # ax2 = Axis(fig[1,2],
    #     xlabel = L"\sigma",
    #     ylabel = L"\rho"
    # )
    # Legend(
    #     fig[1,3],
    #     ax2
    # )
    # xlims!(ax2,0,1700)
    # ylims!(ax2,-20,400)
    Legend(fig[1,2],
        ax1
    )
    savefig(fig,"Tbars_curve")
    fig
end

#-- end Tbars 

#-- begin uni bot
unibot = uni(0.0;
        μ = 0.9,
        e = 0.0,
        z0 = 0.2,
        isbody = true,
)
bot = unibot
plot_traj!(bot;showground=false)

@myshow bot.tg.ndof

plot_rigid(TR.get_rigidbodies(bot)[2])


TR.check_static_equilibrium_output_multipliers(bot.tg)

q = TR.get_q(bot.tg)
q̌ = TR.get_q̌(bot.tg)
Ǎ = TR.make_A(bot.tg)(q)
Ň_ = TR.nullspace(Ǎ)
Ň = TR.modified_gram_schmidt(Ň_)
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
S,D = TR.static_kinematic_determine(ℬᵀ;atol=1e-14)
S 
D
ns = size(S,2)
nk = size(D,2)
δq̌ = [Ň*D[:,i] for i in axes(D,2)]
scaling=0.1
botvis = deepcopy(bot)
for i = 1:nk
    push!(botvis.traj,deepcopy(botvis.traj[end]))
    botvis.traj.t[end] = i
    δq̌i = δq̌[i]
    ratio = norm(δq̌i)/norm(q̌)
    botvis.traj.q̌[end] .= q̌ .+ scaling.*δq̌i/ratio
end

plot_traj!(botvis;showground=false)

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

# ULTRA Spine 2D
spine2d = make_spine(2)
bot = spine2d

#-- planar_parallel
pp = planar_parallel()
bot = pp


