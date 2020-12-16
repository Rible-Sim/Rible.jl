using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using Makie
AbstractPlotting.__init__()
import PyPlot; const plt = PyPlot
plt.pygui(true)
using LaTeXStrings
# using NLsolve
# using DifferentialEquations
using Printf
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot
cd("examples/manipulator")
include("man_define.jl")
include("man_plotting.jl")

k = 1.e2; c = 0.0
function inversekinematics(θs;k,c)
    n = length(θs)
    ndof = 6
    manipulator = man_ndof(ndof;k,c)
    as = Vector{Vector{Float64}}()

    Y = build_Y(manipulator)


    for θ in θs
        refman = man_ndof(ndof;θ,k,c) # reference
        refqi,_,_ = TR.get_initial(refman)
        refλi,_,a= TR.inverse(manipulator,refman,Y)
        man_act = deepcopy(manipulator)
        TR.actuate!(man_act,a)
        rls = [s.state.restlen for s in man_act.strings]
        push!(as,rls)
    end
    manipulator, as
end

α = 0.8
n = 10
θs = [-i*α/n for i = 0:n]

man_inv, as = inversekinematics(θs; k, c)
# man_inv_scene,_ = plotstructure(man_inv)
# man_inv_scene

fig,ax = plt.subplots(1,1,num = 1,figsize=(4,3))
x = abs.(θs)
y = [a[1] for a in as]
z = [a[2] for a in as]
ax.plot(x,y,marker="o")
ax.plot(x,z,marker="o")
ax.set_xlabel(L"$\theta$ (Rad)")
ax.set_ylabel("Rest length (m)")
# ax.set_xlim(0,1.0)
# ax.set_ylim(0,9.0)
ax.grid(true)
fig.savefig("man_inverse.png",dpi=300,bbox_inches="tight")
fig.savefig("man_inverse.eps",dpi=300,bbox_inches="tight")

function pyplot2(θs;k,c,ndof = 6)
    fig,axs_raw = plt.subplots(2,3,figsize=(9,6))
    axs = permutedims(axs_raw,[2,1])
    for (i,θ) in enumerate(θs[1:2:end])
        refman = man_ndof(ndof;θ,k,c)
        bars_segs,strings_segs = bars_and_strings_segs(refman)
        ax = axs[i]
        ax.add_collection(bars_segs)
        ax.add_collection(strings_segs)
        ax.set_ylim(-1.10,0.20)
        ax.set_xlim(-0.20,1.30)
        ax.set_aspect("equal")
        ax.set_title(latexstring("\\theta=",@sprintf "%3.2f" θ))
        ax.grid(true)
        # ax.set_yticks(collect(0:-50:-150))
        if i ∈[1,2,3]
            ax.set_xticklabels([])
            ax.xaxis.label.set_visible(false)
        end
        if i ∈[2,3,5,6]
            ax.set_yticklabels([])
            ax.yaxis.label.set_visible(false)
        end
    end
    fig.tight_layout()
    fig
end
fig = pyplot2(θs;k,c,ndof = 6)
fig.savefig("inv2x3.png",dpi=300,bbox_inches="tight")


ndof = 6
k = 1.e2
manipulator = man_ndof(ndof; k, c)
ncoords = manipulator.ncoords
function eigenanalysis(θs; k, c)
    ndof = 6
    manipulator = man_ndof(ndof; k, c)

    Y = build_Y(manipulator)

    ωs = Vector{Vector{Float64}}()
    Zs = Vector{Matrix{Float64}}()

    for θ in θs
        refman = man_ndof(ndof;θ,k,c) # reference
        refqi,_,_ = TR.get_initial(refman)
        refλi,_,a= TR.inverse(manipulator,refman,Y)
        actmani = deepcopy(manipulator)
        TR.actuate!(actmani,a)
        ωi,Zi = TR.undamped_eigen!(actmani,refqi,refλi)
        push!(ωs,ωi)
        push!(Zs,Zi)
    end
    ωs,Zs
end

ωs,Zs = eigenanalysis(θs; k, c)
ys = [[ω[i] for ω in ωs] for i = 1:length(ωs[1])]

function plotfrequency(x,ys)
    fig,ax = plt.subplots(1,1,figsize=(5.5,3))
    for (i,y) in enumerate(ys)
        ax.plot(x,y,marker="o",fillstyle="none",label="Mode $i")
    end
    ax.set_xlabel(L"\theta\ (\mathrm{rad})")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_xlim(x[1],x[end])
    # ax.set_ylim(0,8)
    ax.grid(true)
    ax.legend(loc="upper right", bbox_to_anchor=(1.4, 1.0))
    fig.tight_layout()
    fig
end
fig = plotfrequency(abs.(θs),ys)
fig.savefig("man_frequency_k=$k.png",dpi=300,bbox_inches="tight")

function plotfrequency_ratio(x,ys)
    fig,ax = plt.subplots(1,1,figsize=(5.5,3))
    ys_ratio = [[y[i]/y[1] for i = eachindex(y)] for y in ys]
    for (i,y) in enumerate(ys_ratio)
        ax.plot(x,y,marker="o",fillstyle="none",label="Mode $i")
    end
    ax.set_xlabel(L"\theta\ (\mathrm{rad})")
    ax.set_ylabel("Ratio")
    ax.set_xlim(x[1],x[end])
    # ax.set_ylim(0,8)
    ax.grid(true)
    ax.legend(loc="upper right", bbox_to_anchor=(1.4, 1.0))
    fig.tight_layout()
    fig
end
fig = plotfrequency_ratio(θs,ys)
fig.savefig("man_frequency_ratio_k=$k.png",dpi=300,bbox_inches="tight")

function pyplot_mode(tgstruct,θ,Z,xlim,ylim)
    function plot2ax!(ax,tgstruct,q;ref=false)
        TR.distribute_q_to_rbs!(tgstruct,q)
        bars_segs,strings_segs = bars_and_strings_segs(tgstruct;ref)
        ax.add_collection(bars_segs)
        ax.add_collection(strings_segs)
        ax.set_xlim(xlim[1],xlim[2])
        ax.set_ylim(ylim[1],ylim[2])
        ax.set_aspect("equal")
        ax.grid(true)
    end
    fig = plt.figure(figsize=(9,5))
    grid = (1,3)
    q0,q̇0 = TR.get_q(tgstruct)
    # ax = plt.subplot2grid(grid,(0,1),1,1,fig=fig)
    # ax.set_title("Initial pose, "*latexstring("\\theta=",@sprintf "%3.2f" θ))
    # plot2ax!(ax,tgstruct,q0;ref=true)
    for i = 1:3
        ax = plt.subplot2grid(grid,(0,i-1),1,1,fig=fig)
        qmode = q0+0.2Z[:,i]
        plot2ax!(ax,tgstruct,q0;ref=true)
        plot2ax!(ax,tgstruct,qmode)
        ax.set_title("Mode $i")
        if i >= 2
            ax.set_yticklabels([])
            ax.yaxis.label.set_visible(false)
        end
    end
    fig.tight_layout()
    fig
end

manipulator = man_ndof(ndof;θ = θs[1], k, c)
fig = pyplot_mode(manipulator,θs[1],Zs[1],(-0.20,1.50),(-1.00,0.50))
fig.savefig("man_pose1_mode_k=$k.png",dpi=300,bbox_inches="tight")

manipulator = man_ndof(ndof;θ = θs[5], k, c)
fig = pyplot_mode(manipulator,θs[5],Zs[5],(-0.20,1.50),(-1.00,0.50))
fig.savefig("man_pose5_mode_k=$k.png",dpi=300,bbox_inches="tight")
