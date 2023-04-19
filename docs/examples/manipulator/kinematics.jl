using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using Makie
import PyPlot as plt
plt.pygui(true)
using LaTeXStrings
using RecursiveArrayTools
# using NLsolve
# using DifferentialEquations
using Printf
using Revise
import TensegrityRobots as TR
cd("examples/manipulator")
include("man_define.jl")
include("man_plotting.jl")
include("../plot_helpers.jl")
set_pyplot()
k = 1.e2; c = 0.0
function inversekinematics(θs;k,c)
    n = length(θs)
    ndof = 6
    manipulator = man_ndof(ndof;k,c)
    as = Vector{Vector{Float64}}()

    Y = build_Y(manipulator)


    for θ in θs
        refman = man_ndof(ndof;θ,k,c) # reference
        refqi,_,_ = TR.get_initial(refman.tg)
        refλi,_,a= TR.inverse(manipulator,refman,Y)
        man_act = deepcopy(manipulator)
        TR.actuate!(man_act,a)
        rls = [s.state.restlen for s in man_act.tg.cables]
        push!(as,rls)
    end
    manipulator, as
end

α = 0.8
n = 10
θs = [-i*α/n for i = 0:n]

man_inv, as = inversekinematics(θs; k, c)

man_inv = man_ndof(6;θ=-0.2,k,c)
Y = build_Y(man_inv)
q,_ = TR.get_q(man_inv.tg)
λ,u,a= TR.inverse(man_inv,deepcopy(man_inv),Y)
TR.actuate!(man_inv,a)
TR.check_static_equilibrium(man_inv.tg,q,λ)



man_inv_scene = plotstructure(man_inv)

fig,ax = plt.subplots(1,1,figsize=(Single_width,2.3))
x = abs.(θs)
y = [a[1] for a in as]
z = [a[2] for a in as]
ax.plot(x,y,marker=markers.xs[1],label=L"u_{2i-1}")
ax.plot(x,z,marker=markers.xs[2],label=L"u_{2i}")
ax.set_xlabel(L"$\theta$ (Rad)")
ax.set_ylabel("Rest length (m)")
# ax.set_xlim(0,1.0)
ax.set_ylim(0.20,0.36)
ax.legend()
ax.grid(true)
fig.savefig("man_inverse.png",dpi=300,bbox_inches="tight")
fig.savefig("man_inverse.eps",dpi=300,bbox_inches="tight")

function pyplot2(θs;k,c,ndof = 6)
    fig,axs_raw = plt.subplots(2,3,figsize=(Full_width,4.5))
    axs = permutedims(axs_raw,[2,1])

    for (i,θ) in enumerate(θs[1:2:end])
        refman = man_ndof(ndof;θ,k,c)
        bars_segs,cables_segs = bars_and_cables_segs(refman.tg)
        ax = axs[i]
        ax.add_collection(bars_segs)
        ax.add_collection(cables_segs)
        ax.set_ylim(-1.00,0.20)
        ax.set_xlim(-0.25,1.35)
        ax.set_aspect("equal")
        ax.set_title("($(alphabet[i])) "*
                    latexstring("\\theta")*@sprintf("=%3.2f",abs(θ)),
                    y = -0.3)
        ax.grid(true)
        # ax.set_yticks(collect(0:-50:-150))
        # if i ∈[1,2,3]
        #     ax.set_xticklabels([])
        #     ax.xaxis.label.set_visible(false)
        # end
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
ncoords = manipulator.tg.ncoords
function eigenanalysis(θs; k, c)
    ndof = 6
    manipulator = man_ndof(ndof; k, c)

    Y = build_Y(manipulator)

    ωs = VectorOfArray(Vector{Vector{Float64}}())
    Zs = VectorOfArray(Vector{Matrix{Float64}}())
    as = VectorOfArray(Vector{Vector{Float64}}())

    for θ in θs
        refman = man_ndof(ndof;θ,k,c) # reference
        refqi,_,_ = TR.get_initial(refman.tg)
        refλi,_,ai= TR.inverse(manipulator,refman,Y)
        actmani = deepcopy(manipulator)
        TR.actuate!(actmani,ai)
        ωi,Zi = TR.undamped_eigen!(actmani.tg,refqi,refλi)
        push!(ωs,ωi)
        push!(Zs.u,Zi)
        push!(as.u,ai)
    end
    ωs,Zs,as
end

ωs,Zs,as = eigenanalysis(θs; k, c)
ωs2,_,_ = eigenanalysis(θs; k=10k, c)

function plotfrequency(x,ωs,ωs2)
    fig,_ = plt.subplots(1,2,figsize=(Full_width,6))
    grid = (2,4)
    axs = [plt.subplot2grid(grid,(0,2i-2),1,2,fig=fig) for i = 1:2]
    push!(axs,plt.subplot2grid(grid,(1,1),1,2,fig=fig))

    for i in 1:length(ωs[1])
        axs[1].plot(x,ωs[i,:],marker=markers.xs[i],label="Mode $i")
    end
    axs[1].set_xlabel(L"\theta"*" (rad)")
    axs[1].set_ylabel("Frequency (Hz)")
    axs[1].set_xlim(x[1],x[end])
    # axs[1].set_ylim(0,8)
    axs[1].grid(true)
    axs[1].set_title("(a)", y = -0.35)

    for i in 1:length(ωs2[1])
        axs[2].plot(x,ωs2[i,:],marker=markers.xs[i],label="Mode $i")
    end
    axs[2].set_xlabel(L"\theta"*" (rad)")
    axs[2].set_ylabel("Frequency (Hz)")
    axs[2].set_xlim(x[1],x[end])
    # axs[2].set_ylim(0,8)
    axs[2].grid(true)
    axs[2].set_title("(b)", y = -0.35)

    ys_ratio = [ωs[i,:]/ωs[i,1] for i = 1:length(ωs[1])]
    for (i,y) in enumerate(ys_ratio)
        axs[3].plot(x,y,marker=markers.xs[i],label="Mode $i")
    end
    axs[3].set_xlabel(L"\theta"*" (rad)")
    axs[3].set_ylabel("Ratio")
    axs[3].set_xlim(x[1],x[end])
    # axs[3].set_ylim(0,8)
    axs[3].grid(true)
    axs[3].legend(loc="upper left", bbox_to_anchor=(1.1, 1.04))
    axs[3].set_title("(c)", y = -0.35)

    fig.tight_layout()
    plt.subplots_adjust(wspace=1)
    fig
end
fig = plotfrequency(abs.(θs),ωs,ωs2)
# plotfrequency(abs.(θs),as)
fig.savefig("man_frequency_k=$k.png",dpi=300,bbox_inches="tight")



function pyplot_mode(tgstruct,θ,Z,xlim,ylim)
    function plot2ax!(ax,tgstruct,q;ref=false)
        TR.distribute_q_to_rbs!(tgstruct,q)
        bars_segs,cables_segs = bars_and_cables_segs(tgstruct;ref)
        ax.add_collection(bars_segs)
        ax.add_collection(cables_segs)
        ax.set_xlim(xlim[1],xlim[2])
        ax.set_ylim(ylim[1],ylim[2])
        ax.set_aspect("equal")
        ax.grid(true)
    end
    fig = plt.figure(figsize=(Full_width,2.5))
    grid = (1,3)
    q0,_ = TR.get_q(tgstruct)

    for i = 1:3
        ax = plt.subplot2grid(grid,(0,i-1),1,1,fig=fig)
        qmode = q0+0.1Z[:,i]
        modestruct = deepcopy(tgstruct)
        plot2ax!(ax,modestruct,q0;ref=true)
        plot2ax!(ax,modestruct,qmode)
        ax.set_title("($(alphabet[i])) Mode $i",y = -0.3)
        if i >= 2
            ax.set_yticklabels([])
            ax.yaxis.label.set_visible(false)
        end
    end
    fig.tight_layout()
    fig
end

manipulator = man_ndof(ndof;θ = θs[1], k, c)
fig = pyplot_mode(manipulator,θs[1],Zs[1],(-0.20,1.50),(-1.00,0.25))
fig.savefig("man_pose1_mode_k=$k.png",dpi=300,bbox_inches="tight")

manipulator = man_ndof(ndof;θ = θs[6], k, c)
fig = pyplot_mode(manipulator,θs[6],Zs[6],(-0.20,1.50),(-1.00,0.25))
fig.savefig("man_pose6_mode_k=$k.png",dpi=300,bbox_inches="tight")
