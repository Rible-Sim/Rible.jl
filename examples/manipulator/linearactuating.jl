using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
# using BenchmarkTools
import PyPlot; const plt = PyPlot
plt.pygui(true)
using LaTeXStrings
# using NLsolve
# using DifferentialEquations
# using Sundials
# using DASKR
using Makie
AbstractPlotting.__init__()
# using MakieLayout
using Revise

using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

cd("examples/manipulator")
include("man_define.jl")
include("man_plotting.jl")
include("../analysis.jl")

function simulate_linearactuating(ndof=6;k,c,target_t=20.0,tend=100.0)
    manipulator = man_ndof(ndof,k=k,c=c)
    Y = build_Y(manipulator)
    q0,q̇0,λ0 = TR.get_initial(manipulator)

    refman = man_ndof(ndof,θ=-π/12,k=k) # reference

    refλ0,_,a = TR.inverse(manipulator,deepcopy(refman),Y)
    refq0,_,_ = TR.get_initial(refman)
    # TR.actuate!(manipulator,a)
    # TR.reset_forces!(manipulator)
    # TR.distribute_q_to_rbs!(manipulator,refq0)
    # TR.update_strings_apply_forces!(manipulator)
    # @show transpose(TR.build_A(manipulator)(refq0))*refλ0 - TR.build_Q̃(manipulator)*TR.fvector(manipulator)
    #
    function linearactuate(tgstruct,q0,target_u,target_t)

        M = TR.build_massmatrix(tgstruct)
        Φ = TR.build_Φ(tgstruct,q0)
        A = TR.build_A(tgstruct)

        function F!(F,q,q̇,t)
            if t < target_t
                TR.actuate!(tgstruct,t/target_t*target_u)
            end
            TR.reset_forces!(tgstruct)
            TR.distribute_q_to_rbs!(tgstruct,q,q̇)
            TR.update_strings_apply_forces!(tgstruct)
            TR.assemble_forces!(F,tgstruct)
        end

        M,Φ,A,F!,nothing

        #A,Φ,∂T∂q̇!,F!,M!,nothing
    end



    dt = 0.01 # Same dt used for PID AND Dynamics solver
    prob = TS.DyProblem(linearactuate(manipulator,q0,a,target_t),q0,q̇0,λ0,(0.0,tend))

    sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-14)
    manipulator, sol
end

# bars_and_strings_segs(man_linear)
k = 1.e2
c = 200.0
target_t = 10
tend = 20.0

man_linear, sol_linear = simulate_linearactuating(;k,c,target_t,tend)
# scene,_ = plotstructure(man_linear)

plotstructure(man_linear,sol_linear,sliderplot)
# plotstructure(man_linear,sol_linear,recordplot;filename="linearactuating.mp4")

tstops = [0,5,10,15,20]
# tstops = [0,24,48,72,96,120]

function pyplot2(man_linear, sol_linear, tstops)

    steps = [findfirst((x)->x>i,sol_linear.ts) for i = tstops]

    fig,axs_raw = plt.subplots(1,1,figsize=(9,6))
    grid = (2,6)
    axs1 = [plt.subplot2grid(grid,(0,2i-2),1,2,fig=fig) for i = 1:3]
    axs2 = [plt.subplot2grid(grid,(1,2i-1),1,2,fig=fig) for i = 1:2]
    axs = vcat(axs1,axs2)
    for (i,step) in enumerate(steps)
        TR.distribute_q_to_rbs!(man_linear,sol_linear.qs[step])
        bars_segs,strings_segs = bars_and_strings_segs(man_linear)
        ax = axs[i]
        ax.add_collection(bars_segs)
        ax.add_collection(strings_segs)
        ax.set_ylim(-1.50,0.20)
        ax.set_xlim(-0.20,1.50)
        ax.set_aspect("equal")
        ax.set_title("t=$(tstops[i])s")
        ax.grid(true)
        # ax.set_yticks(collect(0:-50:-150))
        ax.set_xticks([0,0.5,1])
        ax.set_yticks([-1.5,-1.0,-0.5,0])
        if i ∈[0]
            ax.set_xticklabels([])
            ax.xaxis.label.set_visible(false)
        end
        if i ∈[2,3,5]
            ax.set_yticklabels([])
            ax.yaxis.label.set_visible(false)
        end
    end
    fig
end

fig = pyplot2(man_linear, sol_linear, tstops)
fig.savefig("linearactuating3x2_k=$(k)_c=$(c)_t=$(target_t).png",dpi=300,bbox_inches="tight")

kes,epes,gpes,es,es_err = analyse_energy(man_linear,sol_linear;gravity=false,elasticity=true,factor=1)

function plot_energy()
    fig,ax = plt.subplots(1,1,figsize=(7,5))
    grid = (2,4)
    axs = [plt.subplot2grid(grid,(0,2i-2),1,2,fig=fig) for i = 1:2]
    push!(axs,plt.subplot2grid(grid,(1,1),1,2,fig=fig))
    xmin = sol_linear.ts[1]
    xmax = sol_linear.ts[end]
    axs[1].plot(sol_linear.ts,es)
    axs[1].set_xlim(xmin,xmax)
    axs[1].set_ylim(0.5,0.8)
    axs[1].set_xlabel("Time (s)")
    axs[1].set_title("Total Energy")
    axs[2].plot(sol_linear.ts,epes)
    axs[2].set_xlim(xmin,xmax)
    axs[2].set_ylim(0.5,0.8)
    axs[2].set_xlabel("Time (s)")
    axs[2].set_title("Elastic Potential Energy")
    axs[3].plot(sol_linear.ts,kes)
    axs[3].set_xlim(xmin,xmax)
    axs[3].set_ylim(-0.002,0.025)
    axs[3].set_xlabel("Time (s)")
    axs[3].set_title("Kinetic Energy")
    fig.tight_layout()
    fig
end
fig = plot_energy()
fig.savefig("energy3x1_k=$(k)_c=$(c)_t=$(target_t).png",dpi=300,bbox_inches="tight")
