using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using BenchmarkTools
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


function simulate_linearactuating(ndof=6;k,c)
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
    prob = TS.DyProblem(linearactuate(manipulator,q0,a,20.0),q0,q̇0,λ0,(0.0,100.0))

    sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-9,verbose=true)
    manipulator, sol
end

# bars_and_strings_segs(man_linear)
k = 1.e4; c = 1.e5
man_linear, sol_linear = simulate_linearactuating(;k,c)
plotstructure(man_linear,sol_linear,sliderplot)
# plotstructure(man_linear,sol_linear,recordplot;filename="linearactuating.mp4")

tstops = [0,20,40,60,80,100]

function pyplot2(man_linear, sol_linear, tstops)

    steps = [findfirst((x)->x>i,sol_linear.ts) for i = tstops]

    fig,axs_raw = plt.subplots(2,3,figsize=(9,6))
    axs = permutedims(axs_raw,[2,1])
    for (i,step) in enumerate(steps)
        TR.distribute_q_to_rbs!(man_linear,sol_linear.qs[step])
        bars_segs,strings_segs = bars_and_strings_segs(man_linear)
        ax = axs[i]
        ax.add_collection(bars_segs)
        ax.add_collection(strings_segs)
        ax.set_ylim(-150,20)
        ax.set_xlim(-20,150)
        ax.set_aspect("equal")
        ax.set_title("t=$(tstops[i])s")
        ax.grid(true)
        ax.set_yticks(collect(0:-50:-150))
        if i <= 3
            ax.set_xticklabels([])
            ax.xaxis.label.set_visible(false)
        end
        if !(i ∈[1,4])
            ax.set_yticklabels([])
            ax.yaxis.label.set_visible(false)
        end
    end
    fig
end

fig = pyplot2(man_linear, sol_linear, tstops)
fig.savefig("linearactuating_c=$c.png",dpi=300,bbox_inches="tight")
