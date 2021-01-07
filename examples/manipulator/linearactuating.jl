using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using Makie
AbstractPlotting.__init__()
plot(rand(10))
import PyPlot; const plt = PyPlot
plt.pygui(true)
using LaTeXStrings
using RecursiveArrayTools
using Interpolations
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

cd("examples/manipulator")
include("man_define.jl")
include("man_plotting.jl")
include("../analysis.jl")


function actuation_interpolation(target_t,initial_a, target_a,
                                Alg=Quadratic(Flat(OnGrid())))
    itps = [interpolate([a0,af],
            BSpline(Alg)) for (a0,af) in zip(initial_a,target_a)]
    function itp(t)
        [itp(t) for itp in itps]
    end
    function inner_aitp(t)
        if t < target_t
            itp(1+t/target_t)
        else
            target_a
        end
    end
end

function simulate_linearactuating(ndof=6;k,c,target_t=10.0,tend=20.0,Alg=Linear())
    manipulator = man_ndof(ndof,k=k,c=c)
    Y = build_Y(manipulator)
    q0,q̇0,λ0 = TR.get_initial(manipulator)

    refman = man_ndof(ndof,θ=-π/12,k=k) # reference
    _,_,initial_a = TR.inverse(manipulator,deepcopy(manipulator),Y)
    _,_,target_a = TR.inverse(manipulator,deepcopy(refman),Y)

    aitp = actuation_interpolation(target_t,initial_a,target_a,Alg)

    function linearactuate(tgstruct,q0,aitp)

        M = TR.build_massmatrix(tgstruct)
        Φ = TR.build_Φ(tgstruct,q0)
        A = TR.build_A(tgstruct)

        function F!(F,q,q̇,t)
            TR.actuate!(tgstruct,aitp(t))
            TR.reset_forces!(tgstruct)
            TR.distribute_q_to_rbs!(tgstruct,q,q̇)
            TR.update_strings_apply_forces!(tgstruct)
            TR.assemble_forces!(F,tgstruct)
        end

        M,Φ,A,F!,nothing

        #A,Φ,∂T∂q̇!,F!,M!,nothing
    end



    dt = 0.01 # Same dt used for PID AND Dynamics solver
    prob = TS.DyProblem(linearactuate(manipulator,q0,aitp),q0,q̇0,λ0,(0.0,tend))

    sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-14)
    manipulator, sol, aitp
end

# bars_and_strings_segs(man_linear)
k = 1.e2
c = 0
target_t = 20
tend = 40

# Alg =
# Alg = Linear()

man_linear, sol_linear, aitp_linear = simulate_linearactuating(;k,c,target_t,tend,Alg=Linear())
man_quadra, sol_quadra, aitp_quadra = simulate_linearactuating(;k,c,target_t,tend,Alg=Quadratic(Flat(OnGrid())))
# scene,_ = plotstructure(man_linear)
function plot_restlen(man,sol,aitp_linear,aitp_quadra)
    ts = sol.ts
    u0 = TR.get_original_restlen(man)
    Y = build_Y(man)
    fig, ax = plt.subplots(1,1,figsize=(4,3))
    rl_linear = VectorOfArray([Y*a + u0 for a in aitp_linear.(ts)])
    rl_quadra = VectorOfArray([Y*a + u0 for a in aitp_quadra.(ts)])
    ax.plot(ts,rl_linear[1,:], label="Linear")
    ax.plot(ts,rl_quadra[1,:], label="Quadratic")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Rest Length (m)")
    ax.grid("on")
    ax.legend()
    fig.tight_layout()
    fig
end
fig = plot_restlen(man_linear,sol_linear,aitp_linear,aitp_quadra)
fig.savefig("restlen_k=$(k)_c=$(c)_t=$(target_t).png",dpi=300,bbox_inches="tight")

plotstructure(man_linear,sol_linear,sliderplot)
# plotstructure(man_linear,sol_linear,recordplot;filename="linearactuating.mp4")

tstops = [0,10,20]
# tstops = [0,24,48,72,96,120]

function pyplot2(man_linear, sol_linear, tstops)

    steps = [findfirst((x)->x>i,sol_linear.ts) for i = tstops]

    fig,axs = plt.subplots(1,3,figsize=(12,5))
    alphabet = join('a':'z')
    for (i,step) in enumerate(steps)
        TR.distribute_q_to_rbs!(man_linear,sol_linear.qs[step])
        bars_segs,strings_segs = bars_and_strings_segs(man_linear)
        ax = axs[i]
        ax.add_collection(bars_segs)
        ax.add_collection(strings_segs)
        ax.set_ylim(-1.00,0.20)
        ax.set_xlim(-0.20,1.50)
        ax.set_aspect("equal")
        ax.set_title("($(alphabet[i])) "*
                    latexstring("t=$(tstops[i])\\mathrm{s}"),
                        y = -0.3)
        ax.grid(true)
        # ax.set_xticks([0,0.5,1])
        # ax.set_yticks([-1.5,-1.0,-0.5,0])
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
fig.savefig("kinematics_control_1x3_k=$(k)_c=$(c)_t=$(target_t).png",dpi=300,bbox_inches="tight")


function plot_energy(man_linear,sol_linear,aitp_linear,man_quadra, sol_quadra,aitp_quadra)
    ts = sol_linear.ts
    as_linear = aitp_linear.(ts)
    as_quadra = aitp_quadra.(ts)
    kes_linear,epes_linear,_ = analyse_energy(man_linear,sol_linear,as_linear;gravity=false,elasticity=true,factor=1)
    kes_quadra,epes_quadra,_ = analyse_energy(man_quadra,sol_quadra,as_quadra;gravity=false,elasticity=true,factor=1)
    fig,axs = plt.subplots(1,2,figsize=(8,4))
    xmin = ts[1]
    xmax = ts[end]
    axs[1].plot(ts, kes_linear, label="Linear")
    axs[1].plot(ts, kes_quadra, label="Quadratic")
    axs[1].set_xlim(xmin,xmax)
    axs[1].set_ylim(-0.001,0.010)
    axs[1].set_xlabel("Time (s)")
    axs[1].set_ylabel("Energy (J)")
    axs[1].set_title("(a) Kinetic Energy", y= -0.3)
    axs[1].grid(true)
    axs[1].legend()
    axs[2].plot(ts, epes_linear, label="linear")
    axs[2].plot(ts, epes_quadra, label="Quadratic")
    axs[2].set_xlim(xmin,xmax)
    axs[2].set_ylim(0.5,0.6)
    axs[2].set_xlabel("Time (s)")
    axs[2].set_ylabel("Energy (J)")
    axs[2].set_title("(b) Elastic Potential Energy", y = -0.3)
    axs[2].grid(true)
    axs[2].legend()
    fig.tight_layout()
    fig
end
fig = plot_energy(man_linear,sol_linear,aitp_linear,man_quadra, sol_quadra,aitp_quadra)
fig.savefig("energy_k=$(k)_c=$(c)_t=$(target_t).png",dpi=300,bbox_inches="tight")

T = 5
ω = 1/T
π/12
