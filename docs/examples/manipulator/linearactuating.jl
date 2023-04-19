using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using Makie
import PyPlot as plt
plt.pygui(true)
using LaTeXStrings
using RecursiveArrayTools
using Interpolations
using Revise
import TensegrityRobots as TR

cd("examples/manipulator")
includet("man_define.jl")
includet("man_plotting.jl")
includet("man_compare.jl")
includet("../analysis.jl")
include("../plot_helpers.jl")
set_pyplot()
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

    refman = man_ndof(ndof,θ=-π/5,k=k) # reference
    _,_,initial_a = TR.inverse(manipulator,deepcopy(manipulator),Y)
    _,_,target_a = TR.inverse(manipulator,deepcopy(refman),Y)

    aitp = actuation_interpolation(target_t,initial_a,target_a,Alg)

    function linearactuate(tr,aitp)
        @unpack tg = tr
        M = TR.build_massmatrix(tg)
        Φ = TR.build_Φ(tg)
        A = TR.build_A(tg)
        # Q̃ = TR.build_Q̃(tg)
        function F!(F,q,q̇,t)
            TR.actuate!(tr,aitp(t))
            TR.reset_forces!(tg)
            TR.distribute_q_to_rbs!(tg,q,q̇)
            TR.update_cables_apply_forces!(tg)
            TR.assemble_forces!(F,tg)
        end

        M,Φ,A,F!,nothing
    end

    dt = 0.01 # Same dt used for PID AND Dynamics solver
    prob = TR.DyProblem(linearactuate(manipulator,aitp),manipulator,(0.0,tend))

    TR.solve!(manipulator,prob,TR.Zhong06(),dt=dt,ftol=1e-14)
    actuation_trajs = aitp.(manipulator.traj.ts)
    for actuation_traj in actuation_trajs
        foreach((x,a)->push!(x.traj.us,a),manipulator.hub.actuators,actuation_traj)
    end
    manipulator,aitp
end

# bars_and_cables_segs(man_linear)
k = 1.e2
c = 50.0
target_t = 20.0
tend = 40.0
man_linear, aitp_linear = simulate_linearactuating(;k,c,target_t,tend,Alg=Linear())
man_quadra, aitp_quadra = simulate_linearactuating(;k,c,target_t,tend,Alg=Quadratic(Flat(OnGrid())))
# scene,_ = plotstructure(man_linear)
function plot_restlen(man,aitp_linear,aitp_quadra)
    ts = man.traj.ts
    u0 = TR.get_original_restlen(man)
    Y = build_Y(man)
    fig, ax = plt.subplots(1,1,figsize=(Single_width,2.5))
    rl_linear = VectorOfArray([Y*a + u0 for a in aitp_linear.(ts)])
    rl_quadra = VectorOfArray([Y*a + u0 for a in aitp_quadra.(ts)])
    ax.plot(ts,rl_linear[1,:], ls=linestyles.xs[1], label="Linear")
    # ax.plot(ts,rl_quadra[1,:], ls=linestyles.xs[2], label="Quadratic")
    ax.set_ylim(0.28,0.34)
    ax.set_xlim(0,25)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Rest Length (m)")
    ax.grid("on")
    # ax.legend()
    fig.tight_layout()
    fig
end
fig = plot_restlen(man_linear,aitp_linear,aitp_quadra)
fig.savefig("restlen_k=$(k)_c=$(c)_t=$(target_t).png",dpi=300,bbox_inches="tight")

plotstructure(man_linear)
# plotstructure(man_linear,sol_linear,recordplot;filename="linearactuating.mp4")

tstops = [0,5,10,15,20,25]

function pyplot2(man_linear,tstops)
    @unpack tg, traj = man_linear
    steps = [findfirst((x)->x>i,traj.ts) for i = tstops]

    fig,axs_raw = plt.subplots(2,3,figsize=(Full_width,4.5))
    axs = permutedims(axs_raw,[2,1])
    for (i,step) in enumerate(steps)
        TR.distribute_q_to_rbs!(tg,traj.qs[step])
        bars_segs,cables_segs = bars_and_cables_segs(tg)
        ax = axs[i]
        ax.add_collection(bars_segs)
        ax.add_collection(cables_segs)
        ax.set_ylim(-1.00,0.20)
        ax.set_xlim(-0.25,1.35)
        ax.set_aspect("equal")
        ax.set_title("($(alphabet[i])) "*
                    latexstring("t=$(tstops[i])\\mathrm{s}"),
                        y = -0.35)
        ax.grid(true)
        # ax.set_xticks([0,0.5,1])
        # ax.set_yticks([-1.5,-1.0,-0.5,0])
        if i ∈[0]
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
fig = pyplot2(man_linear,tstops)
fig.savefig("kinematics_control_2x3_k=$(k)_c=$(c)_t=$(target_t).png",dpi=300,bbox_inches="tight")

function plot_trajectory(man_linear,man_quadra)
    rp8_linear = get_trajectory(man_linear,7,2)
    rp8_quadra = get_trajectory(man_quadra,7,2)
    ts = man_linear.traj.ts
    fig,axs = plt.subplots(1,2,figsize=(Full_width,3))
    xmin = ts[1]
    xmax = ts[end]
    axs[1].plot(ts, rp8_linear[1,:], ls=linestyles.xs[1], label="Linear")
    axs[1].plot(ts, rp8_quadra[1,:], ls=linestyles.xs[2], label="Quadratic")
    axs[1].set_xlim(xmin,xmax)
    # axs[1].set_ylim(-0.001,0.010)
    axs[1].set_xticks([0,10,20,30,40])
    axs[1].set_xlabel("Time (s)")
    axs[1].set_ylabel("Coordinate (m)")
    axs[1].set_title("(a) x₈", y= -0.5)
    axs[1].grid(true)
    # axs[1].legend()
    axs[2].plot(ts, rp8_linear[2,:], ls=linestyles.xs[1], label="linear")
    axs[2].plot(ts, rp8_quadra[2,:], ls=linestyles.xs[2], label="Quadratic")
    axs[2].set_xlim(xmin,xmax)
    # axs[2].set_ylim(0.5,0.6)
    axs[2].set_xticks([0,10,20,30,40])
    axs[2].set_xlabel("Time (s)")
    axs[2].set_ylabel("Coordinate (m)")
    axs[2].set_title("(b) y₈", y = -0.5)
    axs[2].grid(true)
    axs[2].legend(loc="upper left", bbox_to_anchor=(1.1, 1.0))
    fig.tight_layout()
    fig
end
fig = plot_trajectory(man_linear,man_quadra)
fig.savefig("trajectory_k=$(k)_c=$(c)_t=$(target_t).png",dpi=300,bbox_inches="tight")

function plot_energy(man_linear,man_quadra)
    ts = man_linear.traj.ts
    kes_linear,epes_linear,_ = analyse_energy(man_linear;actuation=true, gravity=false,elasticity=true,factor=1)
    kes_quadra,epes_quadra,_ = analyse_energy(man_quadra;actuation=true, gravity=false,elasticity=true,factor=1)
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
fig = plot_energy(man_linear,man_quadra)
fig.savefig("energy_k=$(k)_c=$(c)_t=$(target_t).png",dpi=300,bbox_inches="tight")
