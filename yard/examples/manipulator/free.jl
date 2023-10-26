using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using RecursiveArrayTools
using BenchmarkTools
using Makie
import PyPlot as plt
plt.pygui(true)
using LaTeXStrings
using CSV
using Revise
import Rible as RB
cd("examples/manipulator")
includet("man_define.jl")
includet("man_plotting.jl")
includet("man_compare.jl")
includet("../analysis.jl")
includet("../plot_helpers.jl")
set_pyplot()
k = 2.e3; c = 0.0
# dt = 0.01 # Same dt used for PID AND Dynamics solver

manipulator = man_ndof(2;k,c)

fig,ax = plt.subplots(1,1,figsize=(Single_width,2.5))
bars,cables = bars_and_cables_segs(manipulator.st)
ax.add_collection(bars)
ax.add_collection(cables)
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_ylim(-0.3,0.2)
ax.set_xlim(-0.05,0.65)
ax.set_aspect("equal")
ax.grid(true)
fig.tight_layout()
fig.savefig("undamped_initial_position.png",dpi=300,bbox_inches="tight")
function simulate_free(ndof;k,c,unit="mks",dt=0.01,tend=10.0,verbose=false)
    manipulator = man_ndof(ndof;k,c,unit)

    function free(tr)
        st = tr.st
        M = RB.build_massmatrix(st)
        Φ = RB.build_Φ(st)
        A = RB.build_A(st)
        Q̃ = RB.build_Q̃(st)
        function F!(F,q,q̇,t)
            RB.reset_forces!(st)
            RB.distribute_q_to_rbs!(st,q,q̇)
            RB.update_cables_apply_forces!(st)
            RB.apply_gravity!(st)
            RB.assemble_forces!(F,st)
            # @show isapprox(F,Q̃*RB.fvector(st))
        end
        Jac_Γ = RB.build_Jac_Γ(st)
        function Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)
            ∂Γ∂q,∂Γ∂q̇ = Jac_Γ(q,q̇)
            Q̃*∂Γ∂q,Q̃*∂Γ∂q̇
        end
        M,Φ,A,F!,Jac_F!
    end


    prob = RB.DyProblem(free(manipulator),manipulator,(0.0,tend))

    RB.solve!(manipulator,prob,RB.Zhong06();dt,ftol=1e-14,verbose)
    manipulator
end


function plot_trajectory(tr,res)
    rp4 = get_trajectory(tr,3,2)
    ts = tr.traj.ts
    tmin = ts[begin]
    tmax = ts[end]

    t = res("time")[begin:50:end]
    p4x = res("P4X")[begin:50:end]
    p4y = res("P4Y")[begin:50:end]

    fig,axs = plt.subplots(1,2,figsize=(Full_width,2.6))

    axs[1].plot(ts,rp4[1,:],lw=1)
    axs[1].plot(t,p4x,lw=1,linestyle="",marker="x")
    axs[1].set_ylabel("Coordinate (m)")
    axs[1].set_xlabel("Time (s)")
    axs[1].set_xlim(tmin,tmax)
    axs[1].text(-0.35, 1, "(a)", transform=axs[1].transAxes,size = 14)
    axs[1].grid(true)

    axs[2].plot(ts,rp4[2,:],lw=1,label="Proposed")
    axs[2].plot(t,p4y,linestyle="",marker="x",label="Adams")
    axs[2].set_ylabel("Coordinate (m)")
    axs[2].set_xlabel("Time (s)")
    axs[2].set_xlim(tmin,tmax)
    axs[2].text(-0.35, 1, "(b)", transform=axs[2].transAxes,size = 14)
    axs[2].grid(true)
    axs[2].legend(loc="upper left", bbox_to_anchor=(1.05, 1.04))

    fig.tight_layout()
    fig
end

man_free = simulate_free(2;k,c=0,dt=0.001,tend=4.0)
fig = plot_trajectory(man_free,RB.AdamsResults("1AL.xml"))
fig.savefig("undamped_trajectory.png",dpi=300,bbox_inches="tight")

man_free = simulate_free(2;k,c=5,dt=0.001,tend=4.0)
fig = plot_trajectory(man_free,RB.AdamsResults("1ALC.xml"))
fig.savefig("damped_trajectory.png",dpi=300,bbox_inches="tight")


man_free_long = simulate_free(2;k,c=0,dt=0.001,tend=400.0)
@benchmark simulate_free(2;k,c=0,dt=0.001,tend=400.0)
man2 = man_ndof(2;k,c)

plotstructure(man_free)

l0s = RB.get_cables_len(man2)

function pyplot_energy(tr,l0s)
    @unpack st, traj = tr
    # res = RB.AdamsResults("1BWI3E11H3.xml")
    # t = res("time")
    ts = traj.ts
    u0 = RB.get_original_restlen(tr)
    # ke,epe,gpe,energy,energy_error = adams_to_energy(res,tr,l0s,u0)
    kes,epes,gpes,es,es_err = analyse_energy(tr;gravity=true,elasticity=true)

    # @show maximum(energy_error)/maximum(es_err)

    fig,ax = plt.subplots(1,1,figsize=(Onehalf_width,4))
    ax.plot(ts,es_err,label="Zu")
    # ax.plot(t,energy_error,marker="x",label="Adams")
    ax.set_xlim(ts[begin],ts[end])
    ax.set_ylim(-0.0002,0.0013)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Energy Rel. Err.")
    ax.legend(loc="upper right")
    ax.grid(true)
    fig.tight_layout()
    fig
end

fig = pyplot_energy(man_free_long,l0s)

fig.savefig("undamped_energy_error.png",dpi=300,bbox_inches="tight")

function plot_constraint_drift(tr)
    @unpack traj = tr
    @unpack ts = traj
    tmin = ts[begin]
    tmax = ts[end]
    pe, ve = compute_contraint_error(tr)
    # @show minimum(pe)
    fig, axs = plt.subplots(1,2,figsize=(2Single_width,3))
    axs[1].plot(ts,pe,lw=1)
    axs[1].set_ylabel("Error")
    axs[1].set_xlabel("Time (s)")
    # axs[1].set_yscale("log")
    axs[1].set_ylim(-1e-16,3e-16)
    axs[1].set_xlim(tmin,tmax)
    # axs[1].set_title("(a)",y=-0.4)
    axs[1].text(-0.18, 1.0, "(a)", transform=axs[1].transAxes,size = 14)
    axs[1].grid(true)

    axs[2].plot(ts,ve,lw=1)
    axs[2].set_ylabel("Error")
    axs[2].set_xlabel("Time (s)")
    axs[2].set_ylim(-4e-4,2e-3)
    axs[2].set_xlim(tmin,tmax)
    axs[2].text(-0.3, 1.0, "(b)", transform=axs[2].transAxes,size = 14)
    axs[2].grid(true)
    fig.tight_layout()
    fig
end

fig = plot_constraint_drift(man_free)
fig.savefig("constraint_drift.png",dpi=300,bbox_inches="tight")
