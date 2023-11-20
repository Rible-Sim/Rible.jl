using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using Parameters
using RecursiveArrayTools
using GLMakie
import PyPlot as plt
plt.pygui(true)
using LaTeXStrings
#using DifferentialEquations
#using ForwardDiff
# using DynamicPolynomials
# using Printf
using Interpolations
using NLsolve
using CoordinateTransformations
# using HomotopyContinuation
using Revise
import Rible as RB

cd("examples/links")
includet("links_define.jl")
includet("link_u_plotting.jl")
includet("link_inverse.jl")
includet("../analysis.jl")
includet("../plot_helpers.jl")
set_pyplot()
function find_actuation(epe0 = 0.05, epef = 0.10, c = 0.0)
    n = 4
    h = 0.066
    linkn = links(n, h, RotX(0.0); k = 3e1, c)
    Y = build_Y(linkn)
    reflinkn = links(n, h, RotY(π / 12); k = 3e1)
    _, _, a0 = inverse_with_energy_3(linkn, deepcopy(linkn), Y, fill(epe0, 3))
    _, _, af = inverse_with_energy_3(linkn, reflinkn, Y, fill(epef, 3))
    linkn, a0, af
end

function actuation_interpolation(
    target_t,
    initial_a,
    target_a,
    Alg = Quadratic(Flat(OnGrid())),)
    itps = [interpolate([a0, af], BSpline(Alg)) for (a0, af) in zip(initial_a, target_a)]
    function itp(t)
        [itp(t) for itp in itps]
    end
    function inner_aitp(t)
        if t < target_t
            itp(1 + t / target_t)
        else
            target_a
        end
    end
end

function kinematics_control(; epe0, epef, c, target_t = 10.0, tend = 20.01, Alg = Linear())

    linkn, a0, af = find_actuation(epe0, epef, c)
    q0, q̇0, λ0 = RB.get_initial(linkn.st)
    function linearactuate(tr, aitp)
        st = tr.st
        M = RB.build_massmatrix(st)
        Φ = RB.build_Φ(st)
        A = RB.build_A(st)
        Q̃ = RB.build_Q̃(st)
        function F!(F, q, q̇, t)
            RB.actuate!(tr, aitp(t))
            RB.reset_forces!(st)
            RB.distribute_q_to_rbs!(st, q, q̇)
            RB.update_strings_apply_forces!(st)
            RB.assemble_forces!(F, st)
        end
        Jac_Γ = RB.build_Jac_Γ(st)
        function Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)
            ∂Γ∂q,∂Γ∂q̇ = Jac_Γ(q,q̇)
            Q̃*∂Γ∂q,Q̃*∂Γ∂q̇
        end
        M, Φ, A, F!, nothing
    end

    dt = 0.01
    aitp = actuation_interpolation(target_t, a0, af, Alg)
    prob = RB.DynamicsProblem(linkn,(x)->linearactuate(x, aitp), (0.0, tend))
    RB.solve!(prob, RB.Zhong06(), dt = dt, ftol = 1e-14)
    actuation_trajs = aitp.(linkn.traj.ts)
    for actuation_traj in actuation_trajs
        foreach((x,a)->push!(x.traj.us,a),linkn.hub.actuators,actuation_traj)
    end
    linkn
end

epe0, epef = 0.05, 0.05
epe0, epef = 0.10, 0.05
epe0, epef = 0.05, 0.10
c = 0.0
c = 1.0

linkn_linear = kinematics_control(; epe0, epef, c, Alg = Linear())
linkn_quadra = kinematics_control(; epe0, epef, c, Alg = Quadratic(Flat(OnGrid())))

plotstructure(linkn_linear)
plotstructure(linkn_linear;do_record=true)
function pyplot_energy(linkn_linear,linkn_quadra,epe0,epef)
    ts = linkn_linear.traj.ts
    kes_linear, epes_linear, _ =
        analyse_energy(linkn_linear; elasticity = true, actuation = true)
    kes_quadra, epes_quadra, _ =
        analyse_energy(linkn_quadra; elasticity = true, actuation = true)
    fig, ax = plt.subplots(1, 2, figsize = (Full_width, 3))

    ax[1].plot(ts, kes_linear, label = "Linear")
    ax[1].plot(ts, kes_quadra, label = "Quadratic")
    ax[1].set_xlim(ts[1], ts[end])
    ax[1].set_ylim([-0.2, 5] .* 1e-5)
    ax[1].set_xlabel("Time (s)")
    ax[1].set_ylabel("Energy (J)")
    ax[1].set_title("(a)", y = -0.4)
    ax[1].grid("on")
    # ax[1].legend()
    ax[2].plot(ts, epes_linear, label = "Linear")
    ax[2].plot(ts, epes_quadra, label = "Quadratic")
    ax[2].set_xlim(ts[1], ts[end])
    ax[2].set_ylim(0.125, 0.325)
    ax[2].set_xlabel("Time (s)")
    ax[2].set_ylabel("Energy (J)")
    ax[2].set_title("(b)", y = -0.4)
    ax[2].grid("on")
    ax[2].legend(loc="upper left",bbox_to_anchor=(1.0, 1.03))
    fig.tight_layout()
    fig
end

fig_energy = pyplot_energy(linkn_linear,linkn_quadra,epe0,epef)
fig_energy.savefig("kinematics_control_energy_epe0=$(epe0)_epef=$(epef)_c=$c.png",
    dpi = 300,bbox_inches = "tight",)

function plot_trajectory_energy(linkn_linear,linkn_quadra,epe0,epef,)
    ts = linkn_linear.traj.ts
    rp7_linear = get_trajectory(linkn_linear,2,2)
    rp7_quadra = get_trajectory(linkn_quadra,2,2)

    kes_linear, epes_linear, _ =
        analyse_energy(linkn_linear; actuation = true, elasticity = true)
    kes_quadra, epes_quadra, _ =
        analyse_energy(linkn_quadra; actuation = true, elasticity = true)
    fig, ax = plt.subplots(2, 2, figsize = (Full_width, 5))

    ax[1].plot(ts, rp7_linear[1,:], ls = linestyles.xs[1], label = "Linear")
    ax[1].plot(ts, rp7_quadra[1,:], ls = linestyles.xs[2], label = "Quadratic")
    ax[1].set_xlim(ts[1], ts[end])
    ax[1].set_ylim(0.078, 0.105)
    ax[1].set_xlabel("Time (s)")
    ax[1].set_ylabel("Coordinate (M)")
    ax[1].set_title("(a)", y = -0.48)
    ax[1].grid("on")
    # ax[1].legend()
    nstep = length(ts)
    halfrange = Int(floor(nstep/2)):nstep
    ax[3].plot(ts[halfrange], rp7_linear[1,halfrange], ls = linestyles.xs[1], label = "Linear")
    ax[3].plot(ts[halfrange], rp7_quadra[1,halfrange], ls = linestyles.xs[2], label = "Quadratic")
    ax[3].set_xlim(10, ts[end])
    ax[3].set_ylim(0.102, 0.104)
    ax[3].set_xlabel("Time (s)")
    ax[3].set_ylabel("Coordinate (M)")
    ax[3].set_title("(b)", y = -0.48)
    ax[3].grid("on")
    ax[3].legend(loc="upper left",bbox_to_anchor=(1.0, 1.05))


    ax[2].plot(ts, kes_linear, ls = linestyles.xs[1], label = "Linear")
    ax[2].plot(ts, kes_quadra, ls = linestyles.xs[2], label = "Quadratic")
    ax[2].set_xlim(ts[1], ts[end])
    ax[2].set_ylim([-0.2, 5] .* 1e-5)
    ax[2].set_xlabel("Time (s)")
    ax[2].set_ylabel("Energy (J)")
    ax[2].set_title("(c)", y = -0.48)
    ax[2].grid("on")
    # ax[2].legend()
    ax[4].plot(ts, epes_linear, ls = linestyles.xs[1], label = "Linear")
    ax[4].plot(ts, epes_quadra, ls = linestyles.xs[2], label = "Quadratic")
    ax[4].set_xlim(ts[1], ts[end])
    ax[4].set_ylim(0.125, 0.325)
    ax[4].set_xlabel("Time (s)")
    ax[4].set_ylabel("Energy (J)")
    ax[4].set_title("(d)", y = -0.48)
    ax[4].grid("on")
    # ax[4].legend()

    fig.tight_layout()
    plt.subplots_adjust(hspace=0.5)
    fig
end

fig_energy = plot_trajectory_energy(linkn_linear,linkn_quadra,epe0,epef)
fig_energy.savefig("kinematics_control_trajectory_energy_epe0=$(epe0)_epef=$(epef)_c=$c.png",
    dpi = 300,bbox_inches = "tight",)


function plot_actuation(ts, target_t, a0, af, u0)
    aitp_linear = actuation_interpolation(target_t, a0, af, Linear())
    aitp_quadra = actuation_interpolation(target_t, a0, af, Quadratic(Flat(OnGrid())))
    firstact(ts, aitp) = [a[1] for a in aitp.(ts)]
    fig, ax = plt.subplots(1, 1, figsize = (4, 3))
    ax.plot(ts, firstact(ts, aitp_linear) .+ u0[1], label = "Linear")
    ax.plot(ts, firstact(ts, aitp_quadra) .+ u0[1], label = "Qauadratic")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Rest length (m)")
    ax.grid("on")
    ax.legend()
    fig.tight_layout()
    fig
end
u0 = RB.get_original_restlen(linkn_linear)
fig = plot_actuation(sol_linear.ts, 10.0, a0, af, u0)
fig.savefig("linearactuating_actuation.png", dpi = 300, bbox_inches = "tight")

tstops = [0, 5, 10]
function pyplot_linearactuating(tgstruct, sol, tstops)
    steps = [findfirst((x) -> x > i, sol.ts) for i in tstops]
    nsps = length(tstops)

    fig = plt.figure(figsize = (Full_width, 3))
    axs_raw = [fig.add_subplot(1, 3, i, projection = "3d") for i = 1:nsps]
    axs = axs_raw

    for (i, ax) in enumerate(axs)
        RB.distribute_q_to_rbs!(tgstruct, sol.qs[steps[i]])
        bars, strings = bars_and_strings_segs_3D(tgstruct)
        ax.add_collection3d(bars)
        ax.add_collection3d(strings)
        ax.set_title("($(alphabet[i])) t=$(tstops[i]) (s)", y = -0.2)
        set_ax!(ax, 0.32; elev = 11, azim = -76)
    end
    fig.tight_layout()
    fig
end
fig = pyplot_linearactuating(linkn_linear, sol_linear, tstops)
fig.savefig("linearactuating_c=$(c).png", dpi=300, bbox_inches = "tight")
plt.close(fig)
