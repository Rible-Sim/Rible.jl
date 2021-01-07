using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using Parameters
using Makie
AbstractPlotting.__init__()
plot(rand(10))
import PyPlot;
const plt = PyPlot;
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
using TensegritySolvers;
const TS = TensegritySolvers;
using TensegrityRobot
const TR = TensegrityRobot

cd("examples\\links")
include("links_define.jl")
include("link_u_plotting.jl")
include("link_inverse.jl")
include("../analysis.jl")

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
    q0, q̇0, λ0 = TR.get_initial(linkn)
    function linearactuate(tg, q0, aitp)
        M = TR.build_massmatrix(tg)
        Φ = TR.build_Φ(tg, q0)
        A = TR.build_A(tg)
        # Q̃ = TR.build_Q̃(tg)
        function F!(F, q, q̇, t)
            TR.actuate!(tg, aitp(t))
            TR.reset_forces!(tg)
            TR.distribute_q_to_rbs!(tg, q, q̇)
            TR.update_strings_apply_forces!(tg)
            TR.assemble_forces!(F, tg)
        end
        M, Φ, A, F!, nothing
    end

    dt = 0.01
    aitp = actuation_interpolation(target_t, a0, af, Alg)
    prob = TS.DyProblem(linearactuate(linkn, q0, aitp), q0, q̇0, λ0, (0.0, tend))
    sol = TS.solve(prob, TS.Zhong06(), dt = dt, ftol = 1e-14)
    linkn, sol, aitp
end

epe0, epef = 0.05, 0.05
epe0, epef = 0.10, 0.05
epe0, epef = 0.05, 0.10
c = 0.0
c = 1.0

linkn_linear, sol_linear, aitp_linear = kinematics_control(; epe0, epef, c, Alg = Linear())
linkn_quadra, sol_quadra, aitp_quadra = kinematics_control(; epe0, epef, c, Alg = Quadratic(Flat(OnGrid())))

plotstructure(linkn_linear, sol_linear, sliderplot)

function pyplot_energy(
    linkn_linear,
    sol_linear,
    aitp_linear,
    linkn_quadra,
    sol_quadra,
    aitp_quadra,
    epe0,
    epef,
    )
    ts = sol_linear.ts
    as_linear = aitp_linear.(ts)
    as_quadra = aitp_quadra.(ts)
    kes_linear, epes_linear, _ =
        analyse_energy(linkn_linear, sol_linear, as_linear; elasticity = true)
    kes_quadra, epes_quadra, _ =
        analyse_energy(linkn_quadra, sol_quadra, as_quadra; elasticity = true)
    fig, ax = plt.subplots(1, 2, figsize = (8, 4))

    ax[1].plot(ts, kes_linear, label = "Linear")
    ax[1].plot(ts, kes_quadra, label = "Quadratic")
    ax[1].set_xlim(ts[1], ts[end])
    ax[1].set_ylim([-0.2, 5] .* 1e-5)
    ax[1].set_xlabel("Time (s)")
    ax[1].set_ylabel("Energy (J)")
    ax[1].set_title("(a) Kinetic Energy", y = -0.3)
    ax[1].grid("on")
    ax[1].legend()
    ax[2].plot(ts, epes_linear, label = "Linear")
    ax[2].plot(ts, epes_quadra, label = "Quadratic")
    ax[2].set_xlim(ts[1], ts[end])
    ax[2].set_ylim(0.125, 0.325)
    ax[2].set_xlabel("Time (s)")
    ax[2].set_ylabel("Energy (J)")
    ax[2].set_title("(b) Elastic Potential Energy", y = -0.3)
    ax[2].grid("on")
    ax[2].legend()
    fig.tight_layout()
    fig
end

fig_energy = pyplot_energy(
    linkn_linear,
    sol_linear,
    aitp_linear,
    linkn_quadra,
    sol_quadra,
    aitp_quadra,
    epe0,
    epef,)
fig_energy.savefig("kinematics_control_energy_epe0=$(epe0)_epef=$(epef)_c=$c.png",
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
u0 = TR.get_original_restlen(linkn_linear)
fig = plot_actuation(sol_linear.ts, 10.0, a0, af, u0)
fig.savefig("linearactuating_actuation.png", dpi = 300, bbox_inches = "tight")

tstops = [0, 5, 10]
function pyplot_linearactuating(tgstruct, sol, tstops)
    steps = [findfirst((x) -> x > i, sol.ts) for i in tstops]
    nsps = length(tstops)

    fig = plt.figure(figsize = (12, 5))
    axs_raw = [fig.add_subplot(1, 3, i, projection = "3d") for i = 1:nsps]
    axs = axs_raw
    alphabet = join('a':'z')
    for (i, ax) in enumerate(axs)
        TR.distribute_q_to_rbs!(tgstruct, sol.qs[steps[i]])
        bars, strings = bars_and_strings_segs_3D(tgstruct)
        ax.add_collection3d(bars)
        ax.add_collection3d(strings)
        ax.set_title("($(alphabet[i])) t=$(tstops[i]) (s)", y = -0.12)
        set_ax!(ax, 0.32; elev = 11, azim = -76)
    end
    fig.tight_layout()
    fig
end
fig = pyplot_linearactuating(linkn_linear, sol_linear, tstops)
fig.savefig("linearactuating_c=$(c)_t=$tend.png", dpi = 300, bbox_inches = "tight")
plt.close(fig)
