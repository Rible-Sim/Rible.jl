using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using Parameters
import PyPlot; const plt = PyPlot
plt.pygui(true)
#using DifferentialEquations
#using ForwardDiff
# using DynamicPolynomials
# using HomotopyContinuation
using Printf
using CoordinateTransformations
using Makie
AbstractPlotting.__init__()
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

cd("examples\\links")
include("links_define.jl")
include("link_u_plotting.jl")
include("../analysis.jl")

function simulate_linkn(;c=0.0)
    n = 2
    h = 6.6
    R = RotX(0.0)
    linkn = links(n,h,R;c)
    q0,q̇0,λ0 = TR.get_initial(linkn)
    rl,a = inverse2actuation(linkn)
    TR.actuate!(linkn,a)

    function dynfuncs(tg,q0)
        M = TR.build_massmatrix(tg)
        Φ = TR.build_Φ(tg,q0)
        A = TR.build_A(tg)
        # Q̃ = TR.build_Q̃(tg)
        function F!(F,q,q̇,t)
            TR.reset_forces!(tg)
            TR.distribute_q_to_rbs!(tg,q,q̇)
            TR.update_strings_apply_forces!(tg)
            # TR.apply_gravity!(tg)
            # F .= Q̃*TR.fvector(tg)
            TR.assemble_forces!(F,tg)
        end
        M,Φ,A,F!,nothing
    end

    Rx = RotY(π/18)
    reflinkn = links(n,h,Rx)

    refq0,refq̇0,refλ0 = TR.get_initial(reflinkn)

    dt = 0.01
    prob = TS.DyProblem(dynfuncs(linkn,refq0),refq0,refq̇0,refλ0,(0.0,100.0))
    sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-11,verbose=true)

    linkn,sol
end

linkn_undamped, sol_undamped = simulate_linkn()
linkscene = plotstructure(linkn_undamped,sol_undamped,sliderplot)
function pyplot_energy(tgstruct,sol)
    kes,epes,_,es,es_err = analyse_energy(tgstruct,sol,elasticity=true)
    fig, ax = plt.subplots(figsize=(9,6))
    ax.plot(sol.ts,es_err,label="Energy")
    # ax.plot(sol.ts,kes,label="Kinetic Energy")
    # ax.plot(sol.ts,epes,label="Elastic Potential Energy")
    ax.set_xlim(sol.ts[1],sol.ts[end])
    ax.set_ylim(-2e-3,2e-3)
    ax.set_xlabel("Time(s)")
    ax.set_ylabel("Rel. Energy Error")
    fig
end
fig_undamped_energy = pyplot_energy(linkn_undamped,sol_undamped)
fig_undamped_energy.savefig("undamped_energy_error.png",dpi=300,bbox_inches="tight")
plt.close(fig_undamped_energy)

linkn_damped, sol_damped = simulate_linkn(c=5.0)
linkscene = plotstructure(linkn_damped,sol_damped)

function pyplot_energy(tgstruct,sol)
    kes,epes,_,es,es_err = analyse_energy(tgstruct,sol,elasticity=true)
    fig, ax = plt.subplots(figsize=(9,6))
    ax.plot(sol.ts,es,label="Energy")
    # ax.plot(sol.ts,kes,label="Kinetic Energy")
    # ax.plot(sol.ts,epes,label="Elastic Potential Energy")
    ax.set_xlim(sol.ts[1],sol.ts[end])
    # ax.set_ylim(-1e-4,1e-4)
    ax.set_xlabel("Time(s)")
    ax.set_ylabel("Rel. Energy Error")
    fig
end
fig_damped_energy = pyplot_energy(linkn_damped,sol_damped)
fig_damped_energy.savefig("damped_energy.png",dpi=300,bbox_inches="tight")
plt.close(fig_damped_energy)
