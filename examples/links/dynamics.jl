using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using Parameters
using Makie
plot(rand(10))
import PyPlot; const plt = PyPlot
plt.pygui(true)
#using DifferentialEquations
#using ForwardDiff
# using DynamicPolynomials
# using HomotopyContinuation
using Interpolations
using Printf
using CoordinateTransformations
using CSV
AbstractPlotting.__init__()
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

cd("examples\\links")
include("links_define.jl")
include("link_u_plotting.jl")
include("../analysis.jl")

function simulate_linkn(;dt=0.01,c=0.0,tend=10.0)
    n = 2
    h = 6.6e-2
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

    prob = TS.DyProblem(dynfuncs(linkn,refq0),refq0,refq̇0,refλ0,(0.0,tend))
    sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-14,verbose=true)

    linkn,sol
end

linkn_undamped, sol_undamped = simulate_linkn()
linkscene = plotstructure(linkn_undamped,sol_undamped,sliderplot)
function pyplot_energy(tgstruct,sol)
    kes,epes,_,es,es_err = analyse_energy(tgstruct,sol,elasticity=true)
    fig, ax = plt.subplots(3,1,figsize=(6,9))

    ax[1].plot(sol.ts,es)
    ax[1].set_xlim(sol.ts[1],sol.ts[end])
    ax[1].set_ylim(0.05,0.07)
    ax[1].set_xlabel("Time(s)")
    ax[1].set_title("Energy")

    ax[2].plot(sol.ts,kes)
    ax[2].set_xlim(sol.ts[1],sol.ts[end])
    ax[2].set_ylim(-0.005,0.015)
    ax[2].set_yticks([-0.005,0.0,0.005,0.01,0.015])
    ax[2].set_xlabel("Time(s)")
    ax[2].set_title("Elastic Potential Energy")

    ax[3].plot(sol.ts,epes)
    ax[3].set_xlim(sol.ts[1],sol.ts[end])
    ax[3].set_ylim(0.05,0.07)
    ax[3].set_xlabel("Time(s)")
    ax[3].set_title("Kinetic Energy")
    fig.tight_layout()
    fig
end
function pyplot_energy_all(tgstruct,sol)
    kes,epes,_,es,es_err = analyse_energy(tgstruct,sol,elasticity=true)
    fig, ax = plt.subplots(1,1,figsize=(6,3))

    ax.plot(sol.ts,es)
    ax.set_xlim(sol.ts[1],sol.ts[end])
    ax.set_ylim(0.05,0.07)
    ax.set_xlabel("Time(s)")
    ax.set_title("Energy")

    fig.tight_layout()
    fig
end

fig_undamped_energy_all = pyplot_energy_all(linkn_undamped,sol_undamped)
fig_undamped_energy_all.savefig("undamped_energy_all.png",dpi=300,bbox_inches="tight")
fig_undamped_energy = pyplot_energy(linkn_undamped,sol_undamped)
fig_undamped_energy.savefig("undamped_energy.png",dpi=300,bbox_inches="tight")
plt.close(fig_undamped_energy)


c = 0.01
tend = 20.0
linkn_damped_01, sol_damped_01 = simulate_linkn(;dt=0.01,c,tend)
linkscene = plotstructure(linkn_damped_01,sol_damped_01,sliderplot)

fig_damped_energy_01 = pyplot_energy(linkn_damped_01,sol_damped_01)
fig_damped_energy_01.savefig("damped_energy_01.png",dpi=300,bbox_inches="tight")


linkn_damped_001, sol_damped_001 = simulate_linkn(;dt=0.001,c,tend)
linkscene = plotstructure(linkn_damped_001,sol_damped_001,sliderplot)

fig_damped_energy_001 = pyplot_energy(linkn_damped_001,sol_damped_001)
fig_damped_energy_001.savefig("damped_energy_001.png",dpi=300,bbox_inches="tight")
plt.close(fig_damped_energy)


function pyplot_damped_energy(tg1,sol1,tg2,sol2)
    kes1,epes1,_,es1,es_err1 = analyse_energy(tg1,sol1,elasticity=true)
    kes2,epes2,_,es2,es_err2 = analyse_energy(tg2,sol2,elasticity=true)
    fig, ax = plt.subplots(3,1,figsize=(6,9))
    function etp_and_itp(nodes,A)
        extrapolate(interpolate((nodes,),A, Gridded(Linear())), Flat())
    end
    ts1 = sol1.ts
    es2_ts1 = [etp_and_itp(sol2.ts,es2)(t) for t in sol1.ts]
    ax[1].plot(ts1,(es1.-es2_ts1)./abs.(es2_ts1))
    ax[1].set_xlim(ts1[1],ts1[end])
    ax[1].set_xlabel("Time(s)")
    ax[1].set_title("Energy Rel. Err.")

    kes2_ts1 = [etp_and_itp(sol2.ts,kes2)(t) for t in ts1]
    ax[2].set_xlim(ts1[1],ts1[end])
    ax[2].plot(ts1,(kes1.-kes2_ts1))
    ax[2].set_ylim(-0.004,0.004)
    ax[2].set_xlabel("Time(s)")
    ax[2].set_title("Kinetic Energy Abs. Err.")

    epes2_ts1 = [etp_and_itp(sol2.ts,epes2)(t) for t in ts1]
    ax[3].set_xlim(ts1[1],ts1[end])
    ax[3].plot(ts1,(epes1.-epes2_ts1))
    ax[3].set_ylim(-0.004,0.004)
    ax[3].set_xlabel("Time(s)")
    ax[3].set_title("Elastic Potential Energy Abs. Err.")

    # ax.plot(sol.ts,kes,label="Kinetic Energy")
    # ax.plot(sol.ts,epes,label="Elastic Potential Energy")
    # ax.set_xlim(sol.ts[1],sol.ts[end])
    # # ax.set_ylim(-1e-4,1e-4)
    # ax.set_xlabel("Time(s)")
    # ax.set_ylabel("Energy")
    fig.tight_layout()
    fig
end

fig = pyplot_damped_energy(linkn_damped_01,sol_damped_01,linkn_damped_001, sol_damped_001)
fig.savefig("damped_energy_error_c=$c.png",dpi=300,bbox_inches="tight")
