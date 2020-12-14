using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using Parameters
using Makie
AbstractPlotting.__init__()
plot(rand(10))
import PyPlot; const plt = PyPlot
plt.pygui(true)
using LaTeXStrings
#using DifferentialEquations
#using ForwardDiff
# using DynamicPolynomials
# using HomotopyContinuation
using RecursiveArrayTools
using NLsolve
using Interpolations
using Printf
using CoordinateTransformations
using CSV
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

cd("examples\\links")
include("links_define.jl")
include("link_u_plotting.jl")
include("link_inverse.jl")
include("../analysis.jl")

function simulate_linkn(;dt=0.01,c=0.0,tend=10.0)
    n = 2
    h = 6.6e-2
    R = RotX(0.0)
    linkn = links(n,h,R;c)
    q0,q̇0= TR.get_q(linkn)
    _,_,a = inverse_with_energy(linkn,deepcopy(linkn),build_Y(linkn),0.05)
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

    # reflinkn = links(n,h,RotXY(π/24,π/24))
    reflinkn = links(n,h,RotXY(-π/24,π/24))

    fig = plt.figure(figsize=(5,4))
    ax = fig.add_subplot(1,1,1,projection="3d")
    bars,strings = bars_and_strings_segs_3D(reflinkn)
    ax.add_collection3d(strings)
    ax.add_collection3d(bars)
    set_ax!(ax,0.2)
    fig.tight_layout()

    refq0,refq̇0,refλ0 = TR.get_initial(reflinkn)

    prob = TS.DyProblem(dynfuncs(linkn,refq0),refq0,refq̇0,refλ0,(0.0,tend))
    sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-14)

    linkn,sol,fig
end

linkn_undamped, sol_undamped,fig = simulate_linkn(tend=2000.0)
fig.savefig("undamped_initial_position.png",dpi=300,bbox_inches="tight")
fig.show()

linkn_undamped, sol_undamped,fig = simulate_linkn(dt=0.001,tend=100.0)
linkscene = plotstructure(linkn_undamped,sol_undamped,sliderplot)
test_slackness(linkn_undamped,sol_undamped)


function plot_trajectory(tgstruct,sol,step_range=:)
    rp1 = VectorOfArray(Vector{Vector{Float64}}())
    for q in sol.qs
        TR.distribute_q_to_rbs!(tgstruct,q)
        push!(rp1,tgstruct.rigidbodies[2].state.p[1])
    end

    fig,_ = plt.subplots(1,1,figsize=(8,6))
    grid = (2,4)
    axs = [plt.subplot2grid(grid,(0,2i-2),1,2,fig=fig) for i = 1:2]
    push!(axs,plt.subplot2grid(grid,(1,1),1,2,fig=fig))
    trange = sol.ts[step_range]
    tmin = trange[1]
    tmax = trange[end]

    axs[1].plot(trange,rp1[1,step_range])
    axs[1].set_ylabel("Coordinate (m)")
    axs[1].set_xlabel("Time (s)")
    axs[1].set_xlim(tmin,tmax)
    axs[1].set_title("(a) x",y=-0.5)

    axs[2].plot(trange,rp1[2,step_range])
    axs[2].set_ylabel("Coordinate (m)")
    axs[2].set_xlabel("Time (s)")
    axs[2].set_xlim(tmin,tmax)
    axs[2].set_title("(b) y",y=-0.5)

    axs[3].plot(trange,rp1[3,step_range])
    axs[3].set_ylabel("Coordinate (m)")
    axs[3].set_xlabel("Time (s)")
    axs[3].set_xlim(tmin,tmax)
    axs[3].set_title("(c) z",y=-0.5)

    fig.tight_layout()
    fig
end

fig = plot_trajectory(linkn_undamped,sol_undamped,1:200)
fig.savefig("undamped_trajectory.png",dpi=300,bbox_inches="tight")

function pyplot_energy(tgstruct,sol)
    kes,epes,_,es,es_err = analyse_energy(tgstruct,sol,elasticity=true)

    fig,ax = plt.subplots(1,1,figsize=(4,3))
    ax.plot(sol.ts,es_err)
    ax.set_xlim(sol.ts[1],sol.ts[end])
    ax.set_ylim(-0.01,0.01)
    ax.set_ylabel("Energy Rel. Err.")
    ax.set_xlabel("Time (s)")

    fig.tight_layout()
    fig
end

fig_undamped_energy = pyplot_energy(linkn_undamped,sol_undamped)
fig_undamped_energy.savefig("undamped_energy.png",dpi=300,bbox_inches="tight")



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
