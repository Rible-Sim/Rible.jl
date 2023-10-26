using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using Parameters
using BenchmarkTools
using GLMakie
import PyPlot as plt
plt.pygui(true)
using LaTeXStrings
# using TableView
# using Blink
#using DifferentialEquations
#using ForwardDiff
# using DynamicPolynomials
# using HomotopyContinuation
using Cthulhu
using RecursiveArrayTools
using NLsolve
using Interpolations
using Printf
using CoordinateTransformations
using CSV
using Revise
import Rible as RB

cd("examples/links")
includet("links_define.jl")
includet("link_u_plotting.jl")
includet("link_inverse.jl")
includet("links_compare.jl")
includet("../analysis.jl")
includet("../plot_helpers.jl")
set_pyplot()
function simulate_linkn(;dt=0.01,k=3e1,c=0.0,tend=10.0,verbose=false)
    n = 2
    h = 6.0e-2
    R = RotY(0.0)
    linkn = links(n,h,R;k,c)
    # q0,q̇0= RB.get_q(linkn)
    # _,_,a = inverse_with_energy(linkn,deepcopy(linkn),build_Y(linkn),0.05)
    # RB.actuate!(linkn,a)

    function dynfuncs(bot)
        @unpack st = bot
        M = RB.build_massmatrix(st)
        Φ = RB.build_Φ(st)
        A = RB.build_A(st)
        Q̃ = RB.build_Q̃(st)
        function F!(F,q,q̇,t)
            RB.reset_forces!(st)
            RB.distribute_q_to_rbs!(st,q,q̇)
            RB.update_strings_apply_forces!(st)
            RB.apply_gravity!(st)
            # F .= Q̃*RB.fvector(st)
            RB.assemble_forces!(F,st)
        end
        Jac_Γ = RB.build_Jac_Γ(st)
        function Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)
            ∂Γ∂q,∂Γ∂q̇ = Jac_Γ(q,q̇)
            Q̃*∂Γ∂q,Q̃*∂Γ∂q̇
        end
        M,Φ,A,F!,Jac_F!
    end
    # reflinkn = links(n,h,RotXY(π/24,π/24))
    reflinkn = links(n,h,RotZ(π/18))
    # reflinkn = deepcopy(linkn)

    refq0,refq̇0,refλ0 = RB.get_initial(reflinkn.st)
    RB.set_new_initial!(linkn,refq0,refq̇0)

    prob = RB.SimProblem(linkn,dynfuncs,(0.0,tend))
    RB.solve!(prob,RB.Zhong06();dt,ftol=1e-14,verbose)

    linkn,reflinkn
end

linkn_undamped, reflinkn = simulate_linkn(;k=3e1,c=0.0,tend=10.0,dt=0.001,verbose=true)
# Juno.@profiler simulate_linkn(;k=3e1,c=0.0,tend=100.0,dt=0.001)
plotstructure(linkn_undamped)
@descend_code_warntype simulate_linkn(;k=3e1,c=0.0,tend=4.0,dt=0.001,verbose=true)
displacement_constraint_errors, velocity_constraint_errors = compute_contraint_error(linkn_undamped)
plt.plot(displacement_constraint_errors)
plt.plot(velocity_constraint_errors)
kes,epes,_,es,es_err = analyse_energy(linkn_undamped,gravity=true,elasticity=true)
plt.plot(es_err)
test_slackness(linkn_undamped)

n = 2
h = 6.0e-2
R = RotY(0.0)
linkn = links(n,h,RotZ(π/18))
@code_warntype links(n,h,RotZ(π/18))
@code_warntype simulate_linkn(;k=3e1,c=0.0,tend=4.0,dt=0.001,verbose=true)

fig = plt.figure(figsize=(Single_width,3))
ax = fig.add_subplot(1,1,1,projection="3d")
bars,strings = bars_and_strings_segs_3D(reflinkn.st)
ax.add_collection3d(strings)
ax.add_collection3d(bars)
set_ax!(ax,0.2;elev=24,azim=-46)
fig.tight_layout()
fig.savefig("undamped_initial_position.png",dpi=300,bbox_inches="tight")

u0 = RB.get_original_restlen(linkn_undamped)
plotstructure(linkn_undamped)
plotstructure(linkn_undamped,sol_undamped)
test_slackness(linkn_undamped,sol_undamped)

r1 = get_trajectory(linkn_undamped_long,2,1)
plt.plot(r1[3,begin:4000:end])

res = RB.AdamsResults("2A.xml")
t = res("time")
ts = sol_undamped.ts

rp1 = get_trajectory(1,linkn_undamped,sol_undamped)

t = res("time")
p1x = res("P1X")
plt.plot(t,p1x)
plt.plot(ts,rp1[1,:])

p1y = res("P1Y")
plt.plot(t,p1y)
plt.plot(ts,rp1[2,:])

p1z = res("P1Z")
plt.plot(t,p1z)
plt.plot(ts,rp1[3,:])


function plot_trajectory(tr)
    @unpack st,traj = tr
    rp7 = get_trajectory(tr,2,2)
    ts = traj.ts
    res = RB.AdamsResults("2AWI3E14H3.xml")

    t = res("time")
    range = 1:1002
    skip = 1:10:1002
    p3x = res("P3X")
    p3y = res("P3Y")
    p3z = res("P3Z")

    fig,_ = plt.subplots(1,1,figsize=(Full_width,5))
    grid = (2,4)
    axs = [plt.subplot2grid(grid,(0,2i-2),1,2,fig=fig) for i = 1:2]
    push!(axs,plt.subplot2grid(grid,(1,1),1,2,fig=fig))
    trange = traj.ts[range]
    tmin = trange[1]
    tmax = trange[end]

    axs[1].plot(ts[range],rp7[1,range])
    axs[1].plot(t[skip],p3x[skip],linestyle="",marker="x")
    axs[1].set_ylabel("Coordinate (m)")
    axs[1].set_xlabel("Time (s)")
    axs[1].set_xlim(tmin,tmax)
    axs[1].set_title("(a) x₇",y=-0.45)
    axs[1].grid(true)

    axs[2].plot(ts[range],rp7[2,range])
    axs[2].plot(t[skip],p3y[skip],linestyle="",marker="x")
    axs[2].set_ylabel("Coordinate (m)")
    axs[2].set_xlabel("Time (s)")
    axs[2].set_xlim(tmin,tmax)
    axs[2].set_title("(b) y₇",y=-0.45)
    axs[2].grid(true)

    axs[3].plot(ts[range],rp7[3,range],label="Zu")
    axs[3].plot(t[skip],p3z[skip],linestyle="",marker="x",label="Adams")
    axs[3].set_ylabel("Coordinate (m)")
    axs[3].set_xlabel("Time (s)")
    axs[3].set_xlim(tmin,tmax)
    axs[3].set_title("(c) z₇",y=-0.45)
    axs[3].grid(true)
    axs[3].legend(loc="upper left", bbox_to_anchor=(1.05, 1.05))

    fig.tight_layout()
    plt.subplots_adjust(wspace=1.0,hspace=0.5)
    fig
end

fig = plot_trajectory(linkn_undamped)
fig.savefig("undamped_trajectory.png",dpi=300,bbox_inches="tight")


linkn_undamped_long, reflinkn = simulate_linkn(;k=3e1,c=0.0,tend=400.0,dt=0.001)

l0s = RB.get_strings_len(reflinkn)
resE = RB.AdamsResults("2B49G.xml")
tE = resE("time")
resE("P2_CM")
function pyplot_energy(tr,l0s;gravity=false)
    @unpack st, traj = tr
    if gravity
        resH = RB.AdamsResults("2BGWI3E3H3.xml")
        resE = RB.AdamsResults("2B49G.xml")
    else
        resH = RB.AdamsResults("2BWI3E14H3.xml")
        resE = RB.AdamsResults("2B49.xml")
    end
    tH = resH("time")
    tE = resE("time")
    ts = traj.ts
    ss = st.strings
    s_index = [2,3,1,5,6,4]

    kes,epes,_,es,es_err = analyse_energy(tr; gravity, elasticity=true)
    @show mean(es_err)
    energy1 = es[1]
    ke_E,epe_H,energy_H,energy_error_H = adams_to_energy(resH,tr,l0s,s_index,energy1;gravity)
    ke_E,epe_E,energy_E,energy_error_E = adams_to_energy(resE,tr,l0s,s_index,energy1;gravity)
    @show mean(energy_error_H)
    @show mean(energy_error_E)
    fig,axs = plt.subplots(1,2,figsize=(Full_width,3))
    axs[1].plot(ts[1:4000:end],es_err[1:4000:end],label="Zu")
    axs[1].plot(tH,energy_error_H,marker="x",linewidth=1,markersize=4,label="Adams")
    axs[1].set_xlim(traj.ts[1],traj.ts[end])
    # axs[1].set_ylim(-0.01,0.01)
    axs[1].set_ylabel("Energy Rel. Err.")
    axs[1].set_xlabel("Time (s)")
    axs[1].set_title("(a)",y=-0.4)
    axs[1].grid(true)

    axs[2].plot(ts[1:4000:end],es_err[1:4000:end],label="Zu")
    axs[2].plot(tE,energy_error_E,marker="x",linewidth=1,markersize=4,label="Adams")
    axs[2].set_xlim(traj.ts[1],traj.ts[end])
    # axs[2].set_ylim(-0.01,0.01)
    axs[2].set_ylabel("Energy Rel. Err.")
    axs[2].set_xlabel("Time (s)")
    axs[2].set_title("(b)",y=-0.4)
    axs[2].legend(loc="upper left",bbox_to_anchor=(1.05, 1.05))
    axs[2].grid(true)

    fig.tight_layout()
    fig
end
fig = pyplot_energy(linkn_undamped_long,l0s)
fig = pyplot_energy(linkn_undamped_long,l0s;gravity=true)
fig.savefig("undamped_energy_error.png",dpi=300,bbox_inches="tight")
fig.savefig("undamped_energy_error_G.png",dpi=300,bbox_inches="tight")
test_slackness(linkn_undamped_long)
plotstructure(linkn_undamped_long)
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
