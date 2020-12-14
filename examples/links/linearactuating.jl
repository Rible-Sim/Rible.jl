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
# using Printf
using Interpolations
using NLsolve
using CoordinateTransformations
# using HomotopyContinuation
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

cd("examples\\links")
include("links_define.jl")
include("link_u_plotting.jl")
include("link_inverse.jl")
include("../analysis.jl")

function find_actuation(epe0 = 0.05, epef = 0.10, c=0.0)
    n = 4
    h = 0.066
    linkn = links(n,h,RotX(0.0);k=3e1,c)
    Y = build_Y(linkn)
    reflinkn = links(n,h,RotY(π/12);k=3e1)
    _,_,a0 = inverse_with_energy_3(linkn,deepcopy(linkn),Y,fill(epe0,3))
    _,_,af = inverse_with_energy_3(linkn,reflinkn,Y,fill(epef,3))
    linkn, a0, af
end
epe0,epef = 0.05,0.05
epe0,epef = 0.10,0.05
epe0,epef = 0.05,0.10
c = 0.0
linkn,a0,af = find_actuation(epe0,epef)
linkn,a0,af = find_actuation(epe0,epef,5.0)
# TR.actuate!(linkn,a0)
# TR.reset_forces!(linkn)
# TR.update_strings_apply_forces!(linkn)
# TR.elastic_potential_energy(linkn)

q0,q̇0,λ0 = TR.get_initial(linkn)

function actuation_interpolation(target_t,initial_a, target_a)
    itps = [interpolate([a0,af],
            BSpline(Quadratic(Flat(OnGrid())))) for (a0,af) in zip(initial_a,target_a)]
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
# aitp = actuation_interpolation(target_t,[1,2],[3,4])
# ts = LinRange(0,tend,100)
# plt.plot(ts,aitp.(ts))

target_t = 10.0
tend = 20.01

function linearactuate(tg, q0, aitp)
    M = TR.build_massmatrix(tg)
    Φ = TR.build_Φ(tg,q0)
    A = TR.build_A(tg)
    # Q̃ = TR.build_Q̃(tg)
    function F!(F,q,q̇,t)
        TR.actuate!(tg,aitp(t))
        TR.reset_forces!(tg)
        TR.distribute_q_to_rbs!(tg,q,q̇)
        TR.update_strings_apply_forces!(tg)
        TR.assemble_forces!(F,tg)
    end
    M,Φ,A,F!,nothing
end

dt = 0.01
aitp = actuation_interpolation(target_t,a0,af)
prob = TS.DyProblem(linearactuate(linkn,q0,aitp),q0,q̇0,λ0,(0.0,tend))
sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-14)
plotstructure(linkn,sol,sliderplot)

function pyplot_energy(tgstruct,sol,as)
    kes,epes,_,es,_ = analyse_energy(tgstruct,sol,as;elasticity=true)
    fig, ax = plt.subplots(1,2,figsize=(8,4))
    # ax[1].plot(sol.ts,es,label="Energy")
    # ax[1].set_xlim(sol.ts[1],sol.ts[end])
    # ax[1].set_xlabel("Time(s)")
    # ax[1].set_title("Energy")

    ax[1].plot(sol.ts,kes)
    ax[1].set_xlim(sol.ts[1],sol.ts[end])
    ax[1].set_xlabel("Time (s)")
    ax[1].set_ylabel("Energy (J)")
    ax[1].set_title("(a) Kinetic Energy", y = -0.3)

    ax[2].plot(sol.ts,epes)
    ax[2].set_xlim(sol.ts[1],sol.ts[end])
    ax[2].set_ylim(0.04*3,0.06*3)
    ax[2].set_xlabel("Time (s)")
    ax[2].set_ylabel("Energy (J)")
    ax[2].set_title("(b) Elastic Potential Energy", y = -0.3)
    fig.tight_layout()
    fig
end
# function actuation_sequence(a0,af,ts,target_t)
#     target_t_index = findfirst((x)->x>target_t,ts)
#     as_target = LinRange(a0,af,target_t_index)
#     as_rest = LinRange(af,af,length(ts)-target_t_index)
#     as = vcat(as_target,as_rest)
# end
as = aitp.(sol.ts)
fig_energy = pyplot_energy(linkn,sol,aitp.(sol.ts))
fig_energy.savefig("linearactuating_energy_t=$(target_t)_epe0=$(epe0)_epef=$(epef).png",dpi=300,bbox_inches="tight")
plt.close(fig_energy)


tstops = [0,5,10]
function pyplot_linearactuating(tgstruct,sol,tstops)
    steps = [findfirst((x)->x>i,sol.ts) for i = tstops]
    nsps = length(tstops)

    fig = plt.figure(figsize=(12,5))
    axs_raw = [fig.add_subplot(1,3,i,projection="3d") for i = 1:nsps]
    axs = axs_raw
    alphabet = join('a':'z')
    for (i,ax) in enumerate(axs)
        TR.distribute_q_to_rbs!(tgstruct,sol.qs[steps[i]])
        bars,strings = bars_and_strings_segs_3D(tgstruct)
        ax.add_collection3d(bars)
        ax.add_collection3d(strings)
        ax.set_title("($(alphabet[i])) t=$(tstops[i]) (s)",y=-0.12)
        set_ax!(ax,0.32;elev=11,azim=-76)
    end
    fig.tight_layout()
    fig
end
fig = pyplot_linearactuating(linkn,sol,tstops)
fig.savefig("linearactuating_c=$(c)_t=$tend.png",dpi=300,bbox_inches="tight")
plt.close(fig)
