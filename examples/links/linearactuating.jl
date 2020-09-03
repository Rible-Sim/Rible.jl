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
# using Printf
using CoordinateTransformations
using Makie
AbstractPlotting.__init__()
# using HomotopyContinuation
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

cd("examples\\links")
include("links_define.jl")
include("link_u_plotting.jl")
include("../analysis.jl")
n = 4
h = 6.6
R = RotX(0.0)
linkn = links(n,h,R,c=50.0)
plotstructure(linkn)
q0,q̇0,λ0 = TR.get_initial(linkn)

rl,a = inverse2actuation(linkn;n=n,h=h)
# @code_warntype inverse2actuation(linkn;n,h)

function linearactuate(tg,q0,target_u,target_t)
    M = TR.build_massmatrix(tg)
    Φ = TR.build_Φ(tg,q0)
    A = TR.build_A(tg)
    # Q̃ = TR.build_Q̃(tg)
    function F!(F,q,q̇,t)
        if t < target_t
            TR.actuate!(tg,t/target_t*target_u)
        end
        TR.reset_forces!(tg)
        TR.distribute_q_to_rbs!(tg,q,q̇)
        TR.update_strings_apply_forces!(tg)
        # TR.apply_gravity!(tg)
        # F .= Q̃*TR.fvector(tg)
        TR.assemble_forces!(F,tg)
    end
    M,Φ,A,F!,nothing
end

dt = 0.01
prob = TS.DyProblem(linearactuate(linkn,q0,a,25.0),q0,q̇0,λ0,(0.0,25.01))
sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-11,verbose=true)
plotstructure(linkn,sol,sliderplot)

function set_ax!(ax)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    xmarg = ymarg = 0.01
    zmarg = 0.000000
    zup = 35
    x_min = -10+xmarg; x_max = 10-xmarg
    y_min = -10+ymarg; y_max = 10-ymarg
    z_min =  0.0+zmarg; z_max = zup-zmarg
    xspan = x_max - x_min
    yspan = y_max - y_min
    zspan = z_max - z_min

    ax.set_xticks([-10,10])
    ax.set_yticks([-10,10])
    ax.set_zticks([ 0,10,20,30])
    ax.set_xlim3d(x_min,x_max)
    ax.set_ylim3d(y_min,y_max)
    ax.set_zlim3d(z_min,z_max)
    ax.set_box_aspect((xspan, yspan, zspan))
    ax.view_init(elev=9, azim=-79)
end

tstops = [0,5,10,15,20,25]
function pyplot_linearactuating(tgstruct,sol,tstops)
    steps = [findfirst((x)->x>i,sol.ts) for i = tstops]
    nsps = length(tstops)

    fig = plt.figure(figsize=(9,6))
    axs_raw = [fig.add_subplot(2,3,i,projection="3d") for i = 1:nsps]
    axs = axs_raw

    for (i,ax) in enumerate(axs)
        TR.distribute_q_to_rbs!(tgstruct,sol.qs[steps[i]])
        bars,strings = bars_and_strings_segs_3D(tgstruct)
        ax.add_collection3d(bars)
        ax.add_collection3d(strings)
        ax.set_title("t=$(tstops[i])")
        set_ax!(ax)
    end
    fig
end
fig = pyplot_linearactuating(linkn,sol,tstops)
fig.savefig("linearactuating.png",dpi=300,bbox_inches="tight")
plt.close(fig)

function pyplot_energy(tgstruct,sol)
    kes,epes,_,es,_ = analyse_energy(tgstruct,sol,elasticity=true)
    fig, ax = plt.subplots(figsize=(9,6))
    ax.plot(sol.ts,es,label="Energy")
    # ax.plot(sol.ts,kes,label="Kinetic Energy")
    # ax.plot(sol.ts,epes,label="Elastic Potential Energy")
    ax.set_xlim(sol.ts[1],sol.ts[end])
    ax.set_xlabel("Time(s)")
    ax.set_ylabel("Energy")
    fig
end
fig_energy = pyplot_energy(linkn,sol)
fig_energy.savefig("linearactuating_energy.png",dpi=300,bbox_inches="tight")
plt.close(fig_energy)
