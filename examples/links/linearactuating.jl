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
h = 0.066
R = RotX(0.0)
c = 1.0
k = 3e1
linkn = links(n,h,R;k,c)
plotstructure(linkn)
q0,q̇0,λ0 = TR.get_initial(linkn)

reflinkn = links(n,0.065,RotY(π/19);k,c)
rl,a = inverse2actuation(linkn,reflinkn)
# n = 1
# function inverse_angle(tg,α,nangle;n,k,c)
#     rls = Vector{Vector{Float64}}()
#     θs = Vector{Float64}()
#     for i = 1:nangle
#         θ = α/nangle*i
#         push!(θs,θ)
#         reflinkn = links(n,0.065,RotY(θ);k,c)
#         rl,a = inverse2actuation(tg,reflinkn)
#         push!(rls,rl)
#     end
#     θs, rls
# end
# θs, rls = inverse_angle(linkn,π/19,3;n,k,c)

# fig, ax = plt.subplots(1,1,figsize=(6,4))
# for j = 1:length(rl)
#     ax.plot(θs, [rls[i][j] for i=1:3])
# end
# fig

# @code_warntype inverse2actuation(linkn;n,h)
plotstructure(reflinkn)
target_t = 10.0
tend = 25.01
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
prob = TS.DyProblem(linearactuate(linkn,q0,a,target_t),q0,q̇0,λ0,(0.0,tend))
sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-11)
plotstructure(linkn,sol,sliderplot)

function set_ax!(ax)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    xmarg = ymarg = 0.001
    zmarg = 0.000000
    zup = 0.35
    x_min = -0.10+xmarg; x_max = 0.10-xmarg
    y_min = -0.10+ymarg; y_max = 0.10-ymarg
    z_min =  0.0+zmarg; z_max = zup-zmarg
    xspan = x_max - x_min
    yspan = y_max - y_min
    zspan = z_max - z_min

    ax.set_xticks([-0.10,0.10])
    ax.set_yticks([-0.10,0.10])
    ax.set_zticks([ 0,0.10,0.20,0.30])
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

    fig = plt.figure(figsize=(6,9))
    axs_raw = [fig.add_subplot(3,2,i,projection="3d") for i = 1:nsps]
    axs = axs_raw

    for (i,ax) in enumerate(axs)
        TR.distribute_q_to_rbs!(tgstruct,sol.qs[steps[i]])
        bars,strings = bars_and_strings_segs_3D(tgstruct)
        ax.add_collection3d(bars)
        ax.add_collection3d(strings)
        ax.set_title("t=$(tstops[i])")
        set_ax!(ax)
    end
    fig.tight_layout()
    fig
end
fig = pyplot_linearactuating(linkn,sol,tstops)
fig.savefig("linearactuating_c=$(c)_t=$tend.png",dpi=300,bbox_inches="tight")
plt.close(fig)

function pyplot_energy(tgstruct,sol)
    kes,epes,_,es,_ = analyse_energy(tgstruct,sol,elasticity=true)
    fig, ax = plt.subplots(3,1,figsize=(6,9))
    ax[1].plot(sol.ts,es,label="Energy")
    ax[1].set_xlim(sol.ts[1],sol.ts[end])
    ax[1].set_xlabel("Time(s)")
    ax[1].set_title("Energy")

    ax[2].plot(sol.ts,kes,label="Energy")
    ax[2].set_xlim(sol.ts[1],sol.ts[end])
    ax[2].set_xlabel("Time(s)")
    ax[2].set_title("Kinetic Energy")

    ax[3].plot(sol.ts,epes,label="Energy")
    ax[3].set_xlim(sol.ts[1],sol.ts[end])
    ax[3].set_xlabel("Time(s)")
    ax[3].set_title("Elastic Potential Energy")
    fig.tight_layout()
    fig
end
fig_energy = pyplot_energy(linkn,sol)
fig_energy.savefig("linearactuating_energy_t=$(target_t)_k=$(k)_c=$c.png",dpi=300,bbox_inches="tight")
plt.close(fig_energy)
