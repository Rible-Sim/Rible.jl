using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using Parameters
import PyPlot; const plt = PyPlot
plt.pygui(true)
using LaTeXStrings
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


# Δu
# u0 = [s.original_restlen for s in linkn.strings]
# rl = u0 + Δu
# ℓ = [s.state.length for s in linkn.strings]
# ℓ - rl
# s = 1 ./ℓ
# TR.actuate!(linkn,Δu)


function set_ax!(ax)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    xmarg = ymarg = 0.0001
    zmarg = 0.000000
    zup = 0.20
    x_min = -0.10+xmarg; x_max = 0.10-xmarg
    y_min = -0.10+ymarg; y_max = 0.10-ymarg
    z_min =  0.0+zmarg; z_max = zup-zmarg
    xspan = x_max - x_min
    yspan = y_max - y_min
    zspan = z_max - z_min

    ax.set_xticks([-0.10,0,0.10])
    ax.set_yticks([-0.10,0,0.10])
    ax.set_zticks([ 0,0.10,0.20])
    ax.set_xlim3d(x_min,x_max)
    ax.set_ylim3d(y_min,y_max)
    ax.set_zlim3d(z_min,z_max)
    ax.set_box_aspect((xspan, yspan, zspan))
    ax.view_init(elev=9, azim=-79)
end

function eigenanalysis(θs)
    n = length(θs)
    nstage = 2
    h = 0.066
    R = RotY(0.0)
    linkn = links(nstage,h,R;k=3e1)
    Y = Array(build_Y(linkn))

    ωs = Vector{Vector{Float64}}()
    Zs = Vector{Matrix{Float64}}()
    rls = Vector{Vector{Float64}}()
    epes = Vector{Float64}()
    # ω0,Z0 = TR.undamped_eigen(linkn,q0,λ0)

    # push!(ωs,ω0)
    # push!(Zs,Z0)

    fig = plt.figure(figsize=(8,5))
    axs = [fig.add_subplot(2,ceil(Int,n/2),i,projection="3d") for i = 1:n]
    θstrings = [latexstring("\\theta=0")]
    for i = 1:n-1
        raco = 1//12*i/(n-1)
        raco_n = numerator(raco)
        raco_d = denominator(raco)
        if raco_n == 1
            push!(θstrings, latexstring("\\theta=\\pi/$raco_d"))
        else
            push!(θstrings, latexstring("\\theta=$raco_n\\pi/$raco_d"))
        end
    end
    show(θstrings)
    for i = 1:n
        iθ = θs[i]
        R = RotY(iθ)
        reflinkn = links(nstage,h,R)
        refqi,_,_ = TR.get_initial(reflinkn)
        refλi,_,ai= TR.inverse(linkn,reflinkn,Y,gravity=false)
        actlinkn = deepcopy(linkn)
        TR.actuate!(actlinkn,ai)
        rli = [s.state.restlen for s in actlinkn.strings]
        push!(rls,rli)
        TR.reset_forces!(actlinkn)
        TR.distribute_q_to_rbs!(actlinkn,refqi)
        TR.update_strings_apply_forces!(actlinkn)
        @show transpose(TR.build_A(actlinkn)(refqi))*refλi ≈ TR.build_Q̃(actlinkn)*TR.fvector(actlinkn)

        epei = TR.elastic_potential_energy(actlinkn,refqi)
        push!(epes,epei)
        ax = axs[i]
        bars,strings = bars_and_strings_segs_3D(actlinkn)
        ax.add_collection3d(strings)
        ax.add_collection3d(bars)
        ax.set_title(θstrings[i])
        set_ax!(ax)

        ωi,Zi = TR.undamped_eigen(actlinkn,refqi,refλi)
        push!(ωs,ωi)
        push!(Zs,Zi)
    end
    fig.tight_layout()
    ωs,Zs,rls,epes,fig
end
function gen_angles(α1 = 0.0, α2 = π/12, n = 6)
    α = α2 - α1
    α1 .+ α/n.*collect(0:n)
end

θs5 = gen_angles(0.0, π/12, 4)
θs6 = gen_angles(0.0, π/12, 5)
θs7 = gen_angles(0.0, π/12, 6)
θs = θs6
θs = sort(union(θs5,θs6,θs7))

ωs, Zs, rls, epes, fig = eigenanalysis(θs)
fig.savefig("link2_inv.png",dpi=300,bbox_inches="tight")

n = length(θs)
fig, ax = plt.subplots(1,2,figsize=(8,3))
ax[1].plot(θs,[rls[i][1] for i=1:n],marker="o",ls="--",label=L"s_1")
ax[1].plot(θs,[rls[i][2] for i=1:n],marker="o",ls="--",label=L"s_2")
ax[1].plot(θs,[rls[i][3] for i=1:n],marker="o",ls="--",label=L"s_3")
ax[1].set_xlabel(L"\theta")
ax[1].set_ylabel("Rest Length (m)")

ax[2].plot(θs,[rls[i][4] for i=1:n],marker="o",ls="--",label=L"s_4")
ax[2].plot(θs,[rls[i][5] for i=1:n],marker="o",ls="--",label=L"s_5")
ax[2].plot(θs,[rls[i][6] for i=1:n],marker="o",ls="--",label=L"s_6")
ax[2].set_xlabel(L"\theta")
ax[2].set_ylabel("Rest Length (m)")

ax[1].legend()
ax[2].legend()
fig.tight_layout()
fig.savefig("link2_inv_restlen.png",dpi=300,bbox_inches="tight")


fig, ax = plt.subplots(1,1,figsize=(4,3))
ax.plot(θs,epes,marker="o",ls="--")
ax.set_xlabel(L"\theta")
ax.set_ylabel("Potential Energy")
fig.tight_layout()
fig.savefig("link2_inv_potential.png",dpi=300,bbox_inches="tight")

fig, ax = plt.subplots(1,1,figsize=(6,3))
ax.plot(θs,[ωs[i][1] for i=1:n],marker="o",fillstyle="none",ls="--",label="Mode 1")
ax.plot(θs,[ωs[i][2] for i=1:n],marker="o",fillstyle="none",ls="--",label="Mode 2")
ax.plot(θs,[ωs[i][3] for i=1:n],marker="o",fillstyle="none",ls="--",label="Mode 3")
ax.set_xlabel(L"\theta")
ax.set_ylabel("Frequency")
ax.grid(true)
ax.legend(loc="upper right", bbox_to_anchor=(1.4, 1.0))
fig.tight_layout()
fig.savefig("link2_inv_frequency.png",dpi=300,bbox_inches="tight")

fig, ax = plt.subplots(1,1,figsize=(6,3))
