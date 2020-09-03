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
    xmarg = ymarg = 0.01
    zmarg = 0.000000
    zup = 20
    x_min = -10+xmarg; x_max = 10-xmarg
    y_min = -10+ymarg; y_max = 10-ymarg
    z_min =  0.0+zmarg; z_max = zup-zmarg
    xspan = x_max - x_min
    yspan = y_max - y_min
    zspan = z_max - z_min

    ax.set_xticks([-10,0,10])
    ax.set_yticks([-10,0,10])
    ax.set_zticks([ 0,10,20])
    ax.set_xlim3d(x_min,x_max)
    ax.set_ylim3d(y_min,y_max)
    ax.set_zlim3d(z_min,z_max)
    ax.set_box_aspect((xspan, yspan, zspan))
    ax.view_init(elev=9, azim=-79)
end

function eigenanalysis(n)
    nstage = 2
    h = 6.6
    R = RotY(0.0)
    linkn = links(nstage,h,R)
    # return linkn
    Y = Array(build_Y(linkn))
    q0,_,_ = TR.get_initial(linkn)
    λ0,_,a = TR.inverse(linkn,deepcopy(linkn),Y,gravity=false)
    TR.actuate!(linkn,a)
    TR.reset_forces!(linkn)
    TR.update_strings_apply_forces!(linkn)

    @show transpose(TR.build_A(linkn)(q0))*λ0 ≈ TR.build_Q̃(linkn)*TR.fvector(linkn)

    ωs = Vector{Vector{Float64}}()
    Zs = Vector{Matrix{Float64}}()

    ω0,Z0 = TR.undamped_eigen(linkn,q0,λ0)

    push!(ωs,ω0)
    push!(Zs,Z0)

    fig = plt.figure(figsize=(9,6))
    axs = [fig.add_subplot(2,3,i,projection="3d") for i = 1:n]
    α = π/15
    θ = α/(n-1).*collect(0:n-1)
    θstrings = [latexstring("\\theta=0")]
    for i = 1:n-1
        raco = 1//15*i/(n-1)
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
        iθ = θ[i]
        R = RotY(iθ)
        reflinkn = links(nstage,h,R)
        refqi,_,_ = TR.get_initial(reflinkn)
        refλi,_,ai= TR.inverse(linkn,reflinkn,Y,gravity=false)
        actlinkn = deepcopy(linkn)
        TR.actuate!(actlinkn,ai)
        TR.reset_forces!(actlinkn)
        TR.distribute_q_to_rbs!(actlinkn,refqi)
        TR.update_strings_apply_forces!(actlinkn)

        @show transpose(TR.build_A(actlinkn)(refqi))*refλi ≈ TR.build_Q̃(actlinkn)*TR.fvector(actlinkn)

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
    ωs,Zs,fig
end
ωs, Zs, fig = eigenanalysis(6)
fig.savefig("eigenplot.png",dpi=300,bbox_inches="tight")
