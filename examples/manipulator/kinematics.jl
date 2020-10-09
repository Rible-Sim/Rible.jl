using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using BenchmarkTools
import PyPlot; const plt = PyPlot
plt.pygui(true)
using LaTeXStrings
# using NLsolve
# using DifferentialEquations
using Makie
AbstractPlotting.__init__()
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot
cd("examples/manipulator")
include("man_define.jl")
include("man_plotting.jl")

k = 1.e4; c = 0.0
function inversekinematics(θs;k,c)
    n = length(θs)
    ndof = 1
    manipulator = man_ndof(ndof;k,c)
    as = Vector{Vector{Float64}}()

    Y = build_Y(manipulator)


    for θ in θs
        refman = man_ndof(ndof;θ,k,c) # reference
        refqi,_,_ = TR.get_initial(refman)
        refλi,_,a= TR.inverse(manipulator,refman,Y)
        push!(as,a)
    end
    as
end


α = 0.9
n = 10
θs = [-i*α/n for i = 0:n]

as = inversekinematics(θs; k, c)
fig,ax = plt.subplots(1,1,num = 1,figsize=(4,3))
x = abs.(θs)
y = [a for a in as]
ax.plot(x,y,marker="o",ls="--")
ax.set_xlabel(L"$\theta$ (Rad)")
ax.set_ylabel("Rest length (cm)")
ax.set_xlim(0,1.0)
ax.set_ylim(0,9.0)
ax.grid(true)
fig.savefig("man_inverse.png",dpi=300,bbox_inches="tight")

ndof = 6
manipulator = man_ndof(ndof; k, c)
ncoords = manipulator.ncoords
function eigenanalysis(θs; k, c)
    ndof = 6
    manipulator = man_ndof(ndof; k, c)

    Y = build_Y(manipulator)

    ωs = Vector{Vector{Float64}}()
    Zs = Vector{Matrix{Float64}}()

    for θ in θs
        refman = man_ndof(ndof;θ,k,c) # reference
        refqi,_,_ = TR.get_initial(refman)
        refλi,_,_= TR.inverse(manipulator,refman,Y)
        actmani = deepcopy(manipulator)

        ωi,Zi,_ = TR.undamped_eigen(actmani,refqi,refλi)
        M = TR.build_massmatrix(actmani)
        TR.normalize_wrt_mass!(Zi,M)
        push!(ωs,ωi)
        push!(Zs,Zi)
    end
    ωs,Zs
end

ωs,Zs = eigenanalysis(θs; k, c)
ys = [[ω[i] for ω in ωs] for i = 1:length(ωs[1])]

function plotmodes(x,ys)
    fig,ax = plt.subplots(1,1,num = 2,figsize=(4,3))
    for (i,y) in enumerate(ys)
        ax.plot(x,y,marker="o",fillstyle="none",ls="--",label="Mode $i")
    end
    ax.set_xlabel("Rest length (cm)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_xlim(0,8)
    ax.set_ylim(0,8)
    ax.grid(true)
    ax.legend(loc="upper right", bbox_to_anchor=(1.4, 1.0))
    fig
end
fig = plotmodes(y,ys)
fig.savefig("man_mode.png",dpi=300,bbox_inches="tight")

fig,ax = plt.subplots(1,1,figsize=(4,3))
npos = 11
nmode = 2
for i = 1:nmode
# for i = 1:2
    ax.plot(1:ncoords,Zs[npos][1:ncoords,i],label = "Mode $i")
    ax.set_xticks(collect(1:ncoords))
end
ax.legend()
