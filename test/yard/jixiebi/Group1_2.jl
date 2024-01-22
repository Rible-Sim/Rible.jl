using LinearAlgebra
using SparseArrays
using StaticArrays
using BenchmarkTools
using FileIO, MeshIO
using Printf
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using Unitful
using NLsolve
using Revise
using StructArrays
using EponymTuples
using LaTeXStrings
using StaticArrays
using GeometryBasics
using Rotations
using CoordinateTransformations
using Meshing
using TypeSortedCollections
using Match
import TensegrityRobots as TR
# cd("docs/examples/jixiebi")
includet("define.jl")
includet("dyfun.jl")
includet("../vis.jl")

bot = build_jixiebi2(4)
prob1 = TR.SimProblem(bot,dynfuncs1)
dt = 1e-2; T = 7.0

TR.solve!(prob1,TR.FBZhong06(),
        (
            actuate! = actuate12!,
            prescribe! = nothing
            );
        dt=dt,tspan=(0.0,T),ftol=1e-7,verbose=true)


with_theme(theme_pub;
    Axis3 = (
        azimuth = -π/2,
        elevation = 0.0,
    ),
    Mesh = (
        # color = :black,
        transparency = false,
    ),
    ) do
plot_traj!(bot;
    # AxisType=LScene,
    # doslide=false,
    doslide=true,
    showinfo=true,
    axtitle=false,
    # hidezdecorations=true,
    hideydecorations=true,
    fontsize=10,
    # gridsize=(2,2),
    # attimes=[0.0,2.25,3.75,6.0],
    # atsteps=nothing,
    showground=false,
    showlabels=false,
    rigidcolor=:grey,
    xlims=(-100,3400),
    ylims=(-800,800),
    zlims=(-800,800),
)
end