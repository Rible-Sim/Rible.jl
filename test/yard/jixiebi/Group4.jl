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
# using Meshes
cd("docs/examples/jixiebi")
includet("define.jl")
includet("dyfun.jl")
includet("../vis.jl")

delta = 350
bot1 = build_jixiebi(3, [0, delta, delta],ifreverse = false)
bot2 = build_jixiebi(3, [0, -delta, delta], ifreverse = false)
bot3 = build_jixiebi(3, [0, delta, -delta],ifreverse = true)
bot4 = build_jixiebi(3, [0, -delta, -delta], ifreverse = true)


prob1 = TR.SimProblem(bot1,dynfuncs1)
prob2 = TR.SimProblem(bot2,dynfuncs1)
prob3 = TR.SimProblem(bot3,dynfuncs1)
prob4 = TR.SimProblem(bot4,dynfuncs1)

dt = 1e-2; T = 6.0

TR.solve!(prob1, TR.FBZhong06(),
        (
            actuate! = actuate_for_3!,
            prescribe! = nothing
            );
        dt=dt,tspan=(0.0,T),ftol=1e-7,verbose=true)
TR.solve!(prob2, TR.FBZhong06(),
        (
            actuate! = actuate_for_2!,
            prescribe! = nothing
            );
        dt=dt,tspan=(0.0,T),ftol=1e-7,verbose=true)
TR.solve!(prob3, TR.FBZhong06(),
        (
            actuate! = actuate_for_2!,
            prescribe! = nothing
            );
        dt=dt,tspan=(0.0,T),ftol=1e-7,verbose=true)
TR.solve!(prob4, TR.FBZhong06(),
        (
            actuate! = actuate_for_3!,
            prescribe! = nothing
            );
        dt=dt,tspan=(0.0,T),ftol=1e-7,verbose=true)



with_theme(theme_pub;
    Axis3 = (
        # azimuth = -2π/3,
        # elevation = -2π/3,
        azimuth = -2.1043951023931977,
        elevation = 0.7902036732051045
        # azimuth = -pi/2,
        # elevation =  - pi/4
    ),
    Mesh = (
        # color = :black,
        transparency = false,
    ),
    ) do
plot_traj_with_jixiepintai!(bot1,bot2,bot3,bot4;
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
    xlims=(-1500,4000),
    ylims=(-2000,2000),
    zlims=(-2000,2000),
    cylinder_x = 1800,
    cylinder_r = 300,
)
end



with_theme(theme_pub;
    Axis3 = (
        # azimuth = -2π/3,
        # elevation = -2π/3,
        azimuth = -2.1043951023931977,
        elevation = 0.7902036732051045
        # azimuth = -pi/2,
        # elevation =  - pi/4
    ),
    Mesh = (
        # color = :black,
        transparency = false,
    ),
    ) do
plot_traj!(botc1;
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
    xlims=(-1500,4000),
    ylims=(-2000,2000),
    zlims=(-2000,2000),
)
end