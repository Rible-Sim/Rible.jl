using Test
import Rible as RB
using Random
using Test
# arrays
using LinearAlgebra
using SparseArrays
using StaticArrays
using CircularArrays
using OffsetArrays
using BlockDiagonals
using TypeSortedCollections
using RecursiveArrayTools
using Interpolations
using StructArrays
# Optim
## using Optimisers
using ForwardDiff
import FiniteDiff
# data
# using TypedTables
using DataStructures
using Rotations
using CoordinateTransformations
using EponymTuples
using IterTools
using Unitful
using ComponentArrays
# visualize/plot
import GeometryBasics as GB
using Makie
using LaTeXStrings
using MarchingCubes
import Meshes
using Match
# print
# 
# using Printf
# using TexTables
using LoggingExtras
using Logging
using LoggingFormats

using FileIO
using MeshIO
include("../../examples/robots/spinningtop.jl")

@testset "Spinning Top" begin
    # solver
    solver = RB.DynamicsSolver(
        RB.Zhong06(),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        ),
    )

    @testset "First scenario" begin
        # Contact Surfaces
        planes = RB.StaticContactSurfaces(
            [
                RB.HalfSpace([0,0,1.0],[0,0,0.0]),
            ]
        )

        # Initial conditions
        origin_position = [0,0,0.5]
        R = RotX(0.0)
        origin_velocity = [1.0,0.0,0.0]
        Ω = [0.0,0.0,200.0]
        # parameters
        μ = 0.95
        e = 0.5
        # time 
        tspan = (0.0,1.8)
        h = 1e-4

        top = make_top(origin_position,R,origin_velocity,Ω,:NCF;μ,e,loadmesh=true)
        RB.solve!(
            RB.DynamicsProblem(
                top;
                env=planes,
                contact_model=RB.RestitutionFrictionCombined(
                    RB.NewtonRestitution(),
                    RB.CoulombFriction(),
                )
            ),
            solver;
            tspan,
            dt=h,
            ftol=1e-14,
            maxiters=50,exception=false,
        )
        r1 = RB.get_trajectory!(top,1,1).u[end]
        @test r1 ≈  [0.7125525120016437, 0.05778527088190157, 0.04286749115165942]
        ṙ1 = RB.get_velocity!(top,1,1).u[end]
        @test ṙ1 ≈  [2.3383344859230775, -7.523336388327795, -3.677831743208304]
    end

    @testset "Sliding simulation" begin
        planes = RB.StaticContactSurfaces(
            [
                RB.HalfSpace([0,0,1.0],[0,0,0.0]),
            ]
        )
        R = RotX(π/18)
        origin_position = [0,0,0.037]
        origin_velocity = [1.0,0.0,0.0]
        Ω = [0,0,50.0]
        dts = [1e-3,1e-2]
        tops = [
            begin
            bot = make_top(origin_position,R,origin_velocity,Ω,:NCF; μ = 0.01, e = 0.0,loadmesh=true)
            RB.solve!(
                RB.DynamicsProblem(
                    bot;
                    env=planes,
                    contact_model=RB.RestitutionFrictionCombined(
                        RB.NewtonRestitution(),
                        RB.CoulombFriction(),
                    )
                ),
                RB.DynamicsSolver(
                    RB.Zhong06(),
                    RB.InnerLayerContactSolver(
                        RB.InteriorPointMethod()
                    );
                    checkpersist
                );
                tspan=(0.0,2.0),
                dt,ftol=1e-14,maxiters=50,exception=false,
            )
            bot
            end 
            for dt in dts, checkpersist in [true,false]
        ]
        r1 = RB.get_trajectory!(tops[1],1,1).u[end]
        @test r1 ≈ [1.768646474661639, 0.014902561508855354, 0.045264943681779794]
        ṙ1 = RB.get_velocity!(tops[1],1,1).u[end]
        @test ṙ1 ≈ [0.20471706378997628, -2.103694849182572, 0.4033763668782808] 
    end
end
