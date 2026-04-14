using Test
import Rible as RB
using RibleTensegrity
import RibleTensegrity as RT
using LinearAlgebra
using SparseArrays
using StaticArrays
using CircularArrays
import CircularArrays as CA

using TypeSortedCollections
using LoggingExtras
using Rotations
using Unitful
using ElasticArrays
# IO
using FileIO, MeshIO
import GeometryBasics as GB
using EponymTuples
using StructArrays


include("../../examples/bodies/rigidbar_nonsmooth_repro.jl")
include("../../examples/robots/superball.jl")

@testset "Superball" begin
    # Parameters
    l = 1.7/2
    d = l/2

    # Env
    flatplane = RB.StaticContactSurfaces(
        [
            RB.HalfSpace([0,0,1.0],[0,0,0.0]),
        ]
    )

    # Solver
    solver =  RB.DynamicsSolver(
        RB.Zhong06(),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        ),
    )

    @testset "First Scenario (rolling)" begin
        ballbot = superball(
            0.0;
            origin_velocity = SVector(7.0,2.0,-7.0),
            ω = SVector(0.0,0.0,0.0),
            μ = 0.9,
            e = 0.8,
            l,d,
            z0 = l^2/(sqrt(5)*d) + 2.0,
            visible = true,
            constrained = false,
            loadmesh = false,
        )

        # time 
        tspan = (0.0,5.0)
        h = 1e-2

        RB.solve!(
            RB.DynamicsProblem(ballbot;
                env=flatplane,
                contact_model=RB.RestitutionFrictionCombined(
                    RB.NewtonRestitution(),
                    RB.CoulombFriction(),
                )
            ),
            solver;
            tspan,dt=h,ftol=1e-12,maxiters=200,exception=false
        )
        r1 = RB.get_trajectory!(ballbot,1,1).u[end]
        @test r1 ≈ [18.194562968348812, 6.184938174507589, 1.6237835187624035] atol=1e-4
        ṙ1 = RB.get_velocity!(ballbot,1,1).u[end]
        @test ṙ1 ≈ [5.072678933978139, 2.5827950142353484, -0.7467500018097193] atol=1e-4
    end

    @testset "Second Scenario" begin
        ballbot = superball(
            0.0;
            origin_velocity = SVector(2.0,1.0,0),
            ω = SVector(0.0,0.0,1.0),
            μ = 0.05,
            e = 0.0,
            l,d,
            z0 = l^2/(sqrt(5)*d) - 1e-3,
            loadmesh = false,
            constrained = false,
        )

        # time simulations 
        tspan = (0.0,5.0)
        h = 1e-2
        prob = RB.DynamicsProblem(ballbot;
            env=flatplane,
            contact_model=RB.RestitutionFrictionCombined(
                RB.NewtonRestitution(),
                RB.CoulombFriction(),
            )
        )
        RB.solve!(prob,solver;tspan,dt=h,ftol=1e-14,exception=false)
        r1 = RB.get_trajectory!(ballbot,1,1).u[end]
        @test r1 ≈ [4.1654406913447986, 2.3553290347107425, 1.4576338694859803]
        ṙ1 = RB.get_velocity!(ballbot,1,1).u[end]
        @test ṙ1 ≈ [0.01884024343612825, -0.015969393416154276, -0.044877879712985694]
    end
end
