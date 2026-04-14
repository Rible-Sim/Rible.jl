using Test
import Rible as RB
using Random
# arrays
using LinearAlgebra
using SparseArrays
using StaticArrays
using TypeSortedCollections
using RecursiveArrayTools
using StructArrays
# Optim
# using Optimisers
import FiniteDiff
# data
# using TypedTables
using Rotations
using ComponentArrays
# visualize/plot
import GeometryBasics as GB
using MarchingCubes
# print
# 
# using Printf
# using TexTables
using Logging

using FileIO
using MeshIO

include("../../examples/robots/pointmass.jl")

@testset "Point Mass" begin
    @testset "Bouncing" begin
        #mark ground plane
        θ=0.0
        ground_plane = RB.StaticContactSurfaces(
            [
                RB.HalfSpace([-tan(θ),0,1],zeros(3)),
            ]
        )

        # time
        tspan = (0.0,1.5)
        h = 1e-3

        # parameters and initial conditions
        restitution_coefficients = [0.5]
        v0s = [1.0]

        # pointmass
        pm = new_pointmass(;
            e = restitution_coefficients[1],
            μ=0.1,
            origin_velocity = [v0s[1],0,0]
        )

        # frictional contact dynamics problem 
        prob = RB.DynamicsProblem(
            pm;
            env=ground_plane,
            contact_model=RB.RestitutionFrictionCombined(
                RB.NewtonRestitution(),
                RB.CoulombFriction(),
            )
        )

        # simulate
        sol = RB.solve!(
            prob,
            RB.DynamicsSolver(
                RB.Zhong06(),
                RB.InnerLayerContactSolver(
                    RB.InteriorPointMethod()
                )
            );
            tspan,dt=h,
            ftol=1e-14,
            maxiters=50,
            exception=false
        )
        r1 = RB.get_trajectory!(pm,1,1).u[end]
        @test r1[1] ≈ 0.6038355644971142
        ṙ1 = RB.get_velocity!(pm,1,1).u[end]
        @test ṙ1[3] ≈  -0.0033338635879290764
    end

    @testset "Sliding" begin
        #mark inclined plane
        θ = 15 |> deg2rad
        inclined_plane = RB.StaticContactSurfaces(
            [
                RB.HalfSpace([-tan(θ),0,1],zeros(3)),
            ]
        )

        # initial conditions
        origin_position = [0.0,0,-1e-7]
        origin_velocity = [2.0cos(θ),0,2.0sin(θ)]

        # simulation
        tspan = (0.0,0.6)
        μ = 0.3
        pm = new_pointmass(;e=0.0, μ, origin_position, origin_velocity)

        prob = RB.DynamicsProblem(
            pm;
            env=inclined_plane,
            contact_model=RB.RestitutionFrictionCombined(
                RB.NewtonRestitution(),
                RB.CoulombFriction(),
            )
        )
        sol = RB.solve!(
            prob,
            RB.DynamicsSolver(
                RB.Zhong06(),
                RB.InnerLayerContactSolver(
                    RB.InteriorPointMethod()
                )
            );
            tspan,dt=1e-3,ftol=1e-14,maxiters=50,exception=false
        )
        r1 = RB.get_trajectory!(pm,1,1).u[end]
        @test r1 ≈ [0.35895776877373514, 0.0, 0.09618234426055731]
        ṙ1 = RB.get_velocity!(pm,1,1).u[end]
        @test ṙ1[1] ≈ -0.0001165620176319937
    end
end
