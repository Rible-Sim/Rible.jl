using Test
using Rible
import Rible as RB
using RibleQCF
import RibleQCF as QCF
using Random
using Test
# arrays
using LinearAlgebra
using SparseArrays
using StaticArrays
using TypeSortedCollections
using RecursiveArrayTools
using StructArrays
# Optim
## using Optimisers
using ForwardDiff
import FiniteDiff
# data
using Rotations
using EponymTuples
# visualize/plot
import GeometryBasics as GB
using MarchingCubes
# print
using Logging
using LoggingExtras
using FileIO
using MeshIO


try
    @eval using AbbreviatedStackTraces;
    @eval ENV["JULIA_STACKTRACE_ABBREVIATED"] = true;
    @eval ENV["JULIA_STACKTRACE_MINIMAL"] = true;
catch e
    @warn "error while importing AbbreviatedStackTraces" e
end

include("../../examples/robots/spinningtop.jl")

@testset "QCF Spinning Top" begin
    # solver
    solver = RB.DynamicsSolver(
        RB.Zhong06(),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        ),
    )

    #mark First scenario
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
    h = 1e-2

    top = make_top(origin_position,R,origin_velocity,Ω,:QCF;μ,e,loadmesh=true) #src
    body1 = RB.get_bodies(top.structure)[1]
    @test has_constant_mass_matrix(body1.coords) === Val(false)
    @test has_constant_mass_matrix(body1) === Val(false)
    @test has_constant_mass_matrix(top.structure.bodies) === Val(false)
    @test has_constant_mass_matrix(top.structure) === Val(false)
    # top = make_top(origin_position,R,origin_velocity,Ω,:NCF;μ,e,loadmesh=true);
    RB.solve!(
        RB.DynamicsProblem(
            top,
            env=RB.GravityEnv(),
        ),
        solver;
        tspan = (0.0,20h),
        dt=h,
        ftol=1e-14,
        maxiters=50,exception=true,
        # verbose_contact=false,
        # use_fd_jacobian = false,
    );
    @test top.traj.p[end] ≈  [0.5838707000000011, 0.0, -1.1455543134000115, 3.5002561553365563e-7, 0.0, 0.0, 0.12082799999949255]
    me = RB.mechanical_energy!(top, RB.Gravity())

    @test me.E[end] ≈ 9.197221133499966

end
