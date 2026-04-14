#mark deps
#TestEnv.activate(); using Revise;
using Rible
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
# import CairoMakie as CM

include("../examples/robots/pointmass.jl");

include("../examples/robots/spinningtop.jl")

# using AppleAccelerate
cd(@__DIR__)
figdir::String = "figures"
isdir(figdir) || mkdir(figdir)
isfile("private.jl") && include("private.jl") 

tw::Float64 = 455.24411 #src
scalefactor::Float64 = 2 #nb

# coordsType= :QCF
coordsType= :NCF

try
    @eval using AbbreviatedStackTraces;
    @eval ENV["JULIA_STACKTRACE_ABBREVIATED"] = true;
    @eval ENV["JULIA_STACKTRACE_MINIMAL"] = true;
catch e
    @warn "error while importing AbbreviatedStackTraces" e
end

global_logger(ConsoleLogger(stdout, Logging.Warn; show_limited=false))
#mark point mass tests
@testset "Point Mass Gradient Tests" begin
    integrator = RB.Zhong06()
    test_tol = 1e-4
    @testset "Ground Plane - Gradient Accuracy" begin
        # Setup ground plane
        θ = 0.0
        ground_plane_pm = RB.StaticContactSurfaces(
            [RB.HalfSpace([-tan(θ),0,1],zeros(3))]
        )
        
        # Setup point mass
        tend = 1.5
        rovo_pm = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0]
        
        objt_pm = RB.Objective(
            [0], [0], [1], [0], (x)-> 1 
        )
        
        pm = new_pointmass(;
            e = 0.5,
            μ = 0.1,
            origin_position = rovo_pm[1:3],
            origin_velocity = rovo_pm[4:6],
            term_position = [1.0, 0.0, 0.0]
        )
        
        # Create cost function
        pm_cost = RB.make_cost(
            pm, RB.NoPolicy(), ground_plane_pm, integrator, objt_pm, RB.InitialCost();
            tspan = (0.0, tend), dt = 1e-3, ftol = 1e-10
        )
        pm_cost(rovo_pm)
        # Compute finite difference gradient
        G_fd = FiniteDiff.finite_difference_gradient(pm_cost, rovo_pm)
        
        # Compute analytical gradient
        pm_cost_gradient_wrt_initial = RB.make_cost_gradient(
            pm, RB.NoPolicy(), ground_plane_pm, integrator, objt_pm, RB.InitialCost();
            tspan = (0.0, tend), dt = 1e-3, ftol = 1e-10
        )
        G = zeros(length(rovo_pm))
        pm_cost_gradient_wrt_initial(G, rovo_pm)
        
        # Test gradient accuracy
        @test length(G) == length(G_fd)
        @test all(isfinite.(G))
        max_error = maximum(abs.(G .- G_fd))
        @test max_error < test_tol
        @test G ≈ [-0.39616442439125876,0.0,-1.7629079067660314e-6,-0.44746520567044484,0.0,-0.0335455131425537]
        println("  Ground plane gradient max error: $max_error")
        
        # Test combined cost and gradient function
        pm_cost_and_gradient = RB.make_cost_and_gradient(
            pm, RB.NoPolicy(), ground_plane_pm, integrator, objt_pm, RB.InitialCost();
            tspan = (0.0, tend), dt = 1e-3, ftol = 1e-10
        )
        
        # Test cost only
        cost_val = pm_cost_and_gradient(1, nothing, rovo_pm)
        @test isfinite(cost_val)
        @test cost_val ≈ 0.07847312557816533
        # Test gradient only
        G2 = zeros(length(rovo_pm))
        pm_cost_and_gradient(nothing, G2, rovo_pm)
        @test maximum(abs.(G2 .- G)) < test_tol
        
        # Test both cost and gradient
        cost_val2 = pm_cost_and_gradient(1, G2, rovo_pm)
        @test cost_val ≈ cost_val2
        @test maximum(abs.(G2 .- G)) < test_tol
    end
    
    @testset "Inclined Plane - Gradient Accuracy" begin
        # Setup inclined plane
        θ = 15 |> deg2rad
        inclined_plane = RB.StaticContactSurfaces(
            [RB.HalfSpace([-tan(θ), 0, 1], zeros(3))]
        )
        
        # Setup point mass on inclined plane
        origin_position = [0.0, 0.0, -1e-7]
        origin_velocity = [2.0cos(θ), 0.0, 2.0sin(θ)]
        rovo_inclined_pm = vcat(origin_position, origin_velocity)
        
        inclined_pm = new_pointmass(;
            e = 0.0,
            μ = 0.3,
            origin_position,
            origin_velocity
        )
        
        objt_inclined_pm = RB.Objective(
            [0], [0], [1], [0], (x)-> 1 
        )
        
        tend = 0.6
        inclined_pm_cost = RB.make_cost(
            inclined_pm, RB.NoPolicy(), inclined_plane, integrator, objt_inclined_pm, RB.InitialCost();
            tspan = (0.0, tend), dt = 1e-3, ftol = 1e-10
        )
        
        # Compute finite difference gradient
        G_fd = FiniteDiff.finite_difference_gradient(inclined_pm_cost, rovo_inclined_pm)
        
        # Compute analytical gradient
        inclined_pm_cost_gradient = RB.make_cost_gradient(
            inclined_pm, RB.NoPolicy(), inclined_plane, integrator, objt_inclined_pm, RB.InitialCost();
            tspan = (0.0, tend), dt = 1e-3, ftol = 1e-10
        )
        G = zeros(length(rovo_inclined_pm))
        inclined_pm_cost_gradient(G, rovo_inclined_pm)
        
        # Test gradient accuracy
        max_error = maximum(abs.(G .- G_fd))
        @test length(G) == length(G_fd)
        @test all(isfinite.(G))
        @test G ≈ [-2.1410422575185066,0.0,-0.9038176566709808,-0.776788788858103,0.0,-0.46933430798293735]
        @test isapprox(G[4:6],G_fd[4:6];atol=1e-3)
        @test_broken max_error < test_tol
        
        println("  Inclined plane gradient max error: $max_error")
    end
end

#mark top
@testset "Spinning Top Gradient Tests" begin
    integrator = RB.Zhong06()
    test_tol = 1e-3
    
    # Setup contact surfaces
    planes_top = RB.StaticContactSurfaces(
        [
            RB.HalfSpace([0,0,1.0],[0,0,0.0]),
        ]
    )
    
    # Common setup for top
    z0 = 0.5
    origin_position = [0,0,z0]
    R = RotX(0.0)
    origin_velocity = [1.0,0.0,0.0]
    Ω = [0.0,0.0,200.0]
    μ = 0.95
    e = 0.5
    tend = 1.8
    tspan = (0.0,tend)
    h = 1e-3
    
    # Create the top
    top = make_top(origin_position,R,origin_velocity,Ω,:NCF;μ,e,loadmesh=true)
    
    # Setup objective
    objt_top = RB.Objective(
        [0], #trajectory_gauges_weights
        [0], #trajectory_actuators_weights
        [1], #terminal_gauges_weights
        [0], #terminal_actuators_weights
        (x)->1
    )
    
    # Setup gradient transformation matrices
    nq̌ = RB.get_num_of_free_coords(top.structure)
    ∂q0∂y0 = zeros(Int,nq̌,3)
    ∂q0∂y0[1:3,:] .= I(3)
    N = RB.intrinsic_nullspace(top.structure, top.structure.state.system)
    ∂ξ∂ẏ0 = zeros(Int,6,3)
    ∂ξ∂ẏ0[1:3,:] .= I(3)
    ∂q̇0∂ẏ0 = N*∂ξ∂ẏ0
    
    # Test cases with different initial conditions
    test_cases = [
        (name = "Initial position", 
         origin_pos_vel = [0.0, 0.0, 0.5, 1.0, 0.0, 0.0]),
        
        (name = "Offset position", 
         origin_pos_vel = [0.0, 0.0, 0.5+2*0.01897941-1e-5, 1.0, 0.0, 0.0]),
        
        (name = "Perturbed state",
         origin_pos_vel = [-0.011668778495605695, 
                          -0.0004510174113496357, 
                          0.8100598351899901, 
                          0.994208802590762, 
                          -0.00021221362266814233, 
                          0.06565512964076882])
    ]
    
    for test_case in test_cases
        @testset "$(test_case.name)" begin
            origin_pos_vel = test_case.origin_pos_vel
            
            # Create cost function
            top_cost = RB.make_cost(top, RB.NoPolicy(), planes_top, integrator, objt_top, RB.InitialCost();
                idx = [1,2,3,13,14,15],    
                tspan, dt=h, ftol=1e-10, maxiters=50, exception=false,
            )
            
            # Compute cost value
            cost_val = top_cost(origin_pos_vel)
            @test isfinite(cost_val)
            println("  Cost value: $cost_val")
            
            # Compute finite difference gradient
            G_fd = FiniteDiff.finite_difference_gradient(top_cost, origin_pos_vel)
            
            # Create analytical gradient function
            top_cost_gradient_wrt_initial = RB.make_cost_gradient(
                top, RB.NoPolicy(), planes_top, integrator, objt_top, RB.InitialCost();
                idx = [1,2,3,13,14,15],
                ∂q̇0∂ẏ0, ∂q0∂y0, 
                tspan, dt=h, ftol=1e-10, maxiters=50, exception=false,
            )
            
            # Compute analytical gradient
            G = zeros(length(origin_pos_vel))
            top_cost_gradient_wrt_initial(G, origin_pos_vel)
            
            # Test gradient properties
            @test length(G) == length(G_fd)
            @test all(isfinite.(G))
            @test all(isfinite.(G_fd))
            
            # Test gradient accuracy
            max_error = maximum(abs.(G[[1,2,5]] .- G_fd[[1,2,5]]))
            println("  Max gradient error: $max_error")
            println("  Analytical gradient: $G")
            println("  Finite diff gradient: $G_fd")
            println("  Difference: $(G .- G_fd)")
            
            @test max_error < test_tol
        end
    end
end

# Additional debug/development code (commented out - requires manual setup)
# tend = 0.404
# origin_pos_vel_test = [0.0, 0.0, 0.5, 1.0, 0.0, 0.0]
# get_J(tend,)(origin_pos_vel_test)
# import FiniteDifferences
# G_fd = FiniteDiff.finite_difference_gradient(get_J(tend,), origin_pos_vel_test[3:3],)
# G_fd = FiniteDifferences.grad(FiniteDifferences.central_fdm(5, 1), get_J(tend,), origin_pos_vel_test[3:3])[1]
# G_fd .- G

# fd_data = StructArray(
#     [
#         (J=0.017994401560807883, rovo=[0.788363869724089]),
#         (J=0.01861668439551888, rovo=[0.7955958582127227]),
#         (J=0.0192759371291437, rovo=[0.8028278467013564]),
#         (J=0.019030941696306595, rovo=[0.8100598351899901]),
#         (J=0.018337603299557324, rovo=[0.8172918236786239]),
#         (J=0.017230082192872302, rovo=[0.8245238121672576]),
#         (J=0.016174862743719928, rovo=[0.8317558006558913]),
#         (J=0.019034420509015666, rovo=[0.8100382382918155]),
#         (J=0.019032681044357867, rovo=[0.8100490367409029]),
#         (J=0.019030941696306595, rovo=[0.8100598351899901]),
#         (J=0.019495655284648833, rovo=[0.8100706336390774]),
#         (J=0.01949388462219041, rovo=[0.8100814320881647]),
#     ]
# )

# scatter(VectorOfArray(fd_data.rovo)[1,:],fd_data.J,)
# with_logger(FileLogger("debug.log")) do
#     withlevel(Debug; verbosity = 4) do
#         get_dJ!(tend, )(G,origin_pos_vel_test)
#     end
# end


