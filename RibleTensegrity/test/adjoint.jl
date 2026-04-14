#mark deps
#using TestEnv
#TestEnv.activate()
using Rible
import Rible as RB
using RibleTensegrity
import RibleTensegrity as RT
using Test
using FiniteDiff
using LinearAlgebra
using StaticArrays
using Logging
using LoggingExtras
using Rotations
using StaticArrays
using TypeSortedCollections
using CircularArrays
using Random

include("../../examples/bodies/rigidbar.jl")
include("../../examples/robots/superball.jl")

include("../../examples/bodies/make_3d_bar.jl")
include("../examples/robots/five_bars.jl")

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
#mark superball
@testset "Superball Gradient Tests" begin
    integrator = RB.Zhong06()
    test_tol = 1e-3
    
    # Parameters
    l = 1.7/2
    d = l/2
    
    # Contact surfaces
    flatplane_ballbot = RB.StaticContactSurfaces(
        [
            RB.HalfSpace([0,0,1.0],[0,0,0.0]),
        ]
    )
    
    # Common setup parameters
    μ = 0.9
    e = 0.8
    z0 = l^2/(sqrt(5)*d) + 2.0
    ω = SVector(0.0,0.0,0.0)
    tend = 0.27
    tspan = (0.0,tend)
    h = 1e-2
    
    # Objective function
    objt_ballbot = RB.Objective(
        [0], #trajectory_gauges_weights
        [0], #trajectory_actuators_weights
        [1], #terminal_gauges_weights
        [0], #terminal_actuators_weights
        (x)->1
    )
    
    # Test cases with different initial velocities
    test_cases = [
        (name = "Standard rolling", 
         v0 = MVector(7.0, 2.0, -7.0)),
        
        (name = "Fast rolling", 
         v0 = MVector(10.0, 3.0, -10.0)),
        
        (name = "Slow rolling", 
         v0 = MVector(5.0, 1.0, -5.0)),
        
        (name = "Diagonal motion",
         v0 = MVector(6.0, 6.0, -8.0))
    ]

    for test_case in test_cases[[1]]
        @testset "$(test_case.name)" begin
            v0 = test_case.v0
            
            # Create ballbot with specific initial velocity
            ballbot = superball(
                0.0;
                origin_velocity = SVector(v0...),
                ω = ω,
                μ = μ,
                e = e,
                l, d,
                z0 = z0,
                visible = true,
                constrained = false,
                loadmesh = false,
            )
            
            # Setup gradient transformation matrices
            nq = RB.get_num_of_full_coords(ballbot.structure)
            N = RB.intrinsic_nullspace(ballbot.structure)
            ∂ξ∂ẏ0 = zeros(Int, 30, 3)
            for i = 1:6
                ∂ξ∂ẏ0[((1:3).+5(i-1)),:] .= I(3)
            end
            ∂q̇0∂ẏ0 = N*∂ξ∂ẏ0
            ∂q0∂y0 = zeros(nq, 0)
            
            # Create cost function
            ball_cost = RB.make_cost(
                ballbot, RB.NoPolicy(), flatplane_ballbot, integrator, objt_ballbot, RB.InitialCost();
                idx = nq .+ collect(1:nq),
                params_func = (x) -> repeat([x...,0,0,0], outer=6),
                tspan, dt=h, ftol=1e-10, maxiters=50, exception=false
            )
            
            # Compute cost value
            cost_val = ball_cost(v0)
            @test isfinite(cost_val)
            println("  Cost value: $cost_val")
            
            # Compute finite difference gradient
            G_fd = FiniteDiff.finite_difference_gradient(ball_cost, v0)
            
            # Create analytical gradient function
            ball_cost_gradient_wrt_initial = RB.make_cost_gradient(
                ballbot, RB.NoPolicy(), flatplane_ballbot, integrator, objt_ballbot, RB.InitialCost();
                ∂q̇0∂ẏ0, ∂q0∂y0,
                idx = nq .+ collect(1:nq),
                params_func = (x) -> repeat([x...,0,0,0], outer=6),
                tspan, dt=h, ftol=1e-10, exception=false
            )
            
            # Compute analytical gradient
            G = zeros(length(v0))
            ball_cost_gradient_wrt_initial(G, v0)
            
            # Test gradient properties
            @test length(G) == length(G_fd)
            @test all(isfinite.(G))
            @test all(isfinite.(G_fd))
            
            # Test gradient accuracy
            max_error = maximum(abs.(G .- G_fd))
            println("  Max gradient error: $max_error")
            println("  Analytical gradient: $G")
            println("  Finite diff gradient: $G_fd")
            println("  Difference: $(G .- G_fd)")
            
            @test max_error < test_tol
        end
    end
end

#mark five bars
# Reproduces Dopico D, Zhu Y, Sandu A, et al. Direct and Adjoint Sensitivity Analysis of Ordinary Differential Equation Multibody Formulations[J]. Journal of Computational and Nonlinear Dynamics, 2014, 10(1).
@testset "Five Bars Gradient Tests" begin
    integrator = RB.Zhong06()
    fb = five_bars()
    c = RB.get_params(fb.structure)

    objt_fb = RB.Objective(
        [1],   # trajectory_gauges_weights
        [0, 0], # trajectory_actuators_weights
        [0],   # terminal_gauges_weights
        [0, 0], # terminal_actuators_weights
        t -> 1
    )

    tspan = (0.0, 5.0)
    dt = 1e-3
    idx = [24, 27]

    fb_cost = RB.make_cost(
        fb, RB.NoPolicy(), RB.GravityEnv(), integrator, objt_fb, RB.StructureCost();
        tspan, dt,
        idx
    )

    fb_cost_gradient = RB.make_cost_gradient(
        fb, RB.NoPolicy(), RB.GravityEnv(), integrator, objt_fb, RB.StructureCost();
        tspan, dt,
        idx
    )

    G = zeros(length(idx))
    fb_cost_gradient(G, c[idx])
    G .*= 2

    G_fd = FiniteDiff.finite_difference_gradient(fb_cost, c[idx]) .* 2

    max_error = maximum(abs.(G .- G_fd))
    @test max_error < 1e-3

    target_G = [-4.227926063668897, 3.2110065945391253]
    @test isapprox(G, target_G; rtol=1e-4, atol=1e-6)
end


# v2p2 = RB.get_velocity!(fb,2,2)
# fig = Figure()
# ax = Makie.Axis(fig[1,1])
# lines!(ax,fb.traj.t,v2p2[1,:],label="x")
# lines!(ax,fb.traj.t,v2p2[3,:],label="y")
# axislegend(ax)
# xlims!(ax,extrema(fb.traj.t)...)
# fig

# ME = RB.mechanical_energy!(fb,RB.Gravity())
# fig = Figure()
# ax = Makie.Axis(fig[1,1])
# lines!(ax,fb.traj.t,ME.V.-ME.V[begin],label="V")
# lines!(ax,fb.traj.t,ME.T,label="T")
# axislegend(ax)
# xlims!(ax,extrema(fb.traj.t)...)
# fig

#mark contact five bars
@testset "Five Bars Contact Gradient Tests" begin
    integrator = RB.Zhong06()

    objt_fb = RB.Objective(
        [1],    # trajectory_gauges_weights
        [0, 0], # trajectory_actuators_weights
        [0],    # terminal_gauges_weights
        [0, 0], # terminal_actuators_weights
        t -> 1
    )

    θ = 0.0
    ground_plane = RB.StaticContactSurfaces(
        [
            RB.HalfSpace([-tan(θ), 0, 1], [0, 0, -2.35]),
        ]
    )

    fb_song = five_bars(;
        restitution_coefficient = 0.0,
        friction_coefficient = 0.05,
        err_ref = [-0.5, 0.0, -2.0]
    )

    tspan = (0.0, 2.0)
    dt = 1e-3
    idx = [24, 27]
    params = RB.get_params(fb_song.structure)

    fb_song_contact_cost = RB.make_cost(
        fb_song, RB.NoPolicy(), ground_plane, integrator, objt_fb, RB.StructureCost();
        tspan, dt,
        ftol = 1e-10,
        idx
    )

    G_fd = FiniteDiff.finite_difference_gradient(
        fb_song_contact_cost, params[idx]
    )

    fb_song_contact_cost_gradient = RB.make_cost_gradient(
        fb_song, RB.NoPolicy(), ground_plane, integrator, objt_fb, RB.StructureCost();
        tspan, dt,
        ftol = 1e-10,
        idx
    )

    G = zeros(length(idx))
    fb_song_contact_cost_gradient(G, params[idx])

    @test all(isfinite.(G))
    @test all(isfinite.(G_fd))

    max_error = maximum(abs.(G .- G_fd))
    @test max_error < 1e-3
end



# v2p2 = RB.get_velocity!(fb_song,2,2)
# r2p2 = RB.get_trajectory!(fb_song,2,2)
# fig = Figure()
# ax1 = Makie.Axis(fig[1,1])
# lines!(ax1,fb_song.traj.t,v2p2[3,:],label="ẏ")
# ax2 = Makie.Axis(fig[1,2])
# lines!(ax2,fb_song.traj.t,r2p2[3,:],label="y")
# axislegend(ax1)
# axislegend(ax2)
# xlims!(ax1,extrema(fb_song.traj.t)...)
# xlims!(ax2,extrema(fb_song.traj.t)...)
# fig
