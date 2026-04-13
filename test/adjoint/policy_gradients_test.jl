using Test
using Rible
import Rible as RB
using ApproxFun
using BSplineKit
using FiniteDiff
using LinearAlgebra
using Logging
using Random
using Rotations
using StaticArrays
using FileIO
using TypeSortedCollections
using ComponentArrays

# Include the cart pole robot definition
include("../../examples/robots/cartpole.jl")

try
    @eval using AbbreviatedStackTraces;
    @eval ENV["JULIA_STACKTRACE_ABBREVIATED"] = true;
    @eval ENV["JULIA_STACKTRACE_MINIMAL"] = true;
catch e
    @warn "error while importing AbbreviatedStackTraces" e
end


# Helper functions to create policies

function get_bspline_policy(capta_dim, action_dim; num_knots=10, degree=3, rng=Random.default_rng())
    knots = collect(range(0.0, 1.0, length=num_knots))
    basis = RB.BSplineBasis(knots, degree; augment_knots=true)
    dim = RB.basis_dimension(basis)
    weights = randn(rng, action_dim, dim) * 0.1
    return RB.BasisOpenLoop(weights, basis)
end

function get_approxfun_policy(capta_dim, action_dim; order=10, rng=Random.default_rng())
    # Ensure ApproxFun is loaded
    s = Chebyshev(Interval(0.0, 1.0))
    basis = RB.ApproxFunBasis(s, order)
    weights = randn(rng, action_dim, order) * 0.1
    return RB.BasisOpenLoop(weights, basis)
end

function get_linear_policy(capta_dim, action_dim; rng=Random.default_rng())
    K = randn(rng, action_dim, capta_dim) * 0.1
    extractor = RB.GaugesExtractor([1.0, 1.0])
    featurizer = RB.IdentityFeaturizer()
    return RB.LinearFeedbackPolicy(K; extractor, featurizer)
end

function get_pd_policy(capta_dim, action_dim; rng=Random.default_rng())
    K = randn(rng, action_dim, capta_dim) * 0.01
    x_des = randn(rng, capta_dim)
    extractor = RB.GaugesExtractor([1.0, 1.0])
    featurizer = RB.IdentityFeaturizer()
    return RB.DiscreteProportionalDerivativePolicy(K, x_des; extractor, featurizer)
end


# Suppress logging output during tests
global_logger(ConsoleLogger(stdout, Logging.Error; show_limited=false))

# Setup
integrator = RB.Zhong06()
coordsType = :NCF

# Create cart pole robot
ctp = cart_pole(; θ0 = π / 4, y0 = -1.0, coordsType)
num_of_capta = RB.get_num_of_capta(ctp)
num_of_actions = 1  # Cart pole has 1 actuator (force on cart)

# Define objective
objt_ctp = RB.Objective(
    [0,0,1,1,0,0],  # trajectory_gauges_weights 
    [1],  # trajectory_actuators_weights
    [0,0,1,1,0,0],  # terminal_gauges_weights
    [0],  # terminal_actuators_weights
    _ -> 1, 
)

env = RB.GravityEnv()

# Test parameters
tspan = (0.0, 1.0)
dt = 1e-2

# Define policies to test (excluding RBF as requested)
open_loop_policies = [
    ("BSpline", get_bspline_policy(num_of_capta, num_of_actions, num_knots=10, degree=3, rng=Random.MersenneTwister(4321))),
    ("ApproxFun-Chebyshev", get_approxfun_policy(num_of_capta, num_of_actions, order=10, rng=Random.MersenneTwister(4321))),
]

# Define policies to test
feedback_policies = [
    ("Linear", get_linear_policy(num_of_capta, num_of_actions, rng=Random.MersenneTwister(4321))),
    ("PD", get_pd_policy(num_of_capta, num_of_actions, rng=Random.MersenneTwister(1234))),
]

@testset "Cart pole policies gradient tests" begin
    test_tol = 1e-5
    # Test each policy
    for (plc_name, policy) in open_loop_policies
        @testset "$plc_name Policy Gradients" begin
            # Create cost and gradient function
            ctp_cost_and_gradient = RB.make_cost_and_gradient(
                ctp, policy, env, integrator, objt_ctp, RB.ControlCost();
                tspan = tspan, dt = dt,
            )
            
            params = RB.get_params(policy)
            
            # Define cost function for finite difference
            function cost_func(p_vec)
                ctp_cost_and_gradient(1, nothing, p_vec)
            end
            
            # Compute cost (ensure it works)
            cost = cost_func(params)
            @test isfinite(cost)
            
            # Compute finite difference gradient
            G_fd = FiniteDiff.finite_difference_gradient(cost_func, params)
            
            # Compute analytic gradient
            G = zeros(length(params))
            ctp_cost_and_gradient(nothing, G, params)
            
            # Test gradients are finite
            @test all(isfinite.(G))
            @test all(isfinite.(G_fd))
            
            # Test gradient accuracy
            # Use relaxed tolerance for ApproxFun policies
            tol = plc_name == "ApproxFun-Chebyshev" ? 1e-3 : test_tol
            max_error = maximum(abs.(G .- G_fd))
            @test max_error < tol
            println("$plc_name max_error: $max_error")
        end
    end

    test_tol = 1e-4 # Relaxed tolerance slightly for generic Lux

    # Test each policy
    for (plc_name, policy) in feedback_policies
        @testset "$plc_name Policy Gradients" begin
            # Create cost and gradient function
            ctp_cost_and_gradient = RB.make_cost_and_gradient(
                ctp, policy, env, integrator, objt_ctp, RB.ControlCost();
                tspan = tspan, dt = dt,
            )
            
            params = RB.get_params(policy)
            
            # Define cost function for finite difference
            function cost_func(p_vec)
                ctp_cost_and_gradient(1, nothing, p_vec)
            end
            
            # Compute cost (ensure it works)
            cost = cost_func(params)
            @test isfinite(cost)
            
            # Compute finite difference gradient
            G_fd = FiniteDiff.finite_difference_gradient(cost_func, params)
            
            # Compute analytic gradient
            G = zeros(length(params))
            ctp_cost_and_gradient(nothing, G, params)
            
            # Test gradients are finite
            @test all(isfinite.(G))
            @test all(isfinite.(G_fd))
            
            # Test gradient accuracy
            @test G ≈ G_fd atol = 1e-7 rtol = test_tol
            max_error = maximum(abs.(G .- G_fd))
            println("$plc_name max_error: $max_error")
        end
    end
end