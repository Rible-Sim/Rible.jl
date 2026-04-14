using Test
using Rible
import Rible as RB
import BSplineKit
using LinearAlgebra

@testset "BSplineBasis Tests" begin
    # Define time points
    t0 = 0.0
    tf = 1.0
    
    # 1. Test BSplineBasis creation
    knots = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0] # Cubic spline knots
    degree = 3
    basis = RB.BSplineBasis(knots, degree)
    
    @test RB.basis_dimension(basis) == length(knots) - degree - 1 # 9 - 3 - 1 = 5 basis functions
    
    # 2. Test basis_features at specific points
    # t = 0.0
    # First basis function should be 1.0 at start? 
    # For clamped knots, first basis is 1 at t=0.
    
    phi_0 = RB.basis_features(basis, 0.0)
    @test length(phi_0) == RB.basis_dimension(basis)
    @test abs(phi_0[1] - 1.0) < 1e-6
    @test sum(phi_0) ≈ 1.0 # Partition of unity
    
    # t = 0.5
    phi_mid = RB.basis_features(basis, 0.5)
    @test sum(phi_mid) ≈ 1.0
    
    # t = 1.0 (End)
    phi_end = RB.basis_features(basis, 1.0)
    @test sum(phi_end) ≈ 1.0
    @test abs(phi_end[end] - 1.0) < 1e-6

    # 3. Test BasisOpenLoop with BSplineBasis
    dim = RB.basis_dimension(basis) # 5
    n_actuators = 2
    weights = rand(n_actuators, dim)
    
    policy = RB.BasisOpenLoop(weights, basis)
    
    # Manual verify
    phi_val = RB.basis_features(basis, 0.5)
    u_expected = weights * phi_val

end
