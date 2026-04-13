using Test
using Rible
import Rible as RB
using ApproxFun
using LinearAlgebra

# --- Constructors for specific Spaces ---

"""
    ChebyshevBasis(domain_range, order)

Creates a Chebyshev basis on `domain_range` (e.g., `0.0..1.0` or `(0.0, 1.0)`).
Order is the number of coefficients (degree + 1).
"""
function ChebyshevBasis(domain_range, order::Int)
    # Ensure ApproxFun is loaded
    s = Chebyshev(domain_range)
    return RB.ApproxFunBasis(s, order)
end

"""
    FourierBasisApprox(domain_range, order)

Creates a Fourier basis using ApproxFun. Note: `order` is total coefficients.
"""
function ApproxFunFourierBasis(domain_range, order::Int)
    s = Fourier(domain_range)
    return RB.ApproxFunBasis(s, order)
end

"""
    TaylorBasis(domain_range, order::Int)

Creates a Taylor basis (monomials).
Note: We use the canonical `Taylor()` space (Hardy space on unit circle) which corresponds
to monomials `z^k`. We ignores `domain_range` to avoid conformal mapping that would
warp the time axis. This ensures `u(t) = sum c_k t^k`.
"""
function TaylorBasis(domain_range, order::Int)
    s = Taylor()
    return RB.ApproxFunBasis(s, order)
end


@testset "ApproxFun Policies" begin
    # Test Chebyshev Basis
    @testset "ChebyshevBasis" begin
        # Domain [0, 1], Order 5 (T0 to T4)
        basis = ChebyshevBasis(0.0..1.0, 5)
        # Should be dimension 5
        @test RB.basis_dimension(basis) == 5
        
        # Create a policy with weights selection only the first basis (T0 = 1)
        # Weights: 1 x 5
        W = zeros(1, 5)
        W[1, 1] = 1.0 # T0 constant term
        policy = RB.BasisOpenLoop(W, basis)
        
       
        
        # Select T1 (scaled axis)
        # on 0..1, T1(x) maps x->2x-1
        # at x=0.5 -> 2(0.5)-1 = 0
        W2 = zeros(1, 5)
        W2[1, 2] = 1.0
        policy2 = RB.BasisOpenLoop(W2, basis)
    end

    @testset "FourierBasis" begin
        # ApproxFun Fourier default domain is -π..π but we can specify
        basis = ApproxFunFourierBasis(0.0..6.28, 5) # ~0..2π
        @test RB.basis_dimension(basis) == 5
        
        W = zeros(1, 5)
        W[1, 1] = 1.0 # usually constant 1
        policy = RB.BasisOpenLoop(W, basis)
        
        
        # Test sin term
        # Fourier space order: 1, sin(1x), cos(1x), ...
        # depends on ApproxFun implementation.
        # Let's verify by checking value.
        W_sin = zeros(1, 5)
        W_sin[1, 2] = 1.0
        policy_sin = RB.BasisOpenLoop(W_sin, basis)
        
        # We can also check basis_features if implemented
        feats = RB.basis_features(basis, 1.0)
        @test length(feats) == 5
        
    end
    
    @testset "TaylorBasis" begin
        basis = TaylorBasis(0.0..1.0, 3) # 1, x, x^2
        @test RB.basis_dimension(basis) == 3
        
        W = zeros(1, 3)
        W[1, 2] = 1.0 # x term
        policy = RB.BasisOpenLoop(W, basis)
        
    end
end
