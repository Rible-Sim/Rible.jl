
module RibleBSplineKitExt

using LinearAlgebra
using SparseArrays
import Rible as RB
using BSplineKit

import Rible: Robot, AbstractCoordinatesState, AbstractTimeBasis, basis_dimension, basis_features, BasisOpenLoop
import Rible: actuate!, execute!, vjp_wrt_state, accumulate_param_grad!, get_num_of_full_coords, get_num_of_cstr, get_num_of_aux_var

import Rible: BSplineBasis



function BSplineBasis(knots::Vector{T}, degree::Int; augment_knots=false) where T
    # Create BSplineKit basis with exact knots provided
    # We use augment=Val(false) to prevent BSplineKit from adding extra knots
    # This assumes the user provided a full knot vector or intends to.
    
    order = degree + 1
    # BSplineKit takes ownership of knots? copy to be safe
    # If augment_knots is requested, we could pass Val(true), but default to explicit control.
    augment = augment_knots ? Val(true) : Val(false)
    bk_basis = BSplineKit.BSplineBasis(order, copy(knots); augment=augment)
    n_basis = length(bk_basis)
    ϕcache = zeros(T, n_basis)
    t_min, t_max = BSplineKit.boundaries(bk_basis)
    return BSplineBasis{T, typeof(bk_basis)}(bk_basis, knots, ϕcache, t_min, t_max, n_basis)
end

function basis_dimension(basis::BSplineBasis)
    return basis.n_basis
end

function basis_features(basis::BSplineBasis, t)
    # Reuses a single cache vector; mutates and returns it.
    ϕ = basis.ϕcache
    fill!(ϕ, zero(eltype(ϕ)))
    t_min = basis.t_min
    t_max = basis.t_max
    (t < t_min || t > t_max) && return ϕ

    # Handle right endpoint using limit from left (if strictly needed) or rely on BSplineKit
    # BSplineKit intervals are typically [t_i, t_{i+1}).
    # If t == t_max, we might be 0 unless we handle it.
    if t == t_max
        t = t - 1e-9 * max(1.0, abs(t))
    end
    
    # BSplineKit.evaluate(basis, i, t) returns the value of i-th basis function at t.
    # To be efficient, we should only evaluate non-zero ones, but iterating all is safe for now.
    # Optimization: Find active interval and only evaluate relevant support.
    
    # Using 'evaluate' per basis index
    for i in eachindex(ϕ)
        # evaluate returns 0 if outside support usually, or might throw if out of domain?
        # We checked domain above.
        ϕ[i] = BSplineKit.evaluate(basis.bk_basis, i, t)
    end
    
    return ϕ
end

end

