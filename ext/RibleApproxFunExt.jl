
module RibleApproxFunExt

using LinearAlgebra
using ForwardDiff
using SparseArrays
import Rible as RB
using ApproxFun

import Rible: Robot, AbstractCoordinatesState, AbstractTimeBasis, basis_dimension, basis_features, BasisOpenLoop
import Rible: actuate!, execute!, vjp_wrt_state, accumulate_param_grad!, get_num_of_full_coords, get_num_of_cstr, get_num_of_aux_var

import Rible: ApproxFunBasis

function basis_dimension(basis::ApproxFunBasis)
    return basis.order
end

# But we can implement basis_features for completeness or fallback
function basis_features(basis::ApproxFunBasis, t)
    # Evaluates the first `order` basis functions at t
    # This might be slow if called repeatedly without specialization
    # ApproxFun doesn't always expose "evaluate k-th basis" cheaply without creating a Fun
    # A generic way:
    vals = zeros(eltype(t), basis.order)
    for k in 1:basis.order
        # Create a coefficient vector with 1 at k
        coeffs = zeros(basis.order)
        coeffs[k] = 1.0
        f = Fun(basis.space, coeffs)
        vals[k] = f(t)
    end
    return vals
end

# BasisOpenLoop{<:ApproxFunBasis}
function actuate!(bot::Robot,policy::BasisOpenLoop{<:ApproxFunBasis},inst_state::AbstractCoordinatesState)
    (;structure,hub) = bot
    (;t) = inst_state
    (;actuators,state,coalition,) = hub
    (;actid2sys_actions_idx) = coalition
    
    for i in 1:size(policy.weights, 1)
        # Evaluate basis
        f = Fun(policy.basis.space, policy.weights[i, :])
        
        # Compute control: u = W * ϕ
        # weights is (action_dim x n_basis)
        state.u[i] = f(t)
    end
    
    
    foreach(actuators) do actuator
        idx = actid2sys_actions_idx[actuator.id]
        execute!(structure, actuator, (@view state.u[idx]))
    end
end

function vjp_wrt_state(v,policy::BasisOpenLoop{<:ApproxFunBasis},bot::Robot,num_of_actions,solver,solver_state)
    (; tₘ ) = solver_state
    (; structure) = bot
    
    nq = get_num_of_full_coords(structure)
    nλ = get_num_of_cstr(structure)
    ns = get_num_of_aux_var(structure)
    
    ∂ϕ∂qₖᵀ = spzeros(eltype(tₘ), nq)
    ∂ϕ∂qₖ₊₁ᵀ = spzeros(eltype(tₘ), nq)
    ∂ϕ∂pₖᵀ = spzeros(eltype(tₘ), nq)
    ∂ϕ∂pₖ₊₁ᵀ = spzeros(eltype(tₘ), nq)
    ∂ϕ∂λᵀ = spzeros(eltype(tₘ), nλ)
    ∂ϕ∂sₖᵀ = spzeros(eltype(tₘ), ns)
    ∂ϕ∂sₖ₊₁ᵀ = spzeros(eltype(tₘ), ns)
    
    ∂ϕ∂qₖᵀ, ∂ϕ∂qₖ₊₁ᵀ, ∂ϕ∂pₖᵀ, ∂ϕ∂pₖ₊₁ᵀ, ∂ϕ∂λᵀ, ∂ϕ∂sₖᵀ, ∂ϕ∂sₖ₊₁ᵀ
end

function accumulate_param_grad!(grad_storage, policy::BasisOpenLoop{<:ApproxFunBasis}, v_total, solver_state, bot)
    (; tₘ ) = solver_state
    
    function paramfun(a)
        f = Fun(policy.basis.space, a); 
        f(tₘ)
    end
    
    action_dim, n_basis = size(policy.weights)
    
    for k in 1:action_dim
        # Parameter index for W_{k,j} (col-major)
        param_idx = ((k-1)*n_basis + 1) : (k*n_basis)
        grad_storage[param_idx] .+= v_total[k] .* ForwardDiff.gradient(paramfun, policy.weights[k,:])
    end
end

end