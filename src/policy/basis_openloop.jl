"""
    basis_features(basis::AbstractTimeBasis, t)

Returns a vector of basis function values at time `t`.
"""
function basis_features(basis::AbstractTimeBasis, t) end

"""
    basis_dimension(basis::AbstractTimeBasis)

Returns the number of basis functions (dimension of the feature vector).
"""
function basis_dimension(basis::AbstractTimeBasis) end

"""
    BasisOpenLoop{B <: AbstractTimeBasis, T}

A generic open-loop policy that computes control actions as a linear combination 
of time basis functions.

`u(t) = weights * basis_features(basis, t)`

# Fields
- `weights::Matrix{T}`: Weights matrix where each row corresponds to an actuator/DOF 
                        and each column corresponds to a basis function.
- `basis::B`: The time basis object.
"""
struct BasisOpenLoop{B <: AbstractTimeBasis, T} <: AbstractTimePolicy
    weights::Matrix{T}
    basis::B
end

function update!(policy::BasisOpenLoop) end


# Implement parameter interface for BasisOpenLoop
function get_params(p::BasisOpenLoop)
    vec(p.weights)
end

function set_params!(p::BasisOpenLoop, params)
    p.weights .= reshape(params, size(p.weights))
end

function get_num_of_params(p::BasisOpenLoop)
    length(p.weights)
end

function actuate!(bot::Robot,policy::BasisOpenLoop,inst_state::AbstractCoordinatesState)
    (;structure,hub) = bot
    (;t) = inst_state
    
    # Evaluate basis
    ϕ = basis_features(policy.basis, t)
    
    (;actuators,state,coalition,) = hub
    (;actid2sys_actions_idx) = coalition
    
    # Compute control in-place: state.u = W * ϕ
    mul!(state.u, policy.weights, ϕ)
    foreach(actuators) do actuator
        idx = actid2sys_actions_idx[actuator.id]
        execute!(structure, actuator, (@view state.u[idx]))
    end
end

function vjp_wrt_state(v,policy::BasisOpenLoop,bot::Robot,num_of_actions,solver,solver_state)
    (; tₘ ) = solver_state
    (;structure) = bot
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

function accumulate_param_grad!(grad_storage, policy::BasisOpenLoop, v_total, solver_state, bot)
    (; tₘ ) = solver_state
    
    ϕ = basis_features(policy.basis, tₘ)
    n_basis = length(ϕ)
    action_dim = size(policy.weights, 1)
    
    for j in 1:n_basis
        val = ϕ[j]
        for k in 1:action_dim
            param_idx = (j-1)*action_dim + k
            grad_storage[param_idx] += v_total[k] * val
        end
    end
end

function gen_force_state_jacobian!(∂F∂q, ∂F∂q̇, ∂F∂s, ∂F∂u, policy::AbstractTimePolicy,bot::Robot,inst_state::AbstractCoordinatesState;)
    (;structure, hub) = bot
    (;actuators,coalition) = hub
    (;actid2sys_actions_idx) = coalition
    foreach(actuators) do actuator
        u_idx = actid2sys_actions_idx[actuator.id]
        ∂F∂u_view = @view ∂F∂u[:,u_idx]
        u_view = @view hub.state.u[u_idx]
        gen_force_actu_jacobian!(
            ∂F∂u_view,
            structure,
            actuator,
            u_view
        )
    end
    # Open loop: Force does not depend on state, so ∂F∂q and ∂F∂q̇ are zero
end

# 1. B-Spline Basis

"""
    BSplineBasis{T, B}

A basis using B-splines defined by a knot vector, wrapping BSplineKit.
Uses a single reusable cache to avoid allocations.
"""
struct BSplineBasis{T, B} <: AbstractTimeBasis
    bk_basis::B
    knots::Vector{T}
    ϕcache::Vector{T}
    t_min::T
    t_max::T
    n_basis::Int
end


# 2. ApproxFun Basis Wrapper
"""
    ApproxFunBasis{S} <: AbstractTimeBasis

A wrapper for ApproxFun.jl spaces. 
The dimension of the basis is determined by `order` (number of coefficients).
"""
struct ApproxFunBasis{S} <: AbstractTimeBasis
    space::S
    order::Int # Number of coefficients/basis functions
    # Store specialized domain info if needed, but space usually has it.
end

