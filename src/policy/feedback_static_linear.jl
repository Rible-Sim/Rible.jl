
"""
Linear parameter function: `u = W φ + b`.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct LinearParamFun{T} <: AbstractParamFun
    weights::Matrix{T}
    bias::Vector{T}
end

function evaluate_paramfun(p::LinearParamFun, ϕ)
    return p.weights * ϕ .+ p.bias
end


"""
    LinearFeedbackPolicy(K, k=zeros(...); extractor, featurizer)

Convenience constructor returning a `FeedbackPolicy` configured with
`LinearParamFun` (`u = K * x + k`). No dedicated wrapper struct is used; dispatch
relies on the field types of `FeedbackPolicy`.
"""
function LinearFeedbackPolicy(
        K::AbstractMatrix{T},
        k::AbstractVector{T}=zeros(T, size(K, 1));
        extractor::E,
        featurizer::F=IdentityFeaturizer(),
    ) where {T, E<:AbstractStateExtractor, F<:AbstractFeaturizer}
    paramfun = LinearParamFun{T}(K, k)
    return FeedbackPolicy(extractor, featurizer, paramfun)
end

# Alias for dispatch on linear-feedback-configured policies
const LFPolicy = FeedbackPolicy{<:AbstractStateExtractor,<:AbstractFeaturizer,<:LinearParamFun}

# --- Parameter helpers ---

function get_params(p::FeedbackPolicy{<:AbstractStateExtractor,<:AbstractFeaturizer, <:LinearParamFun})
    pf = p.paramfun
    return vcat(vec(pf.weights), pf.bias)
end

function set_params!(p::LFPolicy, params)
    pf = p.paramfun
    nK = length(pf.weights)
    pf.weights .= reshape(params[1:nK], size(pf.weights))
    pf.bias .= params[nK+1:end]
end

function update!(p::LFPolicy) end

function get_num_of_params(p::LFPolicy)
    pf = p.paramfun
    return length(pf.weights) + length(pf.bias)
end

# --- Differentials ---
# Legacy actuator bridge - uses the policy's extractor to get state from gauges
function actuate!(bot::Robot, policy::LFPolicy, inst_state::AbstractCoordinatesState)
    (;structure,hub) = bot
    # Extract state using the policy's extractor (measures from hub.gauges)
    update_bodies!(structure, inst_state)
    update_apparatuses!(structure, inst_state.s)
    x = extract(policy.extractor, bot)
    ϕ = features(policy.featurizer, x)
    u = evaluate_paramfun(policy.paramfun, ϕ)
    hub.state.u .= u
    (;actuators, coalition) = hub
    (;actid2sys_actions_idx) = coalition

    foreach(actuators) do actuator
        idx = actid2sys_actions_idx[actuator.id]
        execute!(structure, actuator, (@view hub.state.u[idx]))
    end
end

"""
    gen_force_state_jacobian!(∂F∂q̌, ∂F∂q̌̇, ∂F∂u, policy::LFPolicy, bot::Robot, inst_state)

Compute generalized force Jacobian F(q,q̇,u(q,q̇)) for a given policy.
This function updates the provided matrices ∂F∂q̌ and ∂F∂q̌̇ in place.
Formula:
∂F∂q̌ = ∂F∂q̌ + ∂F∂u * ∂u∂q̌
∂F∂q̌̇ = ∂F∂q̌̇ + ∂F∂u * ∂u∂q̌̇

This function concerns only the latter part of the formula: 
∂F∂u * ∂u∂q̌ and ∂F∂u * ∂u∂q̌̇
"""
function gen_force_state_jacobian!(∂F∂q̌, ∂F∂q̌̇, ∂F∂s, ∂F∂u, policy::LFPolicy, bot::Robot, inst_state::AbstractCoordinatesState)
    (;structure, hub) = bot
    (;t, q, q̇) = inst_state
    (;actuators, coalition) = hub
    T = get_numbertype(structure)
    nq̌ = get_num_of_free_coords(structure)
    nq = get_num_of_full_coords(structure)
    (;num_of_actions, actid2sys_actions_idx) = coalition

    update_bodies!(structure, inst_state)
    update_apparatuses!(structure, inst_state.s)

    ∂F∂u .= 0
    foreach(actuators) do actuator
        u_idx = actid2sys_actions_idx[actuator.id]
        gen_force_actu_jacobian!(
            (@view ∂F∂u[:, u_idx]),
            structure,
            actuator,
            (@view hub.state.u[u_idx])
        )
    end

    # Get extractor Jacobian: ∂x/∂q, ∂x/∂q̇ (full coords)
    Jq, Jq̇, Js = extract_jacobian(policy.extractor, bot)
    K = policy.paramfun.weights
    
    # Chain rule: ∂u/∂q = K · ∂x/∂q
    ∂u∂q = K * Jq   # (nu × nq)
    ∂u∂q̇ = K * Jq̇   # (nu × nq)
    ∂u∂s = K * Js   # (nu × ns)
    
    # Extract free coordinates portion (first nq̌ columns)
    ∂u∂q̌ = @view ∂u∂q[:, 1:nq̌]
    ∂u∂q̌̇ = @view ∂u∂q̇[:, 1:nq̌]

    mul!(∂F∂q̌, ∂F∂u, ∂u∂q̌)
    mul!(∂F∂q̌̇, ∂F∂u, ∂u∂q̌̇)
    mul!(∂F∂s, ∂F∂u, ∂u∂s)
end


"""
    vjp_wrt_state(v, policy::LFPolicy, bot::Robot, num_of_actions, solver, solver_state)

Vector-Jacobian product for reverse-mode differentiation using capta gauges.
"""
function vjp_wrt_state(v, policy::LFPolicy, bot::Robot, num_of_actions, solver, solver_state)
    (;structure, hub) = bot
    (;
        qₖ, pₖ, qₖ₊₁, pₖ₊₁, λₘ, tₖ, tₖ₊₁, tₘ, qₘ, q̇ₘ, dt,
    ) = solver_state
    nq = get_num_of_full_coords(structure)
    nλ = get_num_of_cstr(structure)
    ns = get_num_of_aux_var(structure)
    T = typeof(tₖ)
    K = policy.paramfun.weights
    qₘ .= (qₖ .+ qₖ₊₁)./2
    q̇ₘ .= (qₖ₊₁ .- qₖ)./dt
    update_bodies!(structure, solver_state.state_mid)
    update_apparatuses!(structure, solver_state.state_mid.s)
    # Get extractor Jacobian via gauge system: ∂x/∂q, ∂x/∂q̇
    Jq, Jq̇, Js = extract_jacobian(policy.extractor, bot)
    
    # Chain rule: ∂u/∂q = K · ∂x/∂q, ∂u/∂q̇ = K · ∂x/∂q̇
    ∂u∂q = K * Jq   # (nu × nq)
    ∂u∂q̇ = K * Jq̇   # (nu × nq)
    ∂u∂s = K * Js   # (nu × ns)
    
    # Midpoint transformation VJP
    ∂u∂qₖ   = ∂u∂q ./ 2 .- ∂u∂q̇ ./ dt
    ∂u∂qₖ₊₁ = ∂u∂q ./ 2 .+ ∂u∂q̇ ./ dt
    ∂u∂sₖ   = ∂u∂s ./ 2
    ∂u∂sₖ₊₁ = ∂u∂s ./ 2
    
    ∂ϕ∂qₖᵀ   = (v' * ∂u∂qₖ)'
    ∂ϕ∂qₖ₊₁ᵀ = (v' * ∂u∂qₖ₊₁)'
    ∂ϕ∂sₖᵀ = (v' * ∂u∂sₖ)'
    ∂ϕ∂sₖ₊₁ᵀ = (v' * ∂u∂sₖ₊₁)'

    ∂ϕ∂pₖᵀ = spzeros(T, nq)
    ∂ϕ∂pₖ₊₁ᵀ = spzeros(T, nq)
    ∂ϕ∂λᵀ = spzeros(T, nλ)

    ∂ϕ∂qₖᵀ, ∂ϕ∂qₖ₊₁ᵀ, ∂ϕ∂pₖᵀ, ∂ϕ∂pₖ₊₁ᵀ, ∂ϕ∂λᵀ, ∂ϕ∂sₖᵀ, ∂ϕ∂sₖ₊₁ᵀ
end

function accumulate_param_grad!(grad_storage, policy::LFPolicy, v_total, solver_state, bot)
    (;structure) = bot
    (; qₖ, qₖ₊₁, qₘ, q̇ₘ, dt ) = solver_state
    qₘ .= (qₖ .+ qₖ₊₁)./2
    q̇ₘ .= (qₖ₊₁ .- qₖ)./dt
    update_bodies!(structure, solver_state.state_mid)
    update_apparatuses!(structure, solver_state.state_mid.s)
    
    x = extract(policy.extractor, bot)
    nx = length(x)
    num_of_actions = size(policy.paramfun.weights, 1)
    
    for i in 1:num_of_actions
        for j in 1:nx
            k = (j-1) * num_of_actions + i
            grad_storage[k] += v_total[i] * x[j]
        end
    end
    
    offset = num_of_actions * nx
    for i in 1:num_of_actions
        grad_storage[offset + i] += v_total[i]
    end
end
