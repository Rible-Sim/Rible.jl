"""
    DiscreteLinearFeedbackPolicy(K, k=zeros(...); extractor, featurizer)

Convenience constructor returning a `DiscreteFeedbackPolicy` configured with
`LinearParamFun` (`u = K * x + k`). No dedicated wrapper struct is used; dispatch
relies on the field types of `DiscreteFeedbackPolicy`.
"""
function DiscreteLinearFeedbackPolicy(
        K::AbstractMatrix{T},
        k::AbstractVector{T}=zeros(T, size(K, 1));
        extractor::E,
        featurizer::F=IdentityFeaturizer(),
    ) where {T, E<:AbstractStateExtractor, F<:AbstractFeaturizer}
    paramfun = LinearParamFun{T}(K, k)
    return DiscreteFeedbackPolicy(extractor, featurizer, paramfun)
end

"""
Alias for dispatch on linear-feedback-configured policies.
$(TYPEDEF)
"""
const DiscreteLFPolicy = DiscreteFeedbackPolicy{<:AbstractStateExtractor,<:AbstractFeaturizer,<:LinearParamFun}

# --- Parameter helpers ---

function get_params(p::DiscreteLFPolicy)
    pf = p.paramfun
    return vcat(vec(pf.weights), pf.bias)
end

function set_params!(p::DiscreteLFPolicy, params)
    pf = p.paramfun
    nK = length(pf.weights)
    pf.weights .= reshape(params[1:nK], size(pf.weights))
    pf.bias .= params[nK+1:end]
end

function update!(p::DiscreteLFPolicy) end

function get_num_of_params(p::DiscreteLFPolicy)
    pf = p.paramfun
    return length(pf.weights) + length(pf.bias)
end

# --- Differentials ---

"""
Execute the control actions on the robot's apparatuses.
$(TYPEDSIGNATURES)
"""
function actuate!(bot::Robot, policy::DiscreteLFPolicy, inst_state::AbstractCoordinatesState)
    (;structure,hub) = bot
    (;actuators, coalition) = hub
    (;actid2sys_actions_idx) = coalition

    foreach(actuators) do actuator
        idx = actid2sys_actions_idx[actuator.id]
        execute!(structure, actuator, (@view hub.state.u[idx]))
    end
end

"""
Compute the control actions and store them in the robot's hub state.
$(TYPEDSIGNATURES)
"""
function control!(bot::Robot,policy::DiscreteLFPolicy,solver_cache,solver_state)
    (;structure,hub) = bot
    inst_state = solver_state.state_k
    update_bodies!(structure, inst_state)
    update_apparatuses!(structure, inst_state.s)
    x = extract(policy.extractor, bot)
    ϕ = features(policy.featurizer, x)
    u = evaluate_paramfun(policy.paramfun, ϕ)

    hub.state.u .= u
end

"""
    control_jacobian!(bot, policy::DiscreteLFPolicy, ...)

For discrete time policy u = u(qₖ).
"""
function control_jacobian!(jacobian_workspace, bot::Robot,policy::DiscreteLFPolicy,
        solver_cache,solver_state
    )
    (;M⁻¹ₘ, ∂Fₘ∂u, ∂C∂qₖ, ∂C∂pₖ) = jacobian_workspace
    control!(bot,policy,solver_cache,solver_state)
    (;structure,hub) = bot
    (;actuators,coalition) = hub
    (;actid2sys_actions_idx) = coalition
   
    Jq, Jq̇, Js = extract_jacobian(policy.extractor, bot)
    K = policy.paramfun.weights
    
    ∂u∂q = K * Jq
    ∂u∂p = K * Jq̇ * M⁻¹ₘ
    
    lazy_update_bodies!(structure,solver_state.state_mid)
    # Must apply the control `u` to the structure (e.g. set rest lengths) BEFORE updating apparatuses
    # because update_apparatuses! calculates derived state (like tension) that depends on both geometry and control parameters.
    actuate!(bot, policy, solver_state.state_mid)
    update_apparatuses!(structure, solver_state.state_mid.s)
    ∂F∂u = ∂Fₘ∂u .= 0
    foreach(actuators) do actuator
        u_idx = actid2sys_actions_idx[actuator.id]
        gen_force_actu_jacobian!(
            (@view ∂F∂u[:, u_idx]),
            structure,
            actuator,
            (@view hub.state.u[u_idx])
        )
    end
    ∂C∂qₖ .= ∂F∂u * ∂u∂q
    ∂C∂pₖ .= ∂F∂u * ∂u∂p
end


"""
    gen_force_state_jacobian!(∂F∂q, ∂F∂q̇, ∂F∂u, policy::DiscreteLFPolicy, bot::Robot, inst_state)

Zero out gradients because control is fixed w.r.t integration variables q_mid or q_{k+1}.
"""
function gen_force_state_jacobian!(∂F∂q, ∂F∂q̇, ∂F∂s, ∂F∂u, policy::DiscreteLFPolicy, bot::Robot, inst_state::AbstractCoordinatesState)
    (;structure, hub) = bot
    (;t, q, q̇) = inst_state
    (;actuators, coalition) = hub
    T = get_numbertype(structure)
    nq = get_num_of_free_coords(structure)
    nq = get_num_of_full_coords(structure)
    (;num_of_actions, actid2sys_actions_idx) = coalition

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
    # Zero contribution. 
    # Usually this function adds to ∂F/∂q and ∂F/∂q̇ via ∂u/∂q.
    # Since ∂u/∂q (w.r.t mid) is 0, we do not modify gradients.
end


"""
    vjp_wrt_state(v, policy::DiscreteLFPolicy, bot::Robot, num_of_actions, solver_cache, solver_state)

Vector-Jacobian product for reverse-mode differentiation using discrete `q_k`.
"""
function vjp_wrt_state(v, policy::DiscreteLFPolicy, bot::Robot, num_of_actions, solver_cache, solver_state)
    (;structure, hub) = bot
    (;
        qₖ, pₖ, qₖ₊₁, pₖ₊₁, λₘ, tₖ, tₖ₊₁, tₘ, qₘ, q̇ₘ, dt,
    ) = solver_state
    (; M⁻¹ₘ) = solver_cache.jacobian_workspace
    nq = get_num_of_full_coords(structure)
    nλ = get_num_of_cstr(structure)
    ns = get_num_of_aux_var(structure)
    T = typeof(tₖ)
    K = policy.paramfun.weights

    inst_state = solver_state.state_k
    update_bodies!(structure, inst_state)
    update_apparatuses!(structure, inst_state.s)

    # Get extractor Jacobian via gauge system: ∂x/∂q, ∂x/∂p
    Jq, Jq̇, Js = extract_jacobian(policy.extractor, bot)
    
    # Chain rule: ∂u/∂q = K · ∂x/∂q, ∂u/∂p = K · ∂x/∂p
    ∂u∂q = K * Jq          # (nu × nq)
    ∂u∂p = K * Jq̇ * M⁻¹ₘ   # (nu × nq)
    
    ∂ϕ∂qₖᵀ = (v' * ∂u∂q)'
    ∂ϕ∂pₖᵀ = (v' * ∂u∂p)'
    
    ∂ϕ∂qₖ₊₁ᵀ = spzeros(T, nq)
    ∂ϕ∂pₖ₊₁ᵀ = spzeros(T, nq)
    ∂ϕ∂λᵀ = spzeros(T, nλ)
    ∂ϕ∂sₖᵀ = spzeros(T, ns)
    ∂ϕ∂sₖ₊₁ᵀ = spzeros(T, ns)

    ∂ϕ∂qₖᵀ, ∂ϕ∂qₖ₊₁ᵀ, ∂ϕ∂pₖᵀ, ∂ϕ∂pₖ₊₁ᵀ, ∂ϕ∂λᵀ, ∂ϕ∂sₖᵀ, ∂ϕ∂sₖ₊₁ᵀ
end

function accumulate_param_grad!(grad_storage, policy::DiscreteLFPolicy, v_total, solver_state, bot)
    (;structure) = bot
    
    inst_state = solver_state.state_k
    update_bodies!(structure, inst_state)
    update_apparatuses!(structure, inst_state.s)
    
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

