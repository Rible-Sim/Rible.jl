
"""
    DiscreteProportionalDerivativePolicy(K, x_des; extractor, featurizer)

Convenience constructor that returns a `DiscreteFeedbackPolicy` configured as a PD
controller using `PDParamFun`. The control law is u = K(x_des - x).

This is equivalent to the linear form u = -K*x + b where b = K*x_des,
but with structured parameterization [K, x_des] instead of [K, b].
"""
function DiscreteProportionalDerivativePolicy(
        K::AbstractArray{T},
        x_des=zeros(T, size(K, 2));
        extractor::E,
        featurizer::F=IdentityFeaturizer(),
    ) where {T, E<:AbstractStateExtractor, F<:AbstractFeaturizer}

    @assert size(K, 2) == length(x_des) "x_des must match K columns"
    paramfun = PDParamFun{T}(K, x_des)
    return DiscreteFeedbackPolicy(extractor, featurizer, paramfun)
end

"""
Local alias for PD-configured feedback policies.
$(TYPEDEF)
"""
const DiscretePDPolicy = DiscreteFeedbackPolicy{<:AbstractStateExtractor,<:AbstractFeaturizer,<:PDParamFun}

# --- Parameter helpers ---

function get_params(p::DiscretePDPolicy)
    pf = p.paramfun
    return vcat(vec(pf.K), pf.x)
end

function set_params!(p::DiscretePDPolicy, params)
    pf = p.paramfun
    nK = length(pf.K)
    nx = length(pf.x)

    pf.K .= reshape(params[1:nK], size(pf.K))
    pf.x .= params[nK+1:nK+nx]
end

function get_num_of_params(p::DiscretePDPolicy)
    pf = p.paramfun
    return length(pf.K) + length(pf.x)
end

function set_reference!(p::DiscretePDPolicy, x_des)
    pf = p.paramfun
    @assert length(x_des) == length(pf.x) "x_des size mismatch"
    pf.x .= x_des
end


function update!(p::DiscretePDPolicy) end


# Legacy actuator bridge - now simplified to just apply the control inputs
"""
Execute the control actions on the robot's apparatuses.
$(TYPEDSIGNATURES)
"""
function actuate!(bot::Robot, policy::DiscretePDPolicy, inst_state::AbstractCoordinatesState)
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
function control!(bot::Robot,policy::DiscretePDPolicy,solver_cache,solver_state)
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
    control_jacobian!(bot, policy::DiscretePDPolicy, ...)

For discrete time policy u = u(qₖ).
We compute ∂C∂qₖ and ∂C∂pₖ.
∂u/∂qₖ₊₁ = 0.
"""
function control_jacobian!(jacobian_workspace, bot::Robot,policy::DiscretePDPolicy,
        solver_cache,
        solver_state
    )
    control!(bot,policy,solver_cache,solver_state)
    (;M⁻¹ₘ, ∂Fₘ∂u, ∂C∂qₖ, ∂C∂pₖ, ∂C∂sₖ) = jacobian_workspace
    (;structure,hub) = bot
    (;actuators,coalition) = hub
    (;actid2sys_actions_idx) = coalition
    
    Jq, Jq̇, Js = extract_jacobian(policy.extractor, bot)
    K = policy.paramfun.K
    
    # ∂u/∂q = -K * Jq
    # ∂C∂qₖ = ∂F∂u * (-K * Jq) = - (∂F∂u * K) * Jq
    
    # We can do this efficiently or just direct mul
    # Avoiding allocs if possible, but compact code first.
    ∂u∂q = -K * Jq
    ∂u∂q̇ = -K * Jq̇
    ∂u∂p = ∂u∂q̇ * M⁻¹ₘ
    ∂u∂s = -K * Js
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
    ∂C∂sₖ .= ∂F∂u * ∂u∂s
end


"""
    gen_force_state_jacobian!(∂F∂q, ∂F∂q̇, ∂F∂u, policy::DiscretePDPolicy, bot::Robot, inst_state)

Zero out gradients because control is fixed w.r.t integration variables q_mid or q_{k+1}.
"""
function gen_force_state_jacobian!(∂F∂q, ∂F∂q̇, ∂F∂s, ∂F∂u, policy::DiscretePDPolicy, bot::Robot, inst_state::AbstractCoordinatesState)
    (;structure, hub) = bot
    (;t, q, q̇) = inst_state
    (;actuators, coalition) = hub
    T = get_numbertype(structure)
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
    # Control is fixed w.r.t midpoint states; Jacobian contribution corresponds to ∂u/∂q_mid which is 0.
end


# --- Differentials ---

"""
    vjp_wrt_state(v, policy::DiscretePDPolicy, bot::Robot, num_of_actions, solver_cache, solver_state)

Vector-Jacobian product for reverse-mode differentiation.
u = K(x_des - x(qₖ))
"""
function vjp_wrt_state(v, policy::DiscretePDPolicy, bot::Robot, num_of_actions, solver_cache, solver_state)
    (;structure, hub) = bot
    (;
        qₖ, pₖ, qₖ₊₁, pₖ₊₁, λₘ, tₖ, tₖ₊₁, tₘ, qₘ, q̇ₘ, dt,
    ) = solver_state
    (; M⁻¹ₘ) = solver_cache.jacobian_workspace
    nq = get_num_of_full_coords(structure)
    nλ = get_num_of_cstr(structure)
    ns = get_num_of_aux_var(structure)
    T = typeof(tₖ)
    K = policy.paramfun.K
    
    inst_state = solver_state.state_k
    update_bodies!(structure, inst_state)
    update_apparatuses!(structure, inst_state.s)

    # Get extractor Jacobian: ∂x/∂qₖ, ∂x/∂pₖ
    Jq, Jq̇, Js = extract_jacobian(policy.extractor, bot)
    
    # Chain rule: ∂u/∂q = -K · ∂x/∂q
    # (negative because u = K(x_des - x))
    ∂u∂q = -K * Jq   # (nu × nq)
    ∂u∂p = -K * Jq̇ * M⁻¹ₘ   # (nu × nq)
    ∂u∂s = -K * Js

    ∂ϕ∂qₖᵀ = (v' * ∂u∂q)'
    ∂ϕ∂pₖᵀ = (v' * ∂u∂p)'
    ∂ϕ∂sₖᵀ = (v' * ∂u∂s)'
    
    ∂ϕ∂qₖ₊₁ᵀ = spzeros(T, nq)
    ∂ϕ∂pₖ₊₁ᵀ = spzeros(T, nq)
    ∂ϕ∂λᵀ = spzeros(T, nλ)
    ∂ϕ∂sₖ₊₁ᵀ = spzeros(T, ns)

    ∂ϕ∂qₖᵀ, ∂ϕ∂qₖ₊₁ᵀ, ∂ϕ∂pₖᵀ, ∂ϕ∂pₖ₊₁ᵀ, ∂ϕ∂λᵀ, ∂ϕ∂sₖᵀ, ∂ϕ∂sₖ₊₁ᵀ
end

function accumulate_param_grad!(grad_storage, policy::DiscretePDPolicy, v_total, solver_state, bot)
    (;structure) = bot
    num_of_actions = size(policy.paramfun.K, 1)
    K = policy.paramfun.K
    x_des = policy.paramfun.x
    
    inst_state = solver_state.state_k
    update_bodies!(structure, inst_state)
    update_apparatuses!(structure, inst_state.s)

    x = extract(policy.extractor, bot)
    error = x_des .- x

    nx = length(x)
    nK = length(K)

    # v_total is size (num_of_actions,)
    # ∇_K = v_total * error'
    for i in 1:num_of_actions
        for j in 1:nx
            k = (j-1) * num_of_actions + i
            grad_storage[k] += v_total[i] * error[j]
        end
    end
    
    # ∇_x_des = K' * v_total
    # Since K is (nu, nx), K' is (nx, nu). K' * v_total is (nx,)
    grad_storage[nK+1:nK+nx] .+= K' * v_total
end
