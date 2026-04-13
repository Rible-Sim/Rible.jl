
"""
    PDParamFun{T}

PD control parametric function: u = K(x_des - x)

This is equivalent to the linear form u = -K*x + b where b = K*x_des,
but with structured parameterization that separates the gain matrix K
and the reference/setpoint x_des.

# Fields
- `K::Matrix{T}`: Gain matrix (action_dim × state_dim)
- `x::Vector{T}`: Reference/desired state (state_dim,)

# Control Law
u = K * (x_des - x_measured)

# Parameters
The parameters are [vec(K); x], giving a bilinear parameterization where
the bias term b = K*x_des couples the gain and reference.

# Notes
- This is a structured but over-parameterized version of LinearParamFun
- Parameters are nonlinearly related: changing K affects how x_des influences u
- Useful when you want to separate gain tuning from setpoint adjustment
"""
mutable struct PDParamFun{T} <: AbstractParamFun
    K::Matrix{T}
    x::Vector{T}
end

function evaluate_paramfun(p::PDParamFun, x_measured)
    return p.K * (p.x .- x_measured)
end

"""
    ProportionalDerivativePolicy(K, x_des; extractor, featurizer)

Convenience constructor that returns a `FeedbackPolicy` configured as a PD
controller using `PDParamFun`. The control law is u = K(x_des - x).

This is equivalent to the linear form u = -K*x + b where b = K*x_des,
but with structured parameterization [K, x_des] instead of [K, b].
"""
function ProportionalDerivativePolicy(
        K::AbstractArray{T},
        x_des=zeros(T, size(K, 2));
        extractor::E,
        featurizer::F=IdentityFeaturizer(),
    ) where {T, E<:AbstractStateExtractor, F<:AbstractFeaturizer}

    @assert size(K, 2) == length(x_des) "x_des must match K columns"
    paramfun = PDParamFun{T}(K, x_des)
    return FeedbackPolicy(extractor, featurizer, paramfun)
end

# Local alias for PD-configured feedback policies
const PDFeedbackPolicy = FeedbackPolicy{<:AbstractStateExtractor,<:AbstractFeaturizer,<:PDParamFun}


# --- Parameter helpers ---

function get_params(p::PDFeedbackPolicy)
    pf = p.paramfun
    return vcat(vec(pf.K), pf.x)
end

function set_params!(p::PDFeedbackPolicy, params)
    pf = p.paramfun
    nK = length(pf.K)
    nx = length(pf.x)

    pf.K .= reshape(params[1:nK], size(pf.K))
    pf.x .= params[nK+1:nK+nx]
end

function get_num_of_params(p::PDFeedbackPolicy)
    pf = p.paramfun
    return length(pf.K) + length(pf.x)
end

function set_reference!(p::PDFeedbackPolicy, x_des)
    pf = p.paramfun
    @assert length(x_des) == length(pf.x) "x_des size mismatch"
    pf.x .= x_des
end


function update!(p::PDFeedbackPolicy) end

# Legacy actuator bridge - uses the policy's extractor to get state from gauges
function actuate!(bot::Robot, policy::PDFeedbackPolicy, inst_state::AbstractCoordinatesState)
    (;structure,hub) = bot
    # Extract state using the policy's extractor (measures from hub.gauges)
    update_bodies!(structure, inst_state)
    update_apparatuses!(structure, inst_state.s)
    x = extract(policy.extractor, bot)

    # u = K*x_des - K*x
    # PD control law: u = K(x_des - x) where x can be [q; v] or any state vector
    # The extractor determines what x is, and K/x_des must match that form
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
    gen_force_state_jacobian!(∂F∂q, ∂F∂q̇, ∂F∂u, policy::PDFeedbackPolicy, bot::Robot, inst_state)

Compute generalized force Jacobian F(q,q̇,u(q,q̇)) for PD policy.

For u = K(x_des - x) = -K*x + K*x_des:
∂F/∂q̌ += ∂F/∂u · (-K · ∂x/∂q̌)
∂F/∂q̌̇ += ∂F/∂u · (-K · ∂x/∂q̌̇)
"""
function gen_force_state_jacobian!(∂F∂q, ∂F∂q̇, ∂F∂s, ∂F∂u, policy::PDFeedbackPolicy, bot::Robot, inst_state::AbstractCoordinatesState)
    (;structure, hub) = bot
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
    K = policy.paramfun.K
    
    # Chain rule: ∂u/∂q = -K · ∂x/∂q (note the negative sign!)
    ∂u∂q = -K * Jq   # (nu × nq)
    ∂u∂q̇ = -K * Jq̇   # (nu × nq)
    ∂u∂s = -K * Js   # (nu × nq)
    
    # Extract free coordinates portion (first nq̌ columns)
    ∂u∂q̌ = @view ∂u∂q[:, 1:nq̌]
    ∂u∂q̌̇ = @view ∂u∂q̇[:, 1:nq̌]

    mul!(∂F∂q, ∂F∂u, ∂u∂q̌, one(T), one(T))
    mul!(∂F∂q̇, ∂F∂u, ∂u∂q̌̇, one(T), one(T))
    mul!(∂F∂s, ∂F∂u, ∂u∂s, one(T), one(T))
end


# --- Differentials ---

"""
    vjp_wrt_state(v, policy::PDFeedbackPolicy, bot::Robot, num_of_actions, solver, solver_state)

Vector-Jacobian product for reverse-mode differentiation using capta gauges.

For PD control u = K(x_des - x), the VJP computes v'·∂u/∂z for all parameters z.
"""
function vjp_wrt_state(v, policy::PDFeedbackPolicy, bot::Robot, num_of_actions, solver, solver_state)
    (;structure, hub) = bot
    (;
        qₖ, pₖ, qₖ₊₁, pₖ₊₁, λₘ, tₖ, tₖ₊₁, tₘ, qₘ, q̇ₘ, dt,
    ) = solver_state
    nq = get_num_of_full_coords(structure)
    nλ = get_num_of_cstr(structure)
    ns = get_num_of_aux_var(structure)
    T = typeof(tₖ)
    K = policy.paramfun.K
    
    update_bodies!(structure, solver_state.state_mid)
    update_apparatuses!(structure, solver_state.state_mid.s)
    
    x = extract(policy.extractor, bot)
    # Get extractor Jacobian: ∂x/∂q, ∂x/∂q̇
    Jq, Jq̇, Js = extract_jacobian(policy.extractor, bot)
    
    # Chain rule: ∂u/∂q = -K · ∂x/∂q, ∂u/∂q̇ = -K · ∂x/∂q̇
    # (negative because u = K(x_des - x) = -K*x + K*x_des)
    ∂u∂q = -K * Jq   # (nu × nq)
    ∂u∂q̇ = -K * Jq̇   # (nu × nq)
    ∂u∂s = -K * Js   # (nu × nq)
    
    # Midpoint transformation VJP
    ∂u∂qₖ   = ∂u∂q ./ 2 .- ∂u∂q̇ ./ dt
    ∂u∂qₖ₊₁ = ∂u∂q ./ 2 .+ ∂u∂q̇ ./ dt
    
    ∂u∂sₖ   = ∂u∂s ./ 2 
    ∂u∂sₖ₊₁ = ∂u∂s ./ 2 
    
    ∂ϕ∂qₖᵀ   = (v' * ∂u∂qₖ)'
    ∂ϕ∂qₖ₊₁ᵀ = (v' * ∂u∂qₖ₊₁)'

    ∂ϕ∂pₖᵀ = spzeros(T, nq)
    ∂ϕ∂pₖ₊₁ᵀ = spzeros(T, nq)
    ∂ϕ∂λᵀ = spzeros(T, nλ)
    ∂ϕ∂sₖᵀ   = (v' * ∂u∂sₖ)'
    ∂ϕ∂sₖ₊₁ᵀ = (v' * ∂u∂sₖ₊₁)'

    ∂ϕ∂qₖᵀ, ∂ϕ∂qₖ₊₁ᵀ, ∂ϕ∂pₖᵀ, ∂ϕ∂pₖ₊₁ᵀ, ∂ϕ∂λᵀ, ∂ϕ∂sₖᵀ, ∂ϕ∂sₖ₊₁ᵀ
end

function accumulate_param_grad!(grad_storage, policy::PDFeedbackPolicy, v_total, solver_state, bot)
    (;structure) = bot
    
    update_bodies!(structure, solver_state.state_mid)
    update_apparatuses!(structure, solver_state.state_mid.s)
    
    x = extract(policy.extractor, bot)
    
    K = policy.paramfun.K
    x_des = policy.paramfun.x
    
    error = x_des .- x
    nx = length(x)
    nK = length(K)
    num_of_actions = size(K, 1)
    
    for i in 1:num_of_actions
        for j in 1:nx
            k = (j-1) * num_of_actions + i
            grad_storage[k] += v_total[i] * error[j]
        end
    end
    
    grad_storage[nK+1:nK+nx] .+= K' * v_total
end
