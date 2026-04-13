

struct StaticLinearPolicy{mat_biasT,matT,biasT} <: AbstractPolicy
    mat_bias::mat_biasT
    mat::matT
    bias::biasT
end

function set_params!(policy::StaticLinearPolicy,params)
    policy.mat_bias[:] .= params
end

function get_params(policy::StaticLinearPolicy)
    vec(policy.mat_bias)
end

function get_num_of_params(policy::StaticLinearPolicy)
    prod(size(policy.mat_bias))
end


function update!(policy::StaticLinearPolicy) end


function actuate!(bot::Robot,policy::StaticLinearPolicy,inst_state::AbstractCoordinatesState)
    (;structure,hub) = bot
    (;q,q̇) = inst_state
    (;actuators,state,coalition,) = hub
    (;actid2sys_actions_idx) = coalition
    # Calculate control actions u = Kx + b, 
    # x = [q; q̇] is the state vector, 
    # K is the policy matrix, 
    # b is the bias vector.
    state.u .= policy.mat * vcat(q,q̇) .+ policy.bias
    foreach(actuators) do actuator
        idx = actid2sys_actions_idx[actuator.id]
        execute!(
            structure,
            actuator,
            (@view state.u[idx])
        )
    end
end

function gen_force_state_jacobian!(∂F∂q,∂F∂q̇,∂F∂s,∂F∂u,policy::StaticLinearPolicy,bot::Robot,inst_state::AbstractCoordinatesState;)
    (;structure, hub) = bot
    (;q,q̇,t) = inst_state
    (;actuators,coalition) = hub
    T = get_numbertype(structure)
    nq = get_num_of_free_coords(structure)
    (;num_of_actions,actid2sys_actions_idx) = coalition
    ∂F∂u .= 0
    foreach(actuators) do actuator
        u_idx = actid2sys_actions_idx[actuator.id]
        gen_force_actu_jacobian!(
            (@view ∂F∂u[:,u_idx]),
            structure,
            actuator,
            (@view hub.state.u[u_idx])
        )
    end
    # u = Kx + b
    # x = [q; q̇] is the state vector, 
    # K is the policy matrix, 
    # b is the bias vector.
    # ∂u/∂x = K
    # ∂u/∂q̌ = K[:,1:nq]
    # ∂u/∂q̌̇ = K[:,nq+1:2nq]
    # ∂u/∂b = I
    ∂u∂x = policy.mat
    ∂u∂q̌ = @view ∂u∂x[:,   1: nq]
    ∂u∂q̌̇ = @view ∂u∂x[:,nq+1:2nq]
    ∂F∂q .= ∂F∂u*∂u∂q̌
    ∂F∂q̇ .= ∂F∂u*∂u∂q̌̇
end



function vjp_wrt_state(v,policy::StaticLinearPolicy,bot::Robot,num_of_actions,solver,solver_state)
    (;structure) = bot
    (;
        qₖ, pₖ, qₖ₊₁, pₖ₊₁, λₘ, tₖ, tₖ₊₁, qₘ, q̇ₘ, dt
    ) = solver_state
    nq = get_num_of_full_coords(structure)
    nλ = get_num_of_cstr(structure)
    ns = get_num_of_aux_var(structure)

    ∂ϕ∂uᵀ = v
    
    ∂ϕ∂qₖᵀ = (v'*(policy.mat[:,1:nq] ./2 .- policy.mat[:,nq+1:2nq] ./dt))'
    ∂ϕ∂qₖ₊₁ᵀ = (v'*(policy.mat[:,1:nq] ./2 .+ policy.mat[:,nq+1:2nq] ./dt))'

    ∂ϕ∂pₖᵀ = spzeros(typeof(tₖ),nq)
    ∂ϕ∂pₖ₊₁ᵀ = spzeros(typeof(tₖ),nq)
    ∂ϕ∂λᵀ = spzeros(typeof(tₖ),nλ)
    ∂ϕ∂sₖᵀ = spzeros(typeof(tₖ),ns)
    ∂ϕ∂sₖ₊₁ᵀ = spzeros(typeof(tₖ),ns)

    ∂ϕ∂qₖᵀ,  ∂ϕ∂qₖ₊₁ᵀ, ∂ϕ∂pₖᵀ, ∂ϕ∂pₖ₊₁ᵀ, ∂ϕ∂λᵀ, ∂ϕ∂sₖᵀ, ∂ϕ∂sₖ₊₁ᵀ
end

function accumulate_param_grad!(grad_storage, policy::StaticLinearPolicy, v_total, solver_state, bot)
    (; qₘ, q̇ₘ ) = solver_state
    x = vcat(qₘ, q̇ₘ)
    # v_total is size (num_of_actions,)
    # policy.mat_bias shape is (num_of_actions, nx + 1)
    # ∇_mat = v_total * x', ∇_bias = v_total
    
    num_of_actions = size(policy.mat_bias, 1)
    nx = length(x)
    
    # Store gradients in a column-major vector as expected by `get_params` / `set_params!`
    for i in 1:num_of_actions
        for j in 1:nx
            k = (j-1) * num_of_actions + i
            grad_storage[k] += v_total[i] * x[j]
        end
    end
    
    # Bias is the last block
    offset = num_of_actions * nx
    for i in 1:num_of_actions
        grad_storage[offset + i] += v_total[i]
    end
end

