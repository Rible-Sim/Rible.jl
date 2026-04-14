

struct TrajectoryPolicy{T,time_nodesT,actuations_trajT} <: AbstractPolicy
    num_of_actions::T
    time_nodes::time_nodesT
    actuations_traj::actuations_trajT
end

function constant_interpolation(time_nodes::AbstractVector,actuation_traj::AbstractVector)
    (t) -> begin
        idx = searchsortedlast(time_nodes, t) 
        if idx == 0
            val = actuation_traj[begin] # Before first node, use first value
        elseif idx == length(time_nodes)  
            val = actuation_traj[end] # After last node, use last value
        else
            val = actuation_traj[idx] # Use value at left endpoint of interval
        end
        val
    end
end

function TrajectoryPolicy(time_nodes::AbstractVector,actuations_traj::AbstractMatrix)
    num_of_actions = size(actuations_traj,2)
    # itp = [
    #     constant_interpolation(time_nodes,actuation_traj)
    #     for actuation_traj in eachcol(actuations_traj)
    # ]
    TrajectoryPolicy(num_of_actions,time_nodes,actuations_traj)
end

function get_num_of_params(policy::TrajectoryPolicy)
    (;actuations_traj) = policy
    prod(size(actuations_traj))
end

function set_params!(policy::TrajectoryPolicy,params)
    policy.actuations_traj[:] .= params
end

function get_params(policy::TrajectoryPolicy)
    (;actuations_traj) = policy
    vec(actuations_traj)
end

function update!(policy::TrajectoryPolicy)
    # for (i,actuation_traj) = zip(eachindex(policy.itp), eachcol(policy.actuations_traj))
    #     policy.itp[i] = constant_interpolation(policy.time_nodes,actuation_traj)
    # end
end

function actuate!(bot::Robot,policy::TrajectoryPolicy,inst_state::AbstractCoordinatesState)
    (;structure,hub) = bot
    (;t) = inst_state
    (;actuators,coalition,) = hub
    (;actid2sys_actions_idx) = coalition
    (;time_nodes) = policy
    k = searchsortedlast(time_nodes, t)
    for(i,j)  in zip(axes(policy.actuations_traj,2), eachindex(hub.state.u))
        hub.state.u[j] = policy.actuations_traj[k,i]
    end
    foreach(actuators) do actuator
        idx = actid2sys_actions_idx[actuator.id]
        execute!(
            structure,
            actuator,
            (@view hub.state.u[idx])
        )
    end
end

function gen_force_state_jacobian!(∂F∂q,∂F∂q̇,∂F∂s,∂F∂u,policy::TrajectoryPolicy,bot::Robot,inst_state::AbstractCoordinatesState;)
    (;structure, hub) = bot
    (;actuators,coalition) = hub
    (;actid2sys_actions_idx) = coalition
    foreach(actuators) do actuator
        u_idx = actid2sys_actions_idx[actuator.id]
        gen_force_actu_jacobian!(
            (@view ∂F∂u[:,u_idx]),
            structure,
            actuator,
            (@view hub.state.u[u_idx])
        )
    end
end


function vjp_wrt_state(v,policy::TrajectoryPolicy,bot::Robot,num_of_actions,solver,solver_state)
    (;structure) = bot
    (;
        qₖ, pₖ, qₖ₊₁, pₖ₊₁, λₘ, tₖ, tₖ₊₁
    ) = solver_state
    nq = get_num_of_full_coords(structure)
    nλ = get_num_of_cstr(structure)
    ns = get_num_of_aux_var(structure)


    ∂ϕ∂qₖᵀ = ∂ϕ∂qₖ₊₁ᵀ = spzeros(typeof(tₖ),nq)
    ∂ϕ∂pₖᵀ = ∂ϕ∂pₖ₊₁ᵀ = spzeros(typeof(tₖ),nq)
    ∂ϕ∂λᵀ = spzeros(typeof(tₖ),nλ)
    ∂ϕ∂sₖᵀ = spzeros(typeof(tₖ),ns)
    ∂ϕ∂sₖ₊₁ᵀ = spzeros(typeof(tₖ),ns)

    ∂ϕ∂qₖᵀ,  ∂ϕ∂qₖ₊₁ᵀ, ∂ϕ∂pₖᵀ, ∂ϕ∂pₖ₊₁ᵀ, ∂ϕ∂λᵀ, ∂ϕ∂sₖᵀ, ∂ϕ∂sₖ₊₁ᵀ
end

function accumulate_param_grad!(grad_storage, policy::TrajectoryPolicy, v_total, solver_state, bot)
    (; tₘ ) = solver_state
    (; time_nodes, num_of_actions ) = policy
    n_steps = size(policy.actuations_traj, 1)
    
    k = searchsortedlast(time_nodes, tₘ)
    
    # Store gradients in a column-major vector as expected by `get_params` / `set_params!`
    for i in eachindex(v_total)
        idx = (i - 1) * n_steps + k
        grad_storage[idx] += v_total[i]
    end
end
