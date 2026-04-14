# For now, the objective is simply the sum of errors
#done weighted sum
#done (partly) user-defined objective
"""
    cost!(bot::Robot, objt::AbstractObjective, inst_state::AbstractCoordinatesState, u; mode=:trajectory)

Compute the cost for a robot system given the current state and control input at a single time point.

# Arguments
- `bot::Robot`: The robot system
- `objt::AbstractObjective`: The objective function containing weights for gauges and actuators
- `inst_state::AbstractCoordinatesState`: Instantaneous state of the robot system
- `u`: Control inputs
- `mode::Symbol`: Cost mode (:trajectory or :terminal) for weights

# Returns
The computed cost value
"""
function cost!(bot::Robot,objt::AbstractObjective, inst_state::AbstractCoordinatesState,u;mode=:trajectory)
    (; structure, hub) = bot
    (;error_gauges,actuators,coalition) = hub
    (; actid2sys_actions_idx, ) = coalition
    T = get_numbertype(structure)
    update!(structure,inst_state)
    ϕ = Ref(zero(T))
    if mode == :trajectory
        error_gauges_weights = objt.trajectory_error_gauges_weights
        actuators_weights    = objt.trajectory_actuators_weights
    elseif mode == :terminal
        error_gauges_weights = objt.terminal_error_gauges_weights
        actuators_weights    = objt.terminal_actuators_weights
    else
        error("Invalid mode")
    end
    foreach(error_gauges) do gauge
        ϕ[] += error_gauges_weights[gauge.id]*measure(structure,gauge)
        # @show gauge.id ϕ[]
    end
    foreach(actuators) do actuator
        act_idx = actid2sys_actions_idx[actuator.id]
        ϕ[] += actuators_weights[actuator.id]*measure(
            structure,actuator,
            (@view u[act_idx])
        )
    end
    ϕ[]
end

"""
    cost_gradient!(bot::Robot, policy::AbstractPolicy, objt::AbstractObjective, field::AbstractField=NoField(); mode=:trajectory)

Compute the gradient of the cost function with respect to the robot's **current** state variables, control inputs, and parameters.

# Arguments
- `bot::Robot`: The robot system
- `policy::AbstractPolicy`: The control policy
- `objt::AbstractObjective`: The objective function
- `field::AbstractField`: The field applied to the system (default: NoField())
- `mode::Symbol`: The mode of operation (:trajectory or :terminal)

# Returns
A tuple containing the gradients with respect to:
- `∂ϕ∂qᵀ`: Generalized coordinates
- `∂ϕ∂q̇ᵀ`: Generalized velocities  
- `∂ϕ∂pᵀ`: Generalized momenta
- `∂ϕ∂uᵀ`: Control inputs
- `∂ϕ∂θᵀ`: Policy parameters
- `∂ϕ∂cᵀ`: System parameters
"""
function cost_gradient!(bot::Robot, policy::AbstractPolicy, env::AbstractEnvironment, objt::AbstractObjective; 
        mode=:trajectory
    )
    (;structure) = bot
    (;q,q̇,p,s,t) = structure.state.system
    (;u) = bot.hub.state
    T = get_numbertype(bot)
    ∂ϕ∂qᵀ = zeros(T,length(q))
    ∂ϕ∂q̇ᵀ = zeros(T,length(q̇))
    ∂ϕ∂pᵀ = zeros(T,length(p))
    ∂ϕ∂uᵀ = zeros(T,length(u))
    ∂ϕ∂sᵀ = zeros(T,length(s))
    cost_gradient!(
        CostGradient(
            ∂ϕ∂qᵀ, ∂ϕ∂q̇ᵀ, ∂ϕ∂pᵀ, ∂ϕ∂uᵀ, 
            ∂ϕ∂sᵀ, zeros(T,0), zeros(T,0) # ∂ϕ∂θᵀ, ∂ϕ∂cᵀ
        ),Zhong06JacobianWorkspace(bot),
        bot, policy, env, objt,
        structure.state.system,u;mode
    )
    ∂ϕ∂qᵀ, ∂ϕ∂q̇ᵀ, ∂ϕ∂pᵀ, ∂ϕ∂uᵀ, ∂ϕ∂sᵀ
end

function cost_gradient!(gradient::CostGradient,workspace,
        bot::Robot, policy::AbstractPolicy, env::AbstractEnvironment, objt::AbstractObjective,
        inst_state::AbstractCoordinatesState, u; mode=:trajectory
    )
    (;structure,hub) = bot
    (;error_gauges,actuators,coalition) = hub
    (;num_of_full_coords) = structure.connectivity
    (; num_of_error_gauges, num_of_actuators, num_of_actions,
        actid2sys_actions_idx, 
        error_gauge_id2sys_errors_idx,
    ) = coalition
    (;
        trajectory_error_gauges_weights, 
        trajectory_actuators_weights, 
        terminal_error_gauges_weights, 
        terminal_actuators_weights, 
    ) = objt
    T = get_numbertype(structure)
    ∂g∂q = workspace.cost_∂g∂q
    ∂g∂q̇ = workspace.cost_∂g∂q̇
    ∂g∂s = workspace.cost_∂g∂s
    ∂g∂u = workspace.cost_∂g∂u
    tmp_vec = workspace.cost_tmp_vec

    @assert size(∂g∂q) == (num_of_error_gauges, num_of_full_coords)
    @assert size(∂g∂q̇) == (num_of_error_gauges, num_of_full_coords)
    @assert size(∂g∂u) == (num_of_actuators, num_of_actions)
    @assert length(tmp_vec) == num_of_full_coords

    fill!(∂g∂q, zero(T))
    fill!(∂g∂q̇, zero(T))
    fill!(∂g∂s, zero(T))
    fill!(∂g∂u, zero(T))
    structure.state.system.t  = inst_state.t
    structure.state.system.q .= inst_state.q
    structure.state.system.q̇ .= inst_state.q̇
    structure.state.system.s .= inst_state.s
    clear_forces!(structure)
    if u isa Nothing
        actuate!(bot,policy,inst_state)
    end
    stretch!(structure)
    lazy_update_bodies!(structure)
    update_apparatuses!(structure)
    apply_field!(structure, env.field)
    assemble_forces!(inst_state, structure)
    update_inertia_cache!(structure)
    (;M⁻¹, ∂M⁻¹p∂q) = structure.cache.system
    foreach(error_gauges) do gauge
        gau_idx = error_gauge_id2sys_errors_idx[gauge.id]
        measure_gradient!(
            (@view ∂g∂q[gau_idx,:]),
            (@view ∂g∂q̇[gau_idx,:]),
            (@view ∂g∂s[gau_idx,:]),
            structure,gauge,
            workspace.gauge_workspaces[gauge.id]
        )
        # @show gauge.id gau_idx ∂g∂q[gau_idx,:] ∂g∂q̇[gau_idx,:] ∂g∂s[gau_idx,:]
    end

    foreach(actuators) do actuator
        act_idx = actid2sys_actions_idx[actuator.id]
        measure_gradient!(
            (@view ∂g∂u[[actuator.id], act_idx]),
            structure,actuator,
            (@view u[act_idx])
        )
    end

    if mode == :trajectory
        error_gauges_weights = trajectory_error_gauges_weights
        actuators_weights    = trajectory_actuators_weights
    elseif mode == :terminal
        error_gauges_weights = terminal_error_gauges_weights
        actuators_weights    = terminal_actuators_weights
    else
        error("Invalid mode")
    end

    (;
        ∂ϕ∂qᵀ, ∂ϕ∂q̇ᵀ, ∂ϕ∂pᵀ, ∂ϕ∂uᵀ,∂ϕ∂sᵀ
    ) = gradient

    mul!(∂ϕ∂qᵀ, transpose(∂g∂q), error_gauges_weights, one(T), zero(T))
    mul!(∂ϕ∂q̇ᵀ, transpose(∂g∂q̇), error_gauges_weights, one(T), zero(T))
    mul!(tmp_vec, transpose(∂M⁻¹p∂q), ∂ϕ∂q̇ᵀ, one(T), zero(T))
    #TODO ∂ϕ∂qᵀ .+= tmp_vec
    mul!(∂ϕ∂pᵀ, transpose(M⁻¹), ∂ϕ∂q̇ᵀ, one(T), zero(T))
    mul!(∂ϕ∂uᵀ, transpose(∂g∂u), actuators_weights, one(T), zero(T))
    mul!(∂ϕ∂sᵀ, transpose(∂g∂s), error_gauges_weights, one(T), zero(T))
end

function cost_hessian!(bot::Robot, policy::AbstractPolicy, env::AbstractEnvironment, objt::AbstractObjective; mode=:trajectory)
    (;structure) = bot
    (;q,q̇,p,t) = structure.state.system
    (;u) = bot.hub.state
    T = get_numbertype(bot)

    ∂ϕ∂qᵀ∂q = zeros(T,length(q),length(q))
    ∂ϕ∂q̇ᵀ∂q̇ = zeros(T,length(q̇),length(q̇))
    ∂ϕ∂pᵀ∂p = zeros(T,length(p),length(p))
    ∂ϕ∂qᵀ∂p = zeros(T,length(q),length(p))

    ∂ϕ∂qᵀ∂u = zeros(T,length(q),length(u))
    ∂ϕ∂q̇ᵀ∂u = zeros(T,length(q),length(u))
    ∂ϕ∂pᵀ∂u = zeros(T,length(p),length(u))
    ∂ϕ∂uᵀ∂u = zeros(T,length(u),length(u))

    cost_hessian!(
        CostHessian(
            ∂ϕ∂qᵀ∂q, ∂ϕ∂q̇ᵀ∂q̇, ∂ϕ∂pᵀ∂p, ∂ϕ∂qᵀ∂p, 
            ∂ϕ∂qᵀ∂u, ∂ϕ∂q̇ᵀ∂u, ∂ϕ∂pᵀ∂u, ∂ϕ∂uᵀ∂u
        ),
        bot,policy,env,objt,structure.state.system,u;mode
    )
    ∂ϕ∂qᵀ∂q, ∂ϕ∂q̇ᵀ∂q̇, ∂ϕ∂pᵀ∂p, ∂ϕ∂qᵀ∂p, 
    ∂ϕ∂qᵀ∂u, ∂ϕ∂q̇ᵀ∂u, ∂ϕ∂pᵀ∂u, ∂ϕ∂uᵀ∂u
end

"""
    cost_hessian!(hessian::CostHessian, bot::Robot, policy::AbstractPolicy, env::AbstractEnvironment, objt::AbstractObjective, inst_state::AbstractCoordinatesState, u)

Compute the Hessian of the cost function w.r.t. the coordinates `q` and `p` for a given robot `bot` at time `t`.
This function updates the provided matrices in `hessian` in place.
"""
function cost_hessian!(
        hessian::CostHessian,
        bot::Robot, policy::AbstractPolicy, env::AbstractEnvironment, objt::AbstractObjective, 
        inst_state::AbstractCoordinatesState, u;
        mode=:trajectory
    )
    (;structure,hub) = bot
    (;error_gauges,actuators,coalition) = hub
    (;num_of_full_coords) = structure.connectivity
    (;num_of_error_gauges, num_of_actuators, error_gauge_id2sys_errors_idx, actid2sys_actions_idx, num_of_actions,) = coalition
    (;
        trajectory_error_gauges_weights, 
        trajectory_actuators_weights, 
        terminal_error_gauges_weights, 
        terminal_actuators_weights, 
    ) = objt
    T = get_numbertype(structure)
    structure.state.system.t = inst_state.t
    structure.state.system.q .= inst_state.q
    structure.state.system.q̇ .= inst_state.q̇
    structure.state.system.s .= inst_state.s
    clear_forces!(structure)
    stretch!(structure)
    update_bodies!(structure)
    update_apparatuses!(structure)
    apply_field!(structure, env.field)
    assemble_forces!(inst_state, structure)
    M⁻¹ = assemble_M⁻¹(structure)
    ∂M⁻¹p∂q = assemble_∂M⁻¹p∂q(structure)

    ∂²g∂q² = zeros(T,num_of_error_gauges,num_of_full_coords,num_of_full_coords)
    ∂²g∂q̇² = zeros(T,num_of_error_gauges,num_of_full_coords,num_of_full_coords)
    ∂²g∂q̇∂q = zeros(T,num_of_error_gauges,num_of_full_coords,num_of_full_coords)

    ∂²g∂q∂u = zeros(T,num_of_actuators,num_of_full_coords,num_of_actions)
    ∂²g∂q̇∂u = zeros(T,num_of_actuators,num_of_full_coords,num_of_actions)
    ∂²g∂u∂u = zeros(T,num_of_actuators,num_of_actions,num_of_actions)

    if mode == :trajectory
        error_gauges_weights = trajectory_error_gauges_weights
        actuators_weights    = trajectory_actuators_weights
    elseif mode == :terminal
        error_gauges_weights = terminal_error_gauges_weights
        actuators_weights    = terminal_actuators_weights
    else
        error("Invalid mode")
    end

    (;
        ∂ϕ∂qᵀ∂q, ∂ϕ∂q̇ᵀ∂q̇, ∂ϕ∂pᵀ∂p, ∂ϕ∂qᵀ∂p, 
        ∂ϕ∂qᵀ∂u, ∂ϕ∂q̇ᵀ∂u, ∂ϕ∂pᵀ∂u, ∂ϕ∂uᵀ∂u,
    ) = hessian

    foreach(error_gauges) do gauge
        measure_hessians!(
            (@view ∂²g∂q²[gauge.id,:,:]),
            (@view ∂²g∂q̇²[gauge.id,:,:]),
            (@view ∂²g∂q̇∂q[gauge.id,:,:]),
            structure,gauge
        )
    end
    
    foreach(actuators) do actuator
        act_idx = actid2sys_actions_idx[actuator.id]
        measure_hessians!(
            (@view ∂²g∂q∂u[actuator.id,:,:]),
            (@view ∂²g∂q̇∂u[actuator.id,:,:]),
            (@view ∂²g∂u∂u[actuator.id,:,:]),
            structure,actuator,
            (@view u[act_idx])
        )
    end

    ∂ϕ∂qᵀ∂q .= sum(error_gauges_weights[i]*∂²g∂q²[i,:,:] for i in 1:num_of_error_gauges)
    ∂ϕ∂q̇ᵀ∂q̇ .= sum(error_gauges_weights[i]*∂²g∂q̇²[i,:,:] for i in 1:num_of_error_gauges)
    ∂ϕ∂q̇ᵀ∂q = sum(error_gauges_weights[i]*∂²g∂q̇∂q[i,:,:] for i in 1:num_of_error_gauges)

    ∂ϕ∂qᵀ∂q .+= transpose(∂M⁻¹p∂q)*∂ϕ∂q̇ᵀ∂q
    
    ∂ϕ∂pᵀ∂p .= transpose(M⁻¹)*∂ϕ∂q̇ᵀ∂q̇*M⁻¹

    ∂ϕ∂qᵀ∂p .= 0.0

    ∂ϕ∂qᵀ∂u .= 0.0
    ∂ϕ∂q̇ᵀ∂u .= 0.0
    ∂ϕ∂pᵀ∂u .= 0.0
    ∂ϕ∂uᵀ∂u .= sum(actuators_weights[i]*∂²g∂u∂u[i,:,:] for i in 1:num_of_actuators)

end

function cost_gradient_and_hessian!(
        cost_workspace, workspace, 
        bot::Robot, policy::AbstractPolicy, env::AbstractEnvironment, objt::AbstractObjective, 
        inst_state::AbstractCoordinatesState, u;
        mode=:trajectory
    )
    cost_gradient!(
        cost_workspace.gradient, workspace, 
        bot, policy, env, objt,
        inst_state, u; mode
    )
    cost_hessian!(
        cost_workspace.hessian,
        bot, policy, env, objt,
        inst_state, u; mode
    )
    nq = get_num_of_full_coords(bot.structure)
    cost_workspace.∂ϕ∂xᵀ .= 0
    cost_workspace.∂ϕ∂xᵀ[   1: nq] .= cost_workspace.gradient.∂ϕ∂qᵀ
    cost_workspace.∂ϕ∂xᵀ[nq+1:2nq] .= cost_workspace.gradient.∂ϕ∂pᵀ
    ns = get_num_of_aux_var(bot.structure)
    nλ = get_num_of_cstr(bot.structure)
    cost_workspace.∂ϕ∂xᵀ[2nq+nλ+1:2nq+nλ+ns] .= cost_workspace.gradient.∂ϕ∂sᵀ
    cost_workspace.∂ϕ∂xᵀ∂x .= 0
    cost_workspace.∂ϕ∂xᵀ∂x[   1: nq,    1: nq] .= cost_workspace.hessian.∂ϕ∂qᵀ∂q
    cost_workspace.∂ϕ∂xᵀ∂x[nq+1:2nq,    1: nq] .= transpose(cost_workspace.hessian.∂ϕ∂qᵀ∂p)
    cost_workspace.∂ϕ∂xᵀ∂x[   1: nq, nq+1:2nq] .= cost_workspace.hessian.∂ϕ∂qᵀ∂p
    cost_workspace.∂ϕ∂xᵀ∂x[nq+1:2nq, nq+1:2nq] .= cost_workspace.hessian.∂ϕ∂pᵀ∂p
end

"""
    cost!(bot::Robot, objt::AbstractObjective)

Compute the total cost over the entire trajectory for a robot system.

Running cost is computed at the **midpoint** of each time step.

# Arguments
- `bot::Robot`: The robot system
- `objt::AbstractObjective`: The objective function containing weights and time scaling factor η

# Returns
The total cost value (trajectory cost + terminal cost)
"""

function cost!(bot::Robot,objt::Objective, solver, dt )
    #trajectory cost
    ηs = get_trajectory_cost_weights(objt, bot.traj.t, dt)
    ϕt = [
        begin
            uₖ = control_stateₖ.u
            ηs[i]*cost!(bot,objt,stateₖ,uₖ;mode=:trajectory)
        end
        for (i,(stateₖ, control_stateₖ)) in enumerate(zip(
            bot.traj[begin:end-1],
            bot.control_traj[begin:end-1]
        ))
    ]
    #terminal cost
    ϕf = cost!(bot,objt,
        bot.traj[end],
        bot.control_traj[end].u;  # should be end-1
        mode=:terminal
    )
    sum(ϕt) + ϕf
end
