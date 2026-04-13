@concrete struct Direct_Sensitivity_Zhong06_Constant_Mass_Cache <: AbstractZhong06Cache
    solver
    bot
    policy
    env
    consts
    jacobian_workspace
    ηs
    Jac_state
    Jac_action
    Jac_control_params
    traj_cost_gradients_wrt_state
    traj_cost_hessians_wrt_state
    traj_cost_gradients_wrt_action
    traj_cost_hessians_wrt_action
    term_cost_gradient_wrt_state
    term_cost_hessian_wrt_state
    term_cost_gradient_wrt_action
    term_cost_hessian_wrt_action
end

function generate_cache(
        prob::DynamicsProblem,
        solver::DirectDynamicsSensitivitySolver{
            <:DynamicsSolver{<:Zhong06},
        },
        has_constant_mass_matrix::Val{true};
        dt,kargs...
    )
    (;bot,policy,env) = prob
    (;structure,hub) = bot
    options = merge(
        (checkpersist=true,), #default
        prob.options,
        solver.forward_solver.options,
    )
    @info "Direct_Sensitivity_Zhong06_Constant_Mass_Cache"

    (;bot,policy,env) = prob
    (;structure,hub) = bot
    options = merge(
        (checkpersist=true,), #default
        prob.options,
        solver.forward_solver.options,
    )
    T = get_numbertype(structure)
    Mₘ = assemble_M(structure)
    M⁻¹ₘ = assemble_M⁻¹(structure)
    ∂Mₘhq̇ₘ∂qₘ = assemble_∂Mq̇∂q(structure)
    consts = Zhong06Constants(bot, policy, structure, norm(Mₘ,Inf), dt)
    (;nq,nλ,nu,nc,nθ) = consts
    ny = 2*nq+nλ
    nx = ny

    ns = consts.ns
    jacobian_workspace = Zhong06JacobianWorkspace(bot)
    
    ηs = get_trajectory_cost_weights(solver.objt, bot.traj.t[begin:end-1], dt)
    Jac_state  = [ zeros(T,nx,nx) for _ in bot.traj[begin+1:end] ]
    Jac_action = [ zeros(T,nx,nu) for _ in bot.traj[begin+1:end] ]
    Jac_control_params  = [ zeros(T,nx,nθ) for _ in bot.traj[begin+1:end] ]
    traj_cost_gradients_wrt_state  = [ zeros(T,nx)    for _ in bot.traj[begin+1:end] ]
    traj_cost_hessians_wrt_state   = [ zeros(T,nx,nx) for _ in bot.traj[begin+1:end] ]
    traj_cost_gradients_wrt_action = [ zeros(T,nu)    for _ in bot.traj[begin+1:end] ]
    traj_cost_hessians_wrt_action  = [ zeros(T,nu,nu) for _ in bot.traj[begin+1:end] ]
    term_cost_gradient_wrt_state  = zeros(T,nx)
    term_cost_hessian_wrt_state   = zeros(T,nx,nx)
    term_cost_gradient_wrt_action = zeros(T,nu)
    term_cost_hessian_wrt_action  = zeros(T,nu,nu)

    Direct_Sensitivity_Zhong06_Constant_Mass_Cache(
        solver,
        bot,policy,env,
        consts,jacobian_workspace,
        ηs,
        Jac_state, 
        Jac_action,
        Jac_control_params,
        traj_cost_gradients_wrt_state,
        traj_cost_hessians_wrt_state,
        traj_cost_gradients_wrt_action,
        traj_cost_hessians_wrt_action,
        term_cost_gradient_wrt_state,
        term_cost_hessian_wrt_state,
        term_cost_gradient_wrt_action,
        term_cost_hessian_wrt_action
    )
end

function solve!(simulator::Simulator,
                forward_cache::Zhong06_Constant_Mass_Cache,
                solver_cache::Direct_Sensitivity_Zhong06_Constant_Mass_Cache;
                dt,ftol=1e-14,verbose=false,maxiters=50,
                progress=true,exception=true)
    (;prob,controller,tspan,restart,) = simulator
    (;hub,traj,control_traj) = prob.bot
    num_of_sim_timestep = simulator.totalstep
    (;
        solver,
        bot,policy,env,
        consts,jacobian_workspace,
        ηs,
        Jac_state, 
        Jac_action,
        Jac_control_params,
        traj_cost_gradients_wrt_state,
        traj_cost_hessians_wrt_state,
        traj_cost_gradients_wrt_action,
        traj_cost_hessians_wrt_action,
        term_cost_gradient_wrt_state,
        term_cost_hessian_wrt_state,
        term_cost_gradient_wrt_action,
        term_cost_hessian_wrt_action
    ) = solver_cache
    (;structure) = bot
    T = get_numbertype(structure)
    (;nq,nλ,ns,nu,nθ,nc) = consts
    nx = ny = 2nq+nλ
    
    # Create cost workspace
    cost_workspace = CostGradientHessianWorkspace(T, nx, nq, ns, nu, nθ, nc)
    
    step = 0
    h = dt
    
    # Create Jacobian blocks
    Jacᵏ⁺¹ₖ_backup = zeros(T,ny,ny)
    jac_blocks = Zhong06JacobianBlocks(T, nx, nu, nc, Jacᵏ⁺¹ₖ_backup)

    prog = Progress(num_of_sim_timestep; dt=1.0, enabled=progress)
    
    uN   = control_traj.u[num_of_sim_timestep+1]

    cost_gradient_and_hessian!(
        cost_workspace, jacobian_workspace, 
        bot, policy, env, solver.objt,
        traj[num_of_sim_timestep+1], uN, ;
        mode=:terminal
    )

    term_cost_gradient_wrt_state .= cost_workspace.∂ϕ∂xᵀ
    term_cost_hessian_wrt_state  .= cost_workspace.∂ϕ∂xᵀ∂x
    term_cost_gradient_wrt_action .= cost_workspace.gradient.∂ϕ∂uᵀ
    term_cost_hessian_wrt_action  .= cost_workspace.hessian.∂ϕ∂uᵀ∂u

    for timestep = (num_of_sim_timestep):-1:1
        #---------Step k Control-----------
        @debug "Timestep $timestep Begin" time = traj.t[timestep]
        #---------Step k Control-----------
        jac_blocks.Jacᵏ⁺¹ₖ_backup .= jac_blocks.Jacᵏ⁺¹ₖ
        uₖ = control_traj.u[timestep]
        
        state_k = traj[timestep]
        state_kp1 = traj[timestep+1]
        state_mid = deepcopy(state_kp1)
        solver_state = Zhong06SolverState(
            state_k, state_kp1, state_mid,
            h
        )
        interpolate!(solver_state)
        Zhong06_constant_mass_Jac!(
            jac_blocks,
            solver_state,
            bot,env,policy,forward_cache,
            consts,jacobian_workspace
        )

        ∂ϕ∂uᵀ = zeros(T,nu)
        vjp_wrt_state(∂ϕ∂uᵀ,policy,bot,nu,forward_cache,solver_state)
    
        Jacᵏ⁺¹ₘθ_T = zeros(T, nθ, nx)
        for i in 1:nx
            v_total = jac_blocks.Jacᵏ⁺¹ₘu[i, :] 
            v_storage = view(Jacᵏ⁺¹ₘθ_T, :, i)
            accumulate_param_grad!(v_storage, policy, v_total, solver_state, bot)
        end
        Jacᵏ⁺¹ₘθ = Matrix(transpose(Jacᵏ⁺¹ₘθ_T))

        luJac = lu(jac_blocks.Jacᵏ⁺¹ₖ₊₁)
        Jac_state[timestep]          .= -(luJac\jac_blocks.Jacᵏ⁺¹ₖ)
        Jac_action[timestep]         .= -(luJac\jac_blocks.Jacᵏ⁺¹ₘu)
        Jac_control_params[timestep] .= -(luJac\Jacᵏ⁺¹ₘθ)

        cost_gradient_and_hessian!(
            cost_workspace, jacobian_workspace, 
            bot, policy, env, solver.objt,
            state_k, uₖ, ;
            mode=:trajectory
        )

        traj_cost_gradients_wrt_state[timestep] .= ηs[timestep]*cost_workspace.∂ϕ∂xᵀ
        traj_cost_hessians_wrt_state[timestep]  .= ηs[timestep]*cost_workspace.∂ϕ∂xᵀ∂x
        traj_cost_gradients_wrt_action[timestep] .= ηs[timestep]*cost_workspace.gradient.∂ϕ∂uᵀ
        traj_cost_hessians_wrt_action[timestep]  .= ηs[timestep]*cost_workspace.hessian.∂ϕ∂uᵀ∂u

        #---------Step k finisher-----------
        step += 1
        #---------Step k finisher-----------
        @debug "Timestep $timestep End  "
        next!(prog)
    end
end
