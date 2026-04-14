# Include shared Jacobian computation module
# Note: CostGradientHessianWorkspace is defined in Zhong06_CCP_constant_mass_mono_jacobians.jl

struct Direct_Sensitivity_Zhong06_CCP_Constant_Mass_Mono_Cache{directcacheType} <: AbstractZhong06Cache
    cache::directcacheType
end

function generate_cache(
        prob::DynamicsProblem,
        solver::DirectDynamicsSensitivitySolver{
            <:DynamicsSolver{
                <:Zhong06,
                <:AbstractBodySolver,
                <:AbstractApparatusSolver,
                <:MonolithicContactSolver
            },
        },
        has_constant_mass_matrix::Val{true};
        dt,totalstep,kargs...
    )
    (;bot,policy,env) = prob
    (;structure,hub) = bot
    options = merge(
        (checkpersist=true,), #default
        prob.options,
        solver.forward_solver.options,
    )
    @info "Direct_Sensitivity_Zhong06_CCP_Constant_Mass_Mono_Cache"

    (;bot,policy,env) = prob
    (;structure,hub,traj,contact_caches_traj) = bot
    options = merge(
        (checkpersist=true,), #default
        prob.options,
        solver.forward_solver.options,
    )
    T = get_numbertype(structure)
    Mₘ = assemble_M(structure)
    M⁻¹ₘ = assemble_M⁻¹(structure)
    ∂Mₘhq̇ₘ∂qₘ = assemble_∂Mq̇∂q(structure)

    zhong06_constants = Zhong06Constants(bot, policy, structure, norm(Mₘ,Inf), dt)
    (;nq, nλ, nu, nc, ns, nθ, n1, n2, n3) = zhong06_constants
    jacobian_workspace = Zhong06JacobianWorkspace(bot)
    direct_sensitivity_workspace = DirectSensitivityWorkspace(T, n3, nu)
    
    # Allocate cost gradient and Hessian workspace
    cost_workspace = CostGradientHessianWorkspace(T, n3, nq, ns, nu, nθ, nc)
    ηs = get_trajectory_cost_weights(solver.objt, traj.t[begin:end-1], dt)      
    
    cache = @eponymtuple(
        solver,
        cost_workspace,
        direct_sensitivity_workspace,
        jacobian_workspace,
        zhong06_constants,
        ηs
    )
    Direct_Sensitivity_Zhong06_CCP_Constant_Mass_Mono_Cache(cache)
end

function solve!(simulator::Simulator,
                forward_cache::Zhong06_CCP_Constant_Mass_Mono_Cache,
                solver_cache::Direct_Sensitivity_Zhong06_CCP_Constant_Mass_Mono_Cache;
                dt,ftol=1e-14,verbose=false,maxiters=50,
                progress=true,exception=true)
    (;prob,controller,tspan,restart,) = simulator
    (;bot,policy,env) = prob
    (;structure,hub,traj,control_traj,contact_caches_traj) = bot
    num_of_sim_timestep = simulator.totalstep
    (;
        solver,
        cost_workspace,
        direct_sensitivity_workspace,
        jacobian_workspace,
        zhong06_constants,
        ηs
    ) = solver_cache.cache
    dsw = direct_sensitivity_workspace
    # Unpack from jacobian_workspace
    (;Mₘ, M⁻¹ₘ, ∂Mₘhq̇ₘ∂qₘ, Fₘ, ∂Fₘ∂u, ∂Fₘ∂c) = jacobian_workspace
    
    T = get_numbertype(structure)
    zhong06_constants = Zhong06Constants(bot, policy, structure, norm(Mₘ,Inf), dt)
    (;nq, nλ, nu, nc, nθ, n1, n2, n3, h) = zhong06_constants
    nx = ny = n3
    step = 0
    prog = Progress(num_of_sim_timestep; dt=1.0, enabled=progress)
    
    state_mid = deepcopy(traj[begin])
    for timestep = 1:num_of_sim_timestep
        #---------Step k Control-----------
        @debug "Timestep $timestep Begin" time = traj.t[timestep]
        #---------Step k Control-----------
        uₖ = control_traj.u[timestep]
        
        Λₖ₊₁ = contact_caches_traj[timestep+1].Λ
        Γₖ₊₁ = contact_caches_traj[timestep+1].Γ

        nΛ = length(Λₖ₊₁)
        nx = ny = n3 + 2nΛ
        state_k = MonoContactCoordinatesState(traj[timestep], Λₖ₊₁, Γₖ₊₁)
        state_kp1 = MonoContactCoordinatesState(traj[timestep+1], Λₖ₊₁, Γₖ₊₁)
        solver_state = Zhong06SolverState(
            state_k,
            state_kp1,
            MonoContactCoordinatesState(
                state_mid,
                Λₖ₊₁, 
                Γₖ₊₁
            ),
            h
        )

        # Prepare velocities
        solver_state.q̇ₖ₊₁ .= M⁻¹ₘ * solver_state.pₖ₊₁
        solver_state.q̇ₖ   .= M⁻¹ₘ * solver_state.pₖ

        interpolate!(solver_state)
        # nx = n3+2*3contact_caches_traj[timestep+1].na
        jac_blocks = Zhong06JacobianBlocks(T, nx, nu, nc, zeros(T, 0, 0))
        Jacᵏ⁺¹ₘθ = zeros(T, nx, nθ)

        
        η_at_timestep = ηs[timestep]


        push!(dsw.Jac_state, zeros(T,nx,nx))
        push!(dsw.Jac_action, zeros(T,nx,nu))
        push!(dsw.Jac_control_params, zeros(T,nx,nθ))
        push!(dsw.traj_cost_gradients_wrt_state, zeros(T,nx))
        push!(dsw.traj_cost_hessians_wrt_state, zeros(T,nx,nx))
        push!(dsw.traj_cost_gradients_wrt_action, zeros(T,nu))
        push!(dsw.traj_cost_hessians_wrt_action, zeros(T,nu,nu))

        
        # Compute Jacobians using shared function
        compute_zhong06_jacobian_blocks!(
            jac_blocks,
            jacobian_workspace,
            solver_state,
            zhong06_constants,
            contact_caches_traj[timestep+1],
            bot, policy, env.field, forward_cache
        )

        ∂ϕ∂uᵀ = zeros(T,nu)
        vjp_wrt_state(∂ϕ∂uᵀ,policy,bot,nu,forward_cache,solver_state)
    
        Jacᵏ⁺¹ₘθ_T = zeros(T, nθ, nx)
        for i in 1:nx
            v_total = jac_blocks.Jacᵏ⁺¹ₘu[i, :] 
            v_storage = view(Jacᵏ⁺¹ₘθ_T, :, i)
            accumulate_param_grad!(v_storage, policy, v_total, solver_state, bot)
        end
        Jacᵏ⁺¹ₘθ .= transpose(Jacᵏ⁺¹ₘθ_T)
        luJac = lu(jac_blocks.Jacᵏ⁺¹ₖ₊₁)
        dsw.Jac_state[timestep]          .= -(luJac\jac_blocks.Jacᵏ⁺¹ₖ)
        dsw.Jac_action[timestep]         .= -(luJac\jac_blocks.Jacᵏ⁺¹ₘu)
        dsw.Jac_control_params[timestep] .= -(luJac\Jacᵏ⁺¹ₘθ)

        # Clear cost workspace
        clear!(cost_workspace)
        
        cost_gradient_and_hessian!(
            cost_workspace, jacobian_workspace, 
            bot, policy, env, solver.objt,
            state_k, uₖ,;
            mode=:trajectory
        )

        dsw.traj_cost_gradients_wrt_state[timestep][1:n3]     .= η_at_timestep.*cost_workspace.∂ϕ∂xᵀ
        dsw.traj_cost_hessians_wrt_state[timestep][1:n3,1:n3] .= η_at_timestep.*cost_workspace.∂ϕ∂xᵀ∂x
        dsw.traj_cost_gradients_wrt_action[timestep]          .= η_at_timestep.*cost_workspace.gradient.∂ϕ∂uᵀ
        dsw.traj_cost_hessians_wrt_action[timestep]           .= η_at_timestep.*cost_workspace.hessian.∂ϕ∂uᵀ∂u

        #---------Step k finisher-----------
        step += 1
        #---------Step k finisher-----------
        @debug "Timestep $timestep End  "
        next!(prog)
    end

    uN   = control_traj.u[num_of_sim_timestep+1]

    # Clear cost workspace
    clear!(cost_workspace)
    
    cost_gradient_and_hessian!(
        cost_workspace, jacobian_workspace, 
        bot, policy, env, solver.objt,
        traj[num_of_sim_timestep+1], uN,;
        mode=:terminal
    )

    dsw.term_cost_gradient_wrt_state .= cost_workspace.∂ϕ∂xᵀ
    dsw.term_cost_hessian_wrt_state  .= cost_workspace.∂ϕ∂xᵀ∂x
    dsw.term_cost_gradient_wrt_action .= cost_workspace.gradient.∂ϕ∂uᵀ
    dsw.term_cost_hessian_wrt_action  .= cost_workspace.hessian.∂ϕ∂uᵀ∂u
end
