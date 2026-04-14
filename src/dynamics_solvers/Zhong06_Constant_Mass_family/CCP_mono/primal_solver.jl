# ============================================================================
# Cache and Setup
# ============================================================================

struct Zhong06_CCP_Constant_Mass_Mono_Cache{
        T,
        RobotType,
        PolicyType,
        FieldType,
        EnvType,
        OptionsType,
        VT,
        MT
    } <: AbstractZhong06Cache
    bot::RobotType
    policy::PolicyType
    field::FieldType
    env::EnvType
    jacobian_workspace::Zhong06JacobianWorkspace{T}
    consts::Zhong06Constants{T}
    options::OptionsType
    cost_workspace::CostGradientHessianWorkspace{VT,MT}
    direct_sensitivity_workspace::DirectSensitivityWorkspace{T}
end

function generate_cache(
        prob::DynamicsProblem{
            <:Robot,
            <:AbstractPolicy,
            <:AbstractEnvironment,
            <:AbstractObjective,
            <:RestitutionFrictionCombined{NewtonRestitution,CoulombFriction}
        },
        solver::DynamicsSolver{
            <:Zhong06,
            <:AbstractBodySolver,
            <:AbstractApparatusSolver,
            <:MonolithicContactSolver
        },
        ::Val{true};
        dt,compute_sensitivities=true,kargs...
    )  
    @info "Zhong06_CCP_Constant_Mass_Mono_Cache" compute_sensitivities
    (;bot,policy,env) = prob
    (;structure,hub) = bot
    (;field) = env
    options = merge(
        (checkpersist=true,), #default
        (compute_sensitivities=compute_sensitivities,),
        prob.options,
        solver.options,
    )
    Mₘ = assemble_M(structure)
    M⁻¹ₘ = assemble_M⁻¹(structure)
    ∂Mₘhq̇ₘ∂qₘ = assemble_∂Mq̇∂q(structure)
    mr = norm(Mₘ,Inf)
    mass_norm = mr
    consts = Zhong06Constants(bot, policy, structure, mass_norm, dt)
    (;nq, nλ, ns, nu, nc, nθ, n1, n2, n3) = consts
    T = get_numbertype(structure)
    
    jacobian_workspace = Zhong06JacobianWorkspace(bot)   
    
    direct_sensitivity_workspace = DirectSensitivityWorkspace(T, n3, nu)
    
    # Allocate cost gradient and Hessian workspace
    cost_workspace = CostGradientHessianWorkspace(T, n3, nq, ns, nu, nθ, nc)
    
    Zhong06_CCP_Constant_Mass_Mono_Cache(bot,
        policy,
        field,
        env,
        jacobian_workspace,
        consts,
        options,
        cost_workspace,
        direct_sensitivity_workspace,
    )
end


# ============================================================================
# Predictor-Corrector Helper Functions
# ============================================================================

"""
Compute predictor step in predictor-corrector interior point method.
$(TYPEDSIGNATURES)
"""
function compute_predictor_step(contact_vars::ContactVariables, na)
    # Duality measure
    μ = [(transpose(contact_vars.Γ_split[i])*contact_vars.Λ_split[i])/3 for i in 1:na]
    
    # Maximum step length
    αp_Λ = find_cone_step_length(contact_vars.Λ_split, contact_vars.ΔΛp_split, contact_vars.J)
    αp_Γ = find_cone_step_length(contact_vars.Γ_split, contact_vars.ΔΓp_split, contact_vars.J)
    αpmax = min(αp_Λ, αp_Γ)
    αp = min(one(αpmax), 0.99αpmax)
    
    return μ, αp
end

"""
Compute centering parameter σ and corrector residual
"""
function compute_centering_parameter(contact_vars::ContactVariables, μ, αp, na)
    # Predicted values
    contact_vars.Λp .= contact_vars.Λ .+ αp.*contact_vars.ΔΛp
    contact_vars.Γp .= contact_vars.Γ .+ αp.*contact_vars.ΔΓp
    
    # Centering parameter
    μp = [(transpose(contact_vars.Γp_split[i])*contact_vars.Λp_split[i])/3 for i in 1:na]
    σ = [(μp[i]/μ[i])^3 for i in 1:na]
    
    # Check for degenerate cases
    if any(isnan.(σ)) || any(iszero.(μ))
        return nothing, nothing, nothing
    end
    
    τ = σ .* μp
    return σ, μp, τ
end

# ============================================================================
# Newton Iteration for Contact Case
# ============================================================================

"""
Perform single Newton iteration with predictor-corrector for contact case (na > 0)

Note: contact_vars.Δxp should already be computed before calling this function
"""
function newton_iteration_with_contact!(
        workspace::NewtonWorkspace,
        contact_vars::ContactVariables,
        solver_state::Zhong06SolverState,
        solver_cache, contact_cache,
        lu𝐉, na, n2, nΛ, ftol
    )
        # Predictor step: compute duality measure and step length
    μ, αp = compute_predictor_step(contact_vars, na)
    
    # Compute centering parameter σ and corrector target τ
    result = compute_centering_parameter(contact_vars, μ, αp, na)
    
    if result === nothing || result[1] === nothing
        return false, Inf, μ, nothing, nothing  # Signal degenerate case
    end
    
    σ, μp, τ = result
    
    # Corrector residual modification
    𝐫𝐞𝐬_c_split = -τ .* contact_vars.𝐞_split
    workspace.Res[n2+nΛ+1:n2+2nΛ] .+= reduce(vcat, 𝐫𝐞𝐬_c_split)
    
    normRes = norm(workspace.Res)
    
    # Check convergence
    if normRes < ftol
        return true, normRes, μ, μp, σ  # Converged
    end
    
    # Corrector direction
    contact_vars.Δxc .= lu𝐉 \ (-workspace.Res)
    
    # Compute step length for corrector
    α_Λ = find_cone_step_length(contact_vars.Λ_split, contact_vars.ΔΛc_split, contact_vars.J)
    α_Γ = find_cone_step_length(contact_vars.Γ_split, contact_vars.ΔΓc_split, contact_vars.J)
    αmax = min(α_Λ, α_Γ)
    α = min(1, 0.99αmax)
    
    # Prepare for line search
    workspace.Δx .= α .* contact_vars.Δxc
    
    # Create line search functions with shared setup (more efficient)
    call_res! = (ws) -> begin
        compute_primal_residual!(ws, solver_state, solver_cache, contact_cache)
        ws.Res[n2+nΛ+1:n2+2nΛ] .+= reduce(vcat, 𝐫𝐞𝐬_c_split)
    end
    call_jac! = (ws) -> compute_primal_jacobian!(ws, solver_state, solver_cache, contact_cache)
    
    ls_workspace = NewtonWorkspace(
        workspace.xₖ, workspace.Res, workspace.Jac,
        workspace.Δx, workspace.x, workspace.𝐰,
        workspace.∂Γ∂x
    )
    ϕ, dϕ, ϕdϕ = create_line_search_functions(workspace, ls_workspace, call_res!, call_jac!)
    
    # Perform backtracking line search
    fx, dϕ_0 = ϕdϕ(0.0)
    β, fx = LS.BackTracking(order=3)(ϕ, dϕ, ϕdϕ, 1.0, fx, dϕ_0)
    
    # Update solution
    normRes = sqrt(2fx)
    workspace.x .+= β .* workspace.Δx
    
    return false, normRes, μ, μp, σ  # Not converged yet, continue iteration
end

# ============================================================================
# Sensitivity Computation
# ============================================================================

"""
Compute direct sensitivities after convergence at timestep.

Computes:
- Jacobian w.r.t. state at k: ∂x_{k+1}/∂x_k
- Jacobian w.r.t. action: ∂x_{k+1}/∂u_{k+1}
- Jacobian w.r.t. control parameters: ∂x_{k+1}/∂θ
- Cost gradients and Hessians w.r.t. state and action
"""
function compute_timestep_sensitivities!(
        solver_cache, forward_cache, solver_state, uₖ, contact_cache, objective, η_at_timestep
    )
    (;
        bot, policy, field, env,
        jacobian_workspace,
        consts,
        cost_workspace,
        direct_sensitivity_workspace,
    ) = solver_cache
    dsw = direct_sensitivity_workspace
    (;structure, hub) = bot
    (;h, nq̌, nλ, nu, nc, nθ, n1, n2, n3) = consts
    (;na) = contact_cache
    
    T = get_numbertype(structure)
    nΛ = 3na
    nx = n3 + 2nΛ
    
    # Allocate Jacobian blocks for this timestep
    jac_blocks = Zhong06JacobianBlocks(T, nx, nu, nc, zeros(T, 0, 0))
    
    Jacᵏ⁺¹ₘθ = zeros(T, nx, nθ)
    
    # Compute Jacobians using shared function
    compute_zhong06_jacobian_blocks!(
        jac_blocks,
        jacobian_workspace,
        solver_state,
        consts,
        contact_cache,
        bot, policy, field, forward_cache
    )
    
    # Compute control parameter Jacobian row by row using VJP
    Jacᵏ⁺¹ₘθ_T = zeros(T, nθ, nx)
    for i in 1:nx
        v_total = jac_blocks.Jacᵏ⁺¹ₘu[i, :] # size nu
        v_storage = view(Jacᵏ⁺¹ₘθ_T, :, i)
        accumulate_param_grad!(v_storage, policy, v_total, solver_state, bot)
    end
    Jacᵏ⁺¹ₘθ = Matrix(transpose(Jacᵏ⁺¹ₘθ_T))
    
    # Solve for sensitivity matrices
    luJac = lu(jac_blocks.Jacᵏ⁺¹ₖ₊₁)
    push!(dsw.Jac_state, -(luJac \ jac_blocks.Jacᵏ⁺¹ₖ))
    push!(dsw.Jac_action, -(luJac \ jac_blocks.Jacᵏ⁺¹ₘu))
    push!(dsw.Jac_control_params, -(luJac \ Jacᵏ⁺¹ₘθ))
    
    # Use cost workspace (clear it first)
    clear!(cost_workspace)
    
    # Compute cost gradient and Hessian for trajectory step
    cost_gradient_and_hessian!(
        cost_workspace, jacobian_workspace, 
        bot, policy, env, objective,
        solver_state.state_k, uₖ,;
        mode=:trajectory
    )
    
    # Store cost gradients and Hessians (scaled by η and h)
    push!(dsw.traj_cost_gradients_wrt_state, η_at_timestep .* copy(cost_workspace.∂ϕ∂xᵀ))
    push!(dsw.traj_cost_hessians_wrt_state, η_at_timestep .* copy(cost_workspace.∂ϕ∂xᵀ∂x))
    push!(dsw.traj_cost_gradients_wrt_action, η_at_timestep .* copy(cost_workspace.gradient.∂ϕ∂uᵀ))
    push!(dsw.traj_cost_hessians_wrt_action, η_at_timestep .* copy(cost_workspace.hessian.∂ϕ∂uᵀ∂u))
    
    return nothing
end

# ============================================================================
# Main Solver
# ============================================================================


function solve!(sim::Simulator,solver_cache::Zhong06_CCP_Constant_Mass_Mono_Cache;
        dt,
        ftol=1e-14,xtol=ftol,
        verbose=false,verbose_contact=false,
        maxiters=50,max_restart=3,
        progress=true,exception=true
    )
    @info "Solving Zhong06_CCP_Constant_Mass_Mono_Cache"
    (;prob,controller,tspan,restart,totalstep) = sim
    (;bot,policy,env) = prob
    # Get objective and η function
    (;objt) = prob
    (;η) = objt
    (;hub,structure,traj,control_traj,contacts_traj,contact_caches_traj) = bot
    (;jacobian_workspace, consts, options) = solver_cache
    (;Mₘ, M⁻¹ₘ) = jacobian_workspace
    (;h, mass_norm, nq, nλ) = consts
    dt = h  # Use h from consts (should match dt parameter)
    (;compute_sensitivities) = options
    (;activated_bits) = structure.contacts_related
    q0 = traj.q[begin]
    state_mid = deepcopy(traj[begin])
    activate_frictional_contacts!(structure,env,q0;checkpersist=true,)
 
    T = eltype(q0)
    nx = nq + nλ

    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Time Step k Control-----------
        @debug "Timestep $timestep Begin" time = traj.t[timestep]
        #---------Time Step k Control-----------
        cₖ = contacts_traj[timestep]
        cₖ₊₁ = contacts_traj[timestep+1]
        qₖ = traj.q[timestep]
        q̇ₖ = traj.q̇[timestep]
        pₖ = traj.p[timestep]
        # λₖ = traj.λ[timestep]
        tₖ = traj.t[timestep]
        tₖ₊₁ = traj.t[timestep+1]
        qₖ₊₁ = traj.q[timestep+1]
        q̇ₖ₊₁ = traj.q̇[timestep+1]
        pₖ₊₁ = traj.p[timestep+1]
        λₘ = traj.λ[timestep+1]
        uₖ = control_traj.u[timestep]
        qˣ = qₖ .+ dt./2 .*q̇ₖ
        qₖ₊₁ .= qₖ .+ dt .*q̇ₖ
        q̇ₖ₊₁ .= q̇ₖ
        contact_cache = activate_frictional_contacts!(structure,env,qˣ;checkpersist=true,)
        (;na) = contact_cache
        (;L,Lv) = contact_cache
        isconverged = false
        normRes = typemax(T)
        iteration_break = 0
        isconverged = false
        nΛ = 3na
        n1 = nq
        n2 = n1 + nλ
        nx = n2 + 2nΛ
        
        # Create Newton workspace using standalone constructor
        workspace = NewtonWorkspace(T, nx, nΛ, n2)
        
        # Create contact variables using standalone constructor
        contact_vars = ContactVariables(workspace.x, n2, na, nx)
        
        # Debugging workspace (not part of NewtonWorkspace)
        Jac_fd = zeros(T, nx, nx)

        state_k = MonoContactCoordinatesState(traj[timestep],contact_vars.Λ,contact_vars.Γ)
        state_kp1 = MonoContactCoordinatesState(traj[timestep+1],contact_vars.Λ,contact_vars.Γ)
        # Note: In solver_state, qₖ corresponds (previous/known), qₖ₊₁ to (current/unknown)
        solver_state = Zhong06SolverState(
            state_k, state_kp1, 
            MonoContactCoordinatesState(state_mid,contact_vars.Λ,contact_vars.Γ), 
            h,
        )
        interpolate!(solver_state)
        get_frictional_directions_and_positions!(structure, contact_cache, state_k, contact_vars.Λ)
        
        restart_count = 0
        Λ_guess = 1.0
        
        # Create unified solver state outside restart loop
        # Note: For primal solver, only qₖ, q̇ₖ, pₖ, tₖ, tₖ₊₁, dt are used
        # Other fields are initialized but not accessed in compute_primal_residual!/jacobian!
        
        while restart_count < max_restart
            @debug "Restart $restart_count Begin" num_of_active_contacts=na
            contact_vars.Λ .= repeat([Λ_guess,0,0],na)
            contact_vars.Γ .= contact_vars.Λ
            workspace.x[      1:nq]          .= qₖ₊₁
            workspace.x[   nq+1:nq+nλ]       .= 0.0
            normRes = typemax(T)
            normRes_last = typemax(T)
            
            for iteration = 1:maxiters
                @debug "Iteration $iteration Begin"
                # Use structured interface for residual and Jacobian computation
                compute_primal_residual!(workspace, solver_state, solver_cache, contact_cache)
                compute_primal_jacobian!(workspace, solver_state, solver_cache, contact_cache)
                # FiniteDiff.finite_difference_jacobian!(Jac_fd,Res_stepk!,workspace.x,Val{:central})
                lu𝐉 = lu(workspace.Jac)
                if na == 0
                    normRes = norm(workspace.Res)
                    if  normRes < ftol
                        isconverged = true
                        iteration_break = iteration-1
                        @debug "Iteration $iteration Break" normRes
                        break
                    elseif normRes == normRes_last
                        isconverged = false
                        iteration_break = iteration-1
                        @debug "Iteration $iteration Stuck" normRes
                        break
                    end
                    workspace.Δx .= lu𝐉\(-workspace.Res)
                    workspace.x .+= workspace.Δx
                else # na!=0
                    # Predictor: compute affine scaling direction
                    contact_vars.Δxp .= lu𝐉\(-workspace.Res)
                    
                    # Call modular helper function for contact iteration
                    converged, normRes, μ, μp, σ = newton_iteration_with_contact!(
                        workspace,
                        contact_vars,
                        solver_state,
                        solver_cache, contact_cache,
                        lu𝐉, na, n2, nΛ, ftol
                    )
                    
                    # Check convergence or termination conditions
                    if converged
                        isconverged = true
                        iteration_break = iteration-1
                        @debug "Iteration $iteration Break" normRes
                        break
                    elseif isnan(normRes) || isinf(normRes) || μ === nothing
                        @debug "Iteration $iteration Degenerate case detected"
                        break
                    elseif normRes == normRes_last
                        iteration_break = iteration-1
                        isconverged = false
                        @debug "Iteration $iteration Stuck" normRes
                        break
                    elseif iteration == maxiters
                        iteration_break = iteration-1
                        isconverged = false
                    end
                    
                end
                @debug "Iteration $iteration End  " normRes normRes_last
                normRes_last = normRes
            end
            if isconverged
                break
            end
            if na > 0
                Λ_guess =  maximum(abs.(contact_vars.Λ[begin:3:end]))
                @debug "Variables" Λ_guess
            end
            @debug "Restart $restart_count End  "
            restart_count += 1
        end
        populate!(solver_state, workspace.x, structure, jacobian_workspace, consts)
        uₖ .= hub.state.u
        Dₘ = contact_cache.Dper
        Dₖ₊₁ = contact_cache.Dimp
        if na != 0
            update_contacts!(cₖ₊₁[activated_bits],cₖ[activated_bits],Dₘ*(qₖ₊₁.-qₖ).+Dₖ₊₁*q̇ₖ₊₁,2*contact_vars.Λ./(mass_norm*dt))
            contact_cache.Λ .= contact_vars.Λ
            contact_cache.Γ .= contact_vars.Γ
            contact_caches_traj[timestep+1] = contact_cache
        end

        if !isconverged
            @warn "Timestep $timestep at time $(tₖ+dt) not converged at $iteration_break Newton iterations (max $maxiters)" normRes restart_count num_of_active_contacts=na
            if exception
                @error "Not converged!"
                break
            else
                # sim.convergence = false
                # break
            end
        end
        
        # Compute sensitivities after convergence if requested
        if compute_sensitivities && isconverged
            η_at_timestep = η(tₖ)
            
            # Update solver_state with converged values
            # (qₖ₊₁, pₖ₊₁, λₘ, q̇ₖ₊₁ were updated after restart loop at lines 493-497)
            interpolate!(solver_state)
            # Contact forces Λ and Γ are views into contact_vars, already updated
            
            compute_timestep_sensitivities!(
                solver_cache, solver_cache, solver_state, uₖ, contact_cache, objt, η_at_timestep
            )
        end
        
        #---------Time Step k finisher-----------
        @debug "Timestep $timestep End  "  normRes restart_count iteration_break num_of_active_contacts=na
        next!(prog)
    end
    
    # Compute terminal cost gradients and Hessians if sensitivities requested
    if compute_sensitivities
        (;
            cost_workspace,
            direct_sensitivity_workspace,
            consts
        ) = solver_cache
        (;nu) = consts
        dsw = direct_sensitivity_workspace
        # Get terminal state
        uN = control_traj.u[totalstep+1]
        
        # Clear cost workspace
        
        (;objt) = prob
        
        cost_gradient_and_hessian!(
            cost_workspace, jacobian_workspace, 
            bot, policy, env, objt,
            traj[totalstep+1], uN;
            mode=:terminal
        )
        
        dsw.term_cost_gradient_wrt_state .= cost_workspace.∂ϕ∂xᵀ
        dsw.term_cost_hessian_wrt_state .= cost_workspace.∂ϕ∂xᵀ∂x
        dsw.term_cost_gradient_wrt_action .= cost_workspace.gradient.∂ϕ∂uᵀ
        dsw.term_cost_hessian_wrt_action .= cost_workspace.hessian.∂ϕ∂uᵀ∂u
    end
    
    bot
end
