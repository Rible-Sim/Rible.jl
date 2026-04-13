

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
            <:InnerLayerContactSolver
        },
        has_constant_mass_matrix::Val{true};
        dt,kargs...
    )  
    (;bot,policy,env) = prob
    (;structure,hub) = bot
    (;field) = env
    options = merge(
        (checkpersist=true,), #default
        prob.options,
        solver.options,
    )
    @info "Zhong06_CCP_Constant_Mass_Inner_Cache"

    Mₘ = Matrix(assemble_M(structure))
    M⁻¹ₘ = inv(Mₘ)
    
    # Use sparse Mₘ for workspace to satisfy type (though we use dense Mₘ in solver)
    M_sparse = assemble_M(structure)
    M⁻¹_sparse = assemble_M⁻¹(structure)
    ∂Mₘhq̇ₘ∂qₘ = assemble_∂Mq̇∂q(structure)
    
    mr = norm(Mₘ,Inf)
    mass_norm = mr
    consts = Zhong06Constants(bot, policy, structure, mass_norm, dt)
    (;nq, nλ, nu, nc, ns) = consts
    T = get_numbertype(structure)
    
    jacobian_workspace = Zhong06JacobianWorkspace(bot)
    
    Zhong06_CCP_Constant_Mass_Inner_Cache(
        bot,
        policy,
        field,
        env,
        jacobian_workspace,
        consts,
        options,
    )
end

function solve!(sim::Simulator,solver_cache::Zhong06_CCP_Constant_Mass_Inner_Cache;
                dt,
                ftol=1e-14,
                xtol=ftol,
                maxiters=50,
                max_restart=3,
                exception=true, 
                show_meter = true,
            )
    (;prob,controller,tspan,restart,totalstep) = sim
    (;bot,env) = prob
    (;structure,traj,contacts_traj,contact_caches_traj) = bot
    (;options,consts,jacobian_workspace) = solver_cache
    (;activated_bits) = structure.contacts_related
    (;Mₘ, M⁻¹ₘ ) = jacobian_workspace
    (;nq, nλ, ns, mass_norm) = consts
    h = dt
    
    q0 = traj.q[begin]
    λ0 = traj.λ[begin]
    s0 = traj.s[begin]
    q̇0 = traj.q̇[begin]
    
    state_mid = deepcopy(traj[begin])

    activate_frictional_contacts!(structure,env,q0)
    
    T = eltype(q0)
    
    nx = nq + nλ + ns
    
    # Create Newton workspace
    # We pass nΛ=0 because we don't use the contact part of NewtonWorkspace for Λ
    workspace = NewtonWorkspace(T, nx, 0, 0)
    
    # Auxiliary workspaces for line search
    Res_α0 = zeros(T, nx)
    Jac_α0 = zeros(T, nx, nx)
    
    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=show_meter)
    
    for timestep = 1:totalstep
        #---------Time Step k Control-----------
        @debug "Timestep $timestep Begin" time = traj.t[timestep]
        cₖ = contacts_traj[timestep]
        cₖ₊₁ = contacts_traj[timestep+1]
        qₖ = traj.q[timestep]
        q̇ₖ = traj.q̇[timestep]
        sₖ = traj.s[timestep]
        qₖ₊₁   = traj.q[timestep+1]
        q̇ₖ₊₁   = traj.q̇[timestep+1]
        
        qˣ = qₖ .+ dt./2 .*q̇ₖ
        qₖ₊₁ .= qₖ .+ dt .*q̇ₖ
        q̇ₖ₊₁ .= q̇ₖ
        
        contact_cache = activate_frictional_contacts!(structure,env,qˣ;checkpersist=options.checkpersist)
        (;na,L) = contact_cache
        activated_bits .= contact_cache.activated_bits
        isconverged = false
        normRes = typemax(T)
        iteration_break = 0
        nΛ = 3na
        contact_workspace = InnerContactWorkspace(T, nx, nΛ)
        
        # Setup contact workspace
        (;Λ, Λʳ, ΔΛ, 𝐁, 𝐁t, 𝐛, 𝐜ᵀ, 𝐍, 𝐲) = contact_workspace
        # Reset contact variables
        Λ .= 0
        state_k = InnerContactCoordinatesState(traj[timestep],T[],)
        state_kp1 = InnerContactCoordinatesState(traj[timestep+1],Λ)
        # Note: In solver_state, qₖ corresponds (previous/known), qₖ₊₁ to (current/unknown)
        solver_state = Zhong06SolverState(
            state_k, state_kp1, 
            InnerContactCoordinatesState(state_mid,T[],), 
            h,
        )
        
        get_frictional_directions_and_positions!(structure, contact_cache, state_k, Λ)
        
       
        precompute!(jacobian_workspace,solver_state,structure)

        qˣ = qₖ .+ dt./2 .*q̇ₖ
        qₖ₊₁ .= qₖ .+ dt .*q̇ₖ
        q̇ₖ₊₁ .= q̇ₖ

        restart_count = 0
        Λ_guess = 0.1
        
        while restart_count < max_restart
            @debug "Restart $restart_count Begin" num_of_active_contacts=na
            
            if na > 0
                Λ .= repeat([Λ_guess,0,0], na)
            end
            
            workspace.x[1:nq] .= qₖ₊₁
            workspace.x[nq+1:nq+nλ] .= 0.0
            workspace.x[nq+nλ+1:nx] .= sₖ
            
            Λʳ .= Λ
            
            Nmax = 50
            
            for iteration = 1:maxiters
                @debug "Iteration $iteration Begin"
                
                # Compute contact linearization if contacts active
                if na != 0
                    compute_inner_contact_linearization!(solver_cache, solver_state, contact_cache, contact_workspace)
                end
                
                # Compute residual and Jacobian
                compute_inner_residual!(workspace, solver_state, solver_cache, contact_cache, contact_workspace)
                compute_inner_jacobian!(workspace, solver_state, solver_cache, contact_cache)
                
                luJac = lu(workspace.Jac)
                
                if na == 0
                    normRes = norm(workspace.Res)
                    if normRes < ftol
                        isconverged = true
                        iteration_break = iteration-1
                        @debug "Iteration $iteration Break" normRes
                        break
                    end
                    workspace.Δx .= luJac \ (-workspace.Res)
                    workspace.x .+= workspace.Δx
                else # na != 0
                    normRes = norm(workspace.Res)
                    if normRes < ftol
                        isconverged = true
                        iteration_break = iteration-1
                        @debug "Iteration $iteration Break" normRes
                        break
                    elseif normRes > 1e10
                        iteration_break = iteration-1
                        isconverged = false
                        break
                    elseif iteration == maxiters
                        iteration_break = iteration-1
                        isconverged = false
                    end
                   
                    populate!(solver_state, workspace.x, structure, jacobian_workspace, consts)
                    get_distribution_law!(structure, contact_cache, solver_state.state_kp1)
                    sd = 1/norm(workspace.Res)^2 * I
                    ϕ0 = 0.5
                    dϕ0 = -1.0
                    c1 = 0.2
                    α0 = 2.0
                    
                    𝐍 .= 𝐜ᵀ * (luJac \ 𝐁)
                    
                    for line_search_step = 1:5
                        α0 /= 2
                        
                        # 𝐡 = 𝐲 .-𝐜ᵀ*(luJac\(Res + (1/α0) .*𝐁*Λʳₖ))
                        temp_vec = workspace.Res + (1/α0) .* 𝐁 * Λʳ
                        𝐡 = 𝐲 .- 𝐜ᵀ * (luJac \ temp_vec)
                        
                        yₖini = 1/α0 .* 𝐍 * Λʳ + 𝐡
                        Λₖini = copy(Λʳ)
                        Λₖini[begin+1:3:end] .= 0.0
                        Λₖini[begin+2:3:end] .= 0.0
                        yₖini .= abs.(yₖini)
                        yₖini[begin+1:3:end] .= 0.0
                        yₖini[begin+2:3:end] .= 0.0
                        
                        IPM!(Λ, na, nΛ, Λₖini, yₖini, (1/α0).*𝐍 .+ L, 𝐡; ftol, Nmax)
                        
                        ΔΛ .= (Λ - Λʳ)
                        minusResΛ = -workspace.Res + 𝐁 * ((1/α0) .* ΔΛ)
                        workspace.Δx .= luJac \ minusResΛ
                        
                        # Check line search condition
                        # We need to compute residual at new state
                        # Use temp workspace for check
                        workspace_α0 = NewtonWorkspace(workspace.x .+ workspace.Δx, Res_α0, Jac_α0, zeros(T, nx), zeros(T, nx), zeros(T, 0), zeros(T, 0, 0))
                        
                        # Update solver state for check
                        solver_state.qₖ₊₁ .= workspace_α0.x[1:nq]
                        solver_state.λₘ   .= workspace_α0.x[nq+1:nq+nλ]
                        solver_state.sₖ₊₁ .= workspace_α0.x[nq+nλ+1:nx]
                        
                        # Compute residual (without contact linearization update)
                        # We need to temporarily set Λ to Λʳ + ΔΛ
                        Λ_temp = Λʳ + ΔΛ
                        # We can modify contact_workspace.Λ temporarily? 
                        # Or pass Λ explicitly to compute_inner_residual?
                        # compute_inner_residual uses contact_workspace.Λ.
                        # So we update it.
                        contact_workspace.Λ .= Λ_temp
                        
                        compute_inner_residual!(workspace_α0, solver_state, solver_cache, contact_cache, contact_workspace)
                        
                        ϕα0 = (transpose(Res_α0) * sd * Res_α0) / 2
                        if ϕα0 <= ϕ0 + c1 * α0 * dϕ0
                            @debug "Variables" line_search_step α0
                            @debug "Line Search Break"
                            break
                        end
                    end
                    
                    Λ .= Λʳ .+ ΔΛ
                    Λʳ .= Λ
                    workspace.x .+= workspace.Δx
                end
                @debug "Iteration $iteration End  " normRes
            end
            
            if isconverged
                break
            end
            
            restart_count += 1
            if na > 0
                Λ_guess = max(Λ_guess/10, maximum(abs.(Λ[begin:3:end])))
            end
            @debug "Restart $restart_count End  "
        end
        
        populate!(solver_state, workspace.x, structure, jacobian_workspace, consts)
        
        Dₘ = contact_cache.Dper
        Dₖ₊₁ = contact_cache.Dimp
        
        if na != 0
            update_contacts!(cₖ₊₁[activated_bits], cₖ[activated_bits], Dₘ*(qₖ₊₁.-qₖ).+Dₖ₊₁*q̇ₖ₊₁, 2mass_norm.*Λ./dt)
            contact_cache.Λ .= Λ
            contact_caches_traj[timestep+1] = contact_cache
        end
        
        if !isconverged
            @warn "Newton max iterations $maxiters, at timestep=$timestep, normRes=$(normRes), restart_count=$(restart_count), num_active_contacts=$(na)"
            if exception
                @error "Not converged!"
                break
            end
        end
        
        @debug "Timestep $timestep End  " normRes restart_count iteration_break num_of_active_contacts=na
        next!(prog)
    end
    sim, solver_cache
end
