function generate_cache(
        prob::DynamicsProblem,
        solver::DynamicsSolver{
            <:Zhong06,
        },
        ::Val{true};
        dt,kargs...
    )
    @info "Zhong06_Constant_Mass_Cache"
    (;bot,policy,env) = prob
    (;traj,structure,hub) = bot
    options = merge(
        (checkpersist=true,), #default
        prob.options,
        solver.options,
    )
    @debug "generate_cache"  policy.ps
    (;M,M⁻¹,M̌,M̌⁻¹,Ḿ,M̄)= build_mass_matrices(structure)
    q0 = traj.q[begin]
    q̇0 = traj.q̇[begin]
    s0 = traj.s[begin]
    p0 = traj.p[begin] .= M*q̇0 # assign initial momentum to make sure it is consistent with initial velocity
    #todo partial p0
    T = get_numbertype(structure)
    
    mr = norm(Ḿ,Inf)
    mass_norm = mr
     
    consts = Zhong06Constants(bot, policy, structure, mass_norm, dt)
    (;nq, nq̌, nλ, nu, nc, nθ, n1, n2, n3, ns) = consts
    ∂Mₘhq̇ₘ∂qₘ = spzeros(T, nq, nq) # Constant mass -> derivative is zero
    
    # We need M⁻¹ₘ. M̌⁻¹ seems to be it.
    
    jacobian_workspace = Zhong06JacobianWorkspace(bot)

    
    state_mid = deepcopy(traj[begin])

    Zhong06_Constant_Mass_Cache(
        solver,
        jacobian_workspace,
        consts,
        options,
        state_mid
    )
end


function solve!(sim::Simulator,solver_cache::Zhong06_Constant_Mass_Cache;
        dt,ftol=1e-14,verbose=false,maxiters=50,
        progress=true,exception=true
    )
    @info "Solving Zhong06_Constant_Mass_Cache"
    (;prob,controller,tspan,restart,totalstep,) = sim
    (;bot,policy,env) = prob
    (;hub,structure,traj,control_traj) = bot
    (;
        jacobian_workspace,
        consts,
        options,
        state_mid
    ) = solver_cache
    #todo use dispatch to select linear solver
    # linear_solver! signature: (ws, solver_state, solver_cache, bot, policy, env)
    linear_solver! = get(options, :linear_solver!, default_linear_solver!) # set to finite_diff_linear_solver! to approximate the Jacobian
    # linear_solver! = finite_diff_linear_solver!
    (;h, mass_norm, nq̌, nλ, ns) = consts
    nx = nq̌ + nλ + ns
    
    # Unpack workspace fields for easy access
    (;M, Ḿ, M̌, M̄, M⁻¹ₘ, M̌⁻¹) = jacobian_workspace
    M⁻¹ = M⁻¹ₘ
    T = get_numbertype(bot)
    newton_workspace = NewtonWorkspace(T, nx, 0, 0)
    ls_workspace = NewtonWorkspace(
        newton_workspace.xₖ, newton_workspace.Res, newton_workspace.Jac,
        newton_workspace.Δx, newton_workspace.x, newton_workspace.𝐰,
        newton_workspace.∂Γ∂x
    )

    iteration_break = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        @debug "Timestep $timestep Begin" time = traj.t[timestep]
        #---------Step k Control-----------
        uₖ = control_traj.u[timestep]
        x = newton_workspace.x
        Res = newton_workspace.Res
        Jac = newton_workspace.Jac
        Δx = newton_workspace.Δx
      
        solver_state = Zhong06SolverState(
            traj[timestep], traj[timestep+1], state_mid,
            h
        )
  
        solver_state.state_mid.t = (traj.t[timestep] + traj.t[timestep+1])/2
        guess_newton!(x, solver_state, consts, bot)
        precompute!(jacobian_workspace, solver_state, structure)
        
        switch!(structure, solver_state.state_k)
    
        isconverged = false
        normRes = typemax(T)
       
        # Newton iteration
        for iteration = 1:maxiters
            @debug "Iteration $iteration Begin"
            compute_constant_mass_residual!(Res, x, solver_state, solver_cache, bot, policy, env)
            normRes = norm(Res)

            iteration_break = iteration - 1
            if normRes < ftol
                isconverged = true
                @debug "Iteration $iteration Break" normRes
                break
            end
            # @show timestep, iteration, norm(Res)
            linear_solver!(
                newton_workspace,
                solver_state, solver_cache, bot, policy, env
            )
            
            ϕ, dϕ, ϕdϕ = create_line_search_functions(
                newton_workspace, ls_workspace, solver_state, solver_cache, bot, policy, env)

            fx = 0.5 * normRes^2
            mul!(newton_workspace.JacΔx, Jac, Δx)
            dϕ_0 = dot(Res, newton_workspace.JacΔx)
            #users should most likely use MoreThuente, HagerZhang or BackTracking
            α = 1.0
            try
                α, fx = LS.BackTracking(order=3)(ϕ, dϕ, ϕdϕ, 1.0, fx, dϕ_0)
                @debug "Line Search" α, fx
            catch e
                if e isa LS.LineSearchException
                    @info "Line Search Failed" e
                    break
                end
                rethrow(e)
            end
            @. x += α * Δx
            @debug "Iteration $iteration End  " normRes
            
        end
        
        if !isconverged
            @warn "Newton max iterations $maxiters, Step=$timestep, normRes=$normRes"
            if exception
                error("Not Converged!")
            else
                sim.convergence = false
                break
            end
        end
        
        populate!(solver_state, x, structure, jacobian_workspace, consts)

        control!(bot, policy, solver_cache, solver_state)
        
        uₖ .= hub.state.u
        #---------Step k finisher-----------
        
        #---------Step k finisher-----------
        @debug "Timestep $timestep End  "
        next!(prog)
    end
    bot
end
