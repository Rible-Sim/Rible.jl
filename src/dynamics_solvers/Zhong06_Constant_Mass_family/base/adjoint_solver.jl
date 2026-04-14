@concrete struct Adjoint_Sensitivity_Zhong06_Constant_Mass_Cache  <: AbstractZhong06Cache
    solver
    jacobian_workspace
    ys
    ηs
    adjoint_ws
    jac_blocks
    state_mid
    cg_term
    cg_traj
    ∂J∂x₀ᵀ
    ∂J∂θᵀ
    ∂J∂cᵀ
    ∂ϕ∂θᵀ_vec
    ∂ϕ∂xᵀₖ_vec
    ∂ϕ∂xᵀₖ₊₁_vec
end

function generate_cache(
        prob::DynamicsProblem,
        solver::AdjointDynamicsSensitivitySolver{
            <:DynamicsSolver{<:Zhong06},
            <:DiscreteAdjointDynamicsSolver{<:Zhong06}
        },
        has_constant_mass_matrix::Val{true};
        dt,totalstep,kwargs...
    )
    @info "Adjoint_Sensitivity_Zhong06_Constant_Mass_Cache"
    (;bot,policy,env) = prob
    (;structure,hub) = bot
    (;coalition) = hub
    options = merge(
        (checkpersist=true,), #default
        prob.options,
        solver.forward_solver.options,
    )
    (;objt) = solver.adjoint_solver
    T = get_numbertype(structure)
    Mₘ = assemble_M(structure)
    ys = [
        ComponentArray(
            a = zero(state.q),
            b = zero(state.p),
            μ = zero(state.λ),
            τ = zero(state.s)
        )
        for state in bot.traj[begin+1:end]
    ]
    
    ηs = get_trajectory_cost_weights(objt,bot.traj.t[begin:end-1],dt)
    consts = Zhong06Constants(bot, policy, structure, norm(Mₘ,Inf), dt)
    (;nq, nλ, nu, nc, nθ, ns) = consts
    ny = 2*nq+nλ+ns
    nx = ny
    # Create jacobian_workspace
    jacobian_workspace = Zhong06JacobianWorkspace(bot)
    
    
    ∂J∂θᵀ = [ zeros(T,nθ) for _ in 1:totalstep ]
    ∂J∂cᵀ = [ zeros(T,nc) for _ in 1:totalstep ]
    ∂J∂x₀ᵀ = zeros(T,nx)

    ∂ϕ∂θᵀ_vec = [ zeros(T,nθ) for _ in 0:totalstep-1 ]
    ∂ϕ∂xᵀₖ_vec = [ zeros(T,nx) for _ in 0:totalstep-1 ]
    ∂ϕ∂xᵀₖ₊₁_vec = [ zeros(T,nx) for _ in 0:totalstep-1 ]

    adjoint_ws = NewtonWorkspace(T, ny, 0, 0)
    
    # Create Jacobian blocks
    Jacᵏ⁺¹ₖ_backup = zeros(T,ny,ny)

    jac_blocks = Zhong06JacobianBlocks(T, nx, nu, nc, Jacᵏ⁺¹ₖ_backup)

    state_mid = deepcopy(bot.traj[begin])
        
    ∂ϕf∂qᵀ = zeros(T,nq)
    ∂ϕf∂q̇ᵀ = zeros(T,nq)
    ∂ϕf∂pᵀ = zeros(T,nq)
    ∂ϕf∂uᵀ = zeros(T,nu)
    ∂ϕf∂sᵀ = zeros(T,ns)

    cg_term = CostGradient(
        ∂ϕf∂qᵀ,∂ϕf∂q̇ᵀ,∂ϕf∂pᵀ,∂ϕf∂uᵀ, 
        ∂ϕf∂sᵀ, zeros(T,0), zeros(T,0) # ∂ϕ∂θᵀ, ∂ϕ∂cᵀ
    )

    ∂ϕ∂qᵀ = zeros(T,nq)
    ∂ϕ∂q̇ᵀ = zeros(T,nq)
    ∂ϕ∂pᵀ = zeros(T,nq)
    ∂ϕ∂uᵀ = zeros(T,nu)
    ∂ϕ∂sᵀ = zeros(T,ns)

    cg_traj = CostGradient(
        ∂ϕ∂qᵀ, ∂ϕ∂q̇ᵀ, ∂ϕ∂pᵀ, ∂ϕ∂uᵀ, 
        ∂ϕ∂sᵀ, zeros(T,0), zeros(T,0) # ∂ϕ∂θᵀ, ∂ϕ∂cᵀ
    )

    Adjoint_Sensitivity_Zhong06_Constant_Mass_Cache(
        solver,
        jacobian_workspace,
        ys,ηs,
        adjoint_ws,jac_blocks, state_mid,
        cg_term,cg_traj,
        ∂J∂x₀ᵀ,∂J∂θᵀ,∂J∂cᵀ,
        ∂ϕ∂θᵀ_vec,∂ϕ∂xᵀₖ_vec,∂ϕ∂xᵀₖ₊₁_vec
    )
end

function solve!(simulator::Simulator,
    forward_cache::Zhong06_Constant_Mass_Cache,
    solver_cache::Adjoint_Sensitivity_Zhong06_Constant_Mass_Cache;
        dt,ftol=1e-14,verbose=false,maxiters=50,
        progress=true,exception=true
    )
    @info "Solving Adjoint_Sensitivity_Zhong06_Constant_Mass_Cache"
    (;prob,controller,tspan,restart,totalstep,) = simulator
    (;bot,policy,env) = prob
    (;hub,traj,control_traj) = bot
    (;structure) = bot
    (;
        solver,
        jacobian_workspace,
        ys,ηs,
        adjoint_ws,
        jac_blocks, state_mid,
        cg_term, cg_traj,
        ∂J∂x₀ᵀ,∂J∂θᵀ,∂J∂cᵀ,
        ∂ϕ∂xᵀₖ_vec,∂ϕ∂xᵀₖ₊₁_vec
    ) = solver_cache
    (;objt) = solver.adjoint_solver
    # Extract from jacobian_workspace container
    (;Mₘ) = jacobian_workspace
    T = get_numbertype(structure)
    
    consts = Zhong06Constants(bot, policy, structure, norm(Mₘ,Inf), dt)
    (;nq, nλ, nu, nc, nθ, ns) = consts
    nx = ny = 2nq+nλ+ns
    ∂ϕf∂xᵀ = zeros(T,nx)
    ∂ϕ∂xᵀ = zeros(T,nx)
    step = 0
    h = dt

    uₖ₊₁ = uₘ = control_traj.u[totalstep+1]

    solver_state = Zhong06SolverState(
        traj[totalstep], traj[totalstep+1], state_mid,
        h
    )
    
    interpolate!(solver_state)
    
    switch!(structure, solver_state.state_k)

    Zhong06_constant_mass_Jac!(
        jac_blocks,
        solver_state,
        bot,env,policy,forward_cache,
        consts,jacobian_workspace
    )
    
    
    cost_gradient!(
        cg_term,jacobian_workspace,
        bot,policy,env,objt,
        traj[totalstep+1],uₖ₊₁;mode=:terminal
    )
    
    ∂ϕf∂xᵀ[   1: nq] .= cg_term.∂ϕ∂qᵀ #.+ ∂u∂qₖ₊₁'*∂ϕf∂uᵀ
    ∂ϕf∂xᵀ[nq+1:2nq] .= cg_term.∂ϕ∂pᵀ #.+ ∂u∂pₖ₊₁'*∂ϕf∂uᵀ
    ∂ϕf∂xᵀ[2nq+nλ+1:2nq+nλ+ns] .= cg_term.∂ϕ∂sᵀ

    cost_gradient!(
        cg_traj,jacobian_workspace,
        bot,policy,env,objt,
        traj[totalstep],control_traj.u[totalstep];mode=:trajectory
    )

    ∂ϕ∂qₖᵀ, ∂ϕ∂qₖ₊₁ᵀ, ∂ϕ∂pₖᵀ, ∂ϕ∂pₖ₊₁ᵀ, ∂ϕ∂λᵀ, ∂ϕ∂sₖᵀ, ∂ϕ∂sₖ₊₁ᵀ = vjp_wrt_state(cg_traj.∂ϕ∂uᵀ,policy,bot,nu,solver_cache,solver_state)

    ∂ϕ∂xᵀ[   1: nq] .= cg_traj.∂ϕ∂qᵀ .+ ∂ϕ∂qₖᵀ
    ∂ϕ∂xᵀ[nq+1:2nq] .= cg_traj.∂ϕ∂pᵀ .+ ∂ϕ∂pₖᵀ
    ∂ϕ∂xᵀ[2nq+nλ+1:2nq+nλ+ns] .= cg_traj.∂ϕ∂sᵀ .+ ∂ϕ∂sₖᵀ
    ∂ϕ∂xᵀₖ_vec[totalstep] .= ∂ϕ∂xᵀ
    ∂ϕ∂xᵀₖ₊₁_vec[totalstep][   1: nq] .= ∂ϕ∂qₖ₊₁ᵀ
    ∂ϕ∂xᵀₖ₊₁_vec[totalstep][nq+1:2nq] .= ∂ϕ∂pₖ₊₁ᵀ
    ∂ϕ∂xᵀₖ₊₁_vec[totalstep][2nq+nλ+1:2nq+nλ+ns] .= ∂ϕ∂sₖ₊₁ᵀ

    yN = ys[totalstep]
    yN .= transpose(jac_blocks.Jacᵏ⁺¹ₖ₊₁)\(-∂ϕf∂xᵀ-ηs[totalstep]*∂ϕ∂xᵀₖ₊₁_vec[totalstep])

    v_total = ηs[totalstep] * cg_traj.∂ϕ∂uᵀ + transpose(jac_blocks.Jacᵏ⁺¹ₘu) * yN
    ∂J∂θᵀ[totalstep] .= 0
    accumulate_param_grad!(∂J∂θᵀ[totalstep], policy, v_total, solver_state, bot)

    ∂J∂cᵀ[totalstep] = transpose(jac_blocks.Jacᵏ⁺¹ₘc)*yN# .+ ηs[totalstep]*∂ϕ∂cᵀ

    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = totalstep:-1:2
        @debug "Timestep $timestep Begin" time = traj.t[timestep]
        jac_blocks.Jacᵏ⁺¹ₖ_backup .= jac_blocks.Jacᵏ⁺¹ₖ
        yₖ₊₁ = ys[timestep]
        yₖ   = ys[timestep-1]
        
        solver_state = Zhong06SolverState(
            traj[timestep-1], traj[timestep], state_mid,
            h
        )

        interpolate!(solver_state)
        
        switch!(structure, solver_state.state_k)

        Zhong06_constant_mass_Jac!(
            jac_blocks,
            solver_state,
            bot,env,policy,forward_cache,
            consts,jacobian_workspace
        )
        
        cost_gradient!(
            cg_traj,jacobian_workspace,
            bot,policy,env,objt,
            traj[timestep-1],hub.state.u;
            mode=:trajectory
        )

        ∂ϕ∂qₖᵀ, ∂ϕ∂qₖ₊₁ᵀ, ∂ϕ∂pₖᵀ, ∂ϕ∂pₖ₊₁ᵀ, ∂ϕ∂λᵀ, ∂ϕ∂sₖᵀ, ∂ϕ∂sₖ₊₁ᵀ = vjp_wrt_state(cg_traj.∂ϕ∂uᵀ,policy,bot,nu,solver_cache,solver_state)

        ∂ϕ∂xᵀ[   1: nq] .= cg_traj.∂ϕ∂qᵀ .+ ∂ϕ∂qₖᵀ 
        ∂ϕ∂xᵀ[nq+1:2nq] .= cg_traj.∂ϕ∂pᵀ .+ ∂ϕ∂pₖᵀ 
        ∂ϕ∂xᵀ[2nq+nλ+1:2nq+nλ+ns] .= cg_traj.∂ϕ∂sᵀ .+ ∂ϕ∂sₖᵀ
   
        ∂ϕ∂xᵀₖ_vec[timestep-1] .= ∂ϕ∂xᵀ
        ∂ϕ∂xᵀₖ₊₁_vec[timestep-1][   1: nq] .= ∂ϕ∂qₖ₊₁ᵀ 
        ∂ϕ∂xᵀₖ₊₁_vec[timestep-1][nq+1:2nq] .= ∂ϕ∂pₖ₊₁ᵀ 
        ∂ϕ∂xᵀₖ₊₁_vec[timestep-1][2nq+nλ+1:2nq+nλ+ns] .= ∂ϕ∂sₖ₊₁ᵀ 

        adjoint_ws.Jac .= transpose(jac_blocks.Jacᵏ⁺¹ₖ₊₁)
        Res_adj = adjoint_ws.Res
        mul!(Res_adj, transpose(jac_blocks.Jacᵏ⁺¹ₖ_backup), yₖ₊₁)
        @. Res_adj += ηs[timestep]*∂ϕ∂xᵀₖ_vec[timestep] + ηs[timestep-1]*∂ϕ∂xᵀₖ₊₁_vec[timestep-1]
        _lu_solve_from_jacobian!(adjoint_ws)
        copyto!(yₖ, adjoint_ws.Δx)

        v_total = ηs[timestep-1] * cg_traj.∂ϕ∂uᵀ + transpose(jac_blocks.Jacᵏ⁺¹ₘu) * yₖ
        ∂J∂θᵀ[timestep-1] .= 0
        accumulate_param_grad!(∂J∂θᵀ[timestep-1], policy, v_total, solver_state, bot)

        ∂J∂cᵀ[timestep-1] = transpose(jac_blocks.Jacᵏ⁺¹ₘc)*yₖ
        step += 1
        @debug "Timestep $timestep End  "
        next!(prog)
    end
    Jac¹₀ = jac_blocks.Jacᵏ⁺¹ₖ
    ∂J∂x₀ᵀ .= transpose(Jac¹₀)*ys[1] .+ ηs[1]*∂ϕ∂xᵀₖ_vec[1]

end
