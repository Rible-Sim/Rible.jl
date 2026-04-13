
struct Adjoint_Sensitivity_Zhong06_CCP_Constant_Mass_Cache_Mono{T, SolverType, WorkspaceType, YType, EtaType} <: AbstractZhong06Cache
    solver::SolverType
    jacobian_workspace::WorkspaceType
    ys::Vector{YType}
    ηs::Vector{EtaType}
    nys::Vector{Int}
    nxs::Vector{Int}
    ∂J∂x₀ᵀ::Vector{T}
    ∂J∂θᵀ::Vector{Vector{T}}
    ∂J∂cᵀ::Vector{Vector{T}}
    ∂ϕ∂θᵀ_vec::Vector{Vector{T}}
    ∂ϕ∂xᵀₖ_vec::Vector{Vector{T}}
    ∂ϕ∂xᵀₖ₊₁_vec::Vector{Vector{T}}
end

function generate_cache(
        prob::DynamicsProblem,
        solver::AdjointDynamicsSensitivitySolver{
            <:DynamicsSolver{<:Zhong06},
            <:DiscreteAdjointDynamicsSolver{
                Zhong06,
                <:AbstractBodySolver,
                <:AbstractApparatusSolver,
                <:MonolithicContactSolver
            }
        },
        has_constant_mass_matrix::Val{true};
        dt,totalstep,kwargs...
    )
    
    @info "Adjoint_Sensitivity_Zhong06_CCP_Constant_Mass_Cache_Mono"
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
    ηs = get_trajectory_cost_weights(solver.adjoint_solver.objt, bot.traj.t[begin:end-1],dt)
    consts = Zhong06Constants(bot, policy, structure, norm(Mₘ,Inf), dt)
    (;nq, nq̌, nλ, nu, nc, nθ, n1, n2, n3) = consts
    # Create Zhong06JacobianWorkspace for caching jacobian_workspace matrices
    ns = consts.ns
    jacobian_workspace = Zhong06JacobianWorkspace(bot)
    
    ∂J∂θᵀ = [
        zeros(T,nθ)
        for _ in  1:totalstep
    ]
    ∂J∂cᵀ = [
        zeros(T,nc) 
        for _ in 1:totalstep
    ]

    ∂J∂x₀ᵀ = zeros(T,2nq̌+nλ)
    
    ∂ϕ∂θᵀ_vec = [ zeros(T,nθ) for _ in 0:totalstep-1 ]
    
    ys = [
        ComponentArray(
            adj_q = zero(state.q),
            adj_p = zero(state.p),
            adj_λ = zero(state.λ),
            adj_Λ = zeros(T,3*cache.na),
            adj_Γ = zeros(T,3*cache.na)
        )
        for (state, cache) in zip(
                            traj[begin+1:end],
            contact_caches_traj[begin+1:end]
        )
    ]
    
    nys = length.(ys)

    nxs = [
        2nq̌+nλ+2*3contact_caches_traj[timestep].na
        for timestep = 1:totalstep+1
    ]
    
    ∂ϕ∂xᵀₖ_vec = [ 
        zeros(T,nxs[timestep]) for timestep in 1:totalstep
    ]
    
    ∂ϕ∂xᵀₖ₊₁_vec = [ 
        zeros(T,nxs[timestep+1]) for timestep in 1:totalstep 
    ]

    Adjoint_Sensitivity_Zhong06_CCP_Constant_Mass_Cache_Mono(
        solver,
        jacobian_workspace,
        ys, ηs,
        nys, nxs,
        ∂J∂x₀ᵀ, ∂J∂θᵀ, ∂J∂cᵀ,
        ∂ϕ∂θᵀ_vec, ∂ϕ∂xᵀₖ_vec, ∂ϕ∂xᵀₖ₊₁_vec
    )
end

function solve!(simulator::Simulator,
        forward_cache,
        solver_cache::Adjoint_Sensitivity_Zhong06_CCP_Constant_Mass_Cache_Mono;
        dt,ftol=1e-14,verbose=false,maxiters=50,
        progress=true,exception=true
    )
    @info "Solving Adjoint_Sensitivity_Zhong06_CCP_Constant_Mass_Cache_Mono"
    (;prob,controller,tspan,restart,totalstep) = simulator
    (;bot,env,policy) = prob
    (;hub,structure,traj,contacts_traj,control_traj,contact_caches_traj) = bot
    (;
        solver,
        jacobian_workspace,
        ys,ηs,
        nys,nxs,
        ∂J∂x₀ᵀ,∂J∂θᵀ,∂J∂cᵀ,
        ∂ϕ∂θᵀ_vec,∂ϕ∂xᵀₖ_vec,∂ϕ∂xᵀₖ₊₁_vec
    ) = solver_cache
    
    # Unpack from jacobian_workspace
    (;Mₘ, M⁻¹ₘ, ∂Mₘhq̇ₘ∂qₘ, ∂Fₘ∂u, ∂Fₘ∂c) = jacobian_workspace
    T = get_numbertype(structure)
    consts = Zhong06Constants(bot, policy, structure, norm(Mₘ,Inf), dt)
    (;nq, nq̌, nλ, nu, nc, nθ, n1, n2, n3) = consts
    step = 0
    # Create context bundle (includes all data needed for closures)
    (;objt) = solver.adjoint_solver
    
    
    #---- Adjoint 
    # terminal step
    nxₖ₊₁ = nxs[totalstep+1]
    nxₖ   = nxs[totalstep]
    ∂S∂xᵀ = zeros(T,nxₖ₊₁)

    Jacᵏ⁺¹ₖ_backup = zeros(T,nxₖ₊₁,nxₖ)
    
    jacobians = Zhong06JacobianBlocks(T, nxₖ₊₁, nu, nc, Jacᵏ⁺¹ₖ_backup)
    
    yN = ys[totalstep]
    compute_adjoint_step!(
        yN, nothing,  # yₖ, yₖ₊₁ (no yₖ₊₁ for terminal)
        ∂J∂θᵀ[totalstep], ∂J∂cᵀ[totalstep],
        ∂ϕ∂θᵀ_vec[totalstep],
        ∂ϕ∂xᵀₖ_vec[totalstep], ∂ϕ∂xᵀₖ₊₁_vec[totalstep],
        totalstep+1, #timestep
        totalstep, #totalstep
        traj, control_traj, contact_caches_traj,
        jacobians,
        jacobian_workspace,
        consts,
        bot, policy, forward_cache, solver, env, objt, ηs, nxs;
        is_terminal=true,
        ∂S∂xᵀ
    )
    Jacᵏ⁺¹ₖ_backup[:,1:(2nq+nλ)] = jacobians.Jacᵏ⁺¹ₖ[:,1:(2nq+nλ)]

    # intermiate steps
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    @debug "Variables" yN, ∂S∂xᵀ
    for timestep = totalstep:-1:2
        #---------Step k Control-----------
        @debug "Timestep $timestep Begin" time = traj.t[timestep]
        #---------Step k Control---------
                
        # Allocate matrices for this step
        nxₖ₊₁ = nxs[timestep]
        nxₖ   = nxs[timestep-1]

        yₖ₊₁ = ys[timestep  ]
        yₖ   = ys[timestep-1]
        
        @debug "Variables" nxₖ₊₁ nxₖ size(yₖ) size(yₖ₊₁)
        @debug "Variables" size(Jacᵏ⁺¹ₖ_backup)
        @debug "Variables" yₖ norm(yₖ)
        
        # Update jacobians structure with new matrices
        jacobians = Zhong06JacobianBlocks(T, nxₖ₊₁, nu, nc, Jacᵏ⁺¹ₖ_backup)
        
        # Solve intermediate adjoint step
        compute_adjoint_step!(
            yₖ, yₖ₊₁,
            ∂J∂θᵀ[timestep-1], ∂J∂cᵀ[timestep-1],
            ∂ϕ∂θᵀ_vec[timestep-1],
            ∂ϕ∂xᵀₖ_vec[timestep], ∂ϕ∂xᵀₖ₊₁_vec[timestep-1],
            timestep, totalstep,
            traj, control_traj, contact_caches_traj,
            jacobians,
            jacobian_workspace,
            consts,
            bot, policy, forward_cache, solver, env, objt, ηs, nxs;
            is_terminal=false
        )
        
        # Backup Jacobian for next iteration
        Jacᵏ⁺¹ₖ_backup = zeros(T,nxₖ₊₁,nxₖ)
        Jacᵏ⁺¹ₖ_backup[:,1:(2nq+nλ)] = jacobians.Jacᵏ⁺¹ₖ[:,1:(2nq+nλ)]

        #---------Step k finisher-----------
        step += 1
        #---------Step k finisher-----------
        @debug "Timestep $timestep End  " num_of_active_contacts=contact_caches_traj[timestep  ].na
        next!(prog)
    end
    @debug "Variables" size(jacobians.Jacᵏ⁺¹ₖ) 
    Jac¹₀ = jacobians.Jacᵏ⁺¹ₖ
    ∂J∂x₀ᵀ .= (transpose(Jac¹₀)*ys[1])[1:(2nq+nλ)].+ ηs[1]*∂ϕ∂xᵀₖ_vec[1]
 
    @debug "Variables" ys
end


"""
Compute adjoint step: Jacobians, cost gradients, VJP, and solve for adjoint variables.
Returns updated adjoint variable yₖ and parameter gradients.
"""
function compute_adjoint_step!(
        yₖ, yₖ₊₁,  # adjoint variables (yₖ is modified in-place)
        ∂J∂θᵀ_k, ∂J∂cᵀ_k,  # parameter gradients (modified in-place)
        ∂ϕ∂θᵀ_k,  # cost gradient w.r.t. params (modified in-place)
        ∂ϕ∂xᵀₖ, ∂ϕ∂xᵀₖ₊₁,  # cost gradients w.r.t. state
        timestep, totalstep,
        traj, control_traj, contact_caches_traj,
        jacobians::Zhong06JacobianBlocks,
        jacobian_workspace::Zhong06JacobianWorkspace,
        consts::Zhong06Constants,
        bot::Robot,
        policy::AbstractPolicy,
        forward_cache,
        solver,
        env,
        objt::AbstractObjective,
        ηs::Vector,
        nxs::Vector{Int};
        is_terminal=false,
        ∂S∂xᵀ=nothing
    )
    
    # Unpack structures for readability
    (;Jacᵏ⁺¹ₖ₊₁, Jacᵏ⁺¹ₖ, Jacᵏ⁺¹ₖ_backup, Jacᵏ⁺¹ₘu, Jacᵏ⁺¹ₘc) = jacobians
    (;Fₘ, ∂F∂q, ∂F∂q̇, ∂C∂qₖ₊₁, ∂C∂pₖ₊₁, ∂Fₘ∂u, ∂Fₘ∂c, ∂ϕ∂qᵀ, ∂ϕ∂q̇ᵀ, ∂ϕ∂pᵀ, ∂ϕ∂uᵀ, ∂ϕ∂sᵀ, ∂ϕf∂qᵀ, ∂ϕf∂q̇ᵀ, ∂ϕf∂pᵀ, ∂ϕf∂uᵀ, ∂ϕf∂sᵀ) = jacobian_workspace

    (;h, mass_norm, nq, nλ, nu, nc, nθ, n1, n2, n3) = consts
    (;structure) = bot

    # Control
    uₖ = control_traj.u[timestep-1]

    # Get contact forces
    Λₖ₊₁ = contact_caches_traj[timestep].Λ
    Γₖ₊₁ = contact_caches_traj[timestep].Γ
    
    # Extract trajectory data
    state_kp1 = MonoContactCoordinatesState(traj[timestep],Λₖ₊₁,Γₖ₊₁)
    state_k = MonoContactCoordinatesState(traj[timestep-1],Λₖ₊₁,Γₖ₊₁)
    state_mid = MonoContactCoordinatesState(
        deepcopy(traj[timestep-1]),Λₖ₊₁,Γₖ₊₁
    )
    solver_state = Zhong06SolverState(state_k,state_kp1,state_mid,h)
    # Midpoint quantities
    interpolate!(solver_state)

    # Constraint Jacobians (store in jacobian_workspace)
    cstr_jacobian!(jacobian_workspace.Aₖ, structure,state_k)
    cstr_jacobian!(jacobian_workspace.Aₖ₊₁, structure,state_kp1)
    
    # Call shared Jacobian computation
    compute_zhong06_jacobian_blocks!(
        jacobians,
        jacobian_workspace,
        solver_state,
        consts,
        contact_caches_traj[timestep],
        bot, policy, env.field, forward_cache
    )
    T = get_numbertype(bot)
    # Compute cost gradients
    cost_gradient!(
        CostGradient(
            ∂ϕ∂qᵀ, ∂ϕ∂q̇ᵀ, ∂ϕ∂pᵀ, ∂ϕ∂uᵀ, 
            ∂ϕ∂sᵀ,zeros(T,0), zeros(T,0) # ∂ϕ∂θᵀ, ∂ϕ∂cᵀ
        ),jacobian_workspace,
        bot, policy, env, objt,
        state_k, uₖ,;
        mode=:trajectory
    )
    
    # Compute VJP w.r.t. state
    ∂ϕ∂qₖᵀ, ∂ϕ∂qₖ₊₁ᵀ, ∂ϕ∂pₖᵀ, ∂ϕ∂pₖ₊₁ᵀ, ∂ϕ∂λᵀ, ∂ϕ∂sₖᵀ, ∂ϕ∂sₖ₊₁ᵀ = 
        vjp_wrt_state(∂ϕ∂uᵀ, policy, bot, nu, forward_cache, solver_state)
    
    # Assemble cost gradient vectors
    ∂ϕ∂xᵀₖ[   1: nq] .= ∂ϕ∂qᵀ .+ ∂ϕ∂qₖᵀ
    ∂ϕ∂xᵀₖ[nq+1:2nq] .= ∂ϕ∂pᵀ .+ ∂ϕ∂pₖᵀ
    
    ∂ϕ∂xᵀₖ₊₁[   1: nq] .= ∂ϕ∂qₖ₊₁ᵀ
    ∂ϕ∂xᵀₖ₊₁[nq+1:2nq] .= ∂ϕ∂pₖ₊₁ᵀ
    
    # Store cost gradient w.r.t. parameters is avoided because we use accumulate_param_grad!
    
    # Solve for adjoint variable
    if is_terminal
        cost_gradient!(
            CostGradient(
                ∂ϕf∂qᵀ,∂ϕf∂q̇ᵀ,∂ϕf∂pᵀ,∂ϕf∂uᵀ,
                ∂ϕf∂sᵀ,zeros(T,0), zeros(T,0) # ∂S∂θᵀ, ∂S∂cᵀ
            ),jacobian_workspace,
            bot, policy, env, objt,
            state_kp1, uₖ,; mode=:terminal
        )
        
        ∂S∂xᵀ[   1: nq] .= ∂ϕf∂qᵀ 
        ∂S∂xᵀ[nq+1:2nq] .= ∂ϕf∂pᵀ 

        # Terminal step: yₖ = -inv(Jacᵏ⁺¹ₖ₊₁ᵀ) * (∂S∂xᵀ + η*∂ϕ∂xᵀ)
        yₖ .= transpose(Jacᵏ⁺¹ₖ₊₁)\(-∂S∂xᵀ - ηs[timestep-1]*∂ϕ∂xᵀₖ₊₁)
    else
        # Internal step: yₖ = -inv(Jacᵏ⁺¹ₖ₊₁ᵀ) * (η*∂ϕ∂xᵀₖ + η*∂ϕ∂xᵀₖ₊₁ + Jacᵏ⁺¹ₖᵀ*yₖ₊₁)
        yₖ .= transpose(Jacᵏ⁺¹ₖ₊₁)\(
            -ηs[timestep  ]*∂ϕ∂xᵀₖ
            -ηs[timestep-1]*∂ϕ∂xᵀₖ₊₁
            -transpose(Jacᵏ⁺¹ₖ_backup)*yₖ₊₁
        )
    end
    
    # Compute parameter gradients via VJP
    v_total = transpose(Jacᵏ⁺¹ₘu)*yₖ .+ ηs[timestep-1]*∂ϕ∂uᵀ
    ∂J∂θᵀ_k .= 0
    accumulate_param_grad!(∂J∂θᵀ_k, policy, v_total, solver_state, bot)
    ∂J∂cᵀ_k .= transpose(Jacᵏ⁺¹ₘc)*yₖ
    
    return nothing
end
