function generate_cache(
        prob::AbstractDynamicsProblem,
        solver::DynamicsSolver{
            <:RKIntegrator,
        },
        ::Val{true};
        dt, kargs...
    )
    @info "RungeKutta_Constant_Mass_Cache"
    (; bot, policy, env) = prob
    (; traj, structure) = bot
    (;
        tableau
    ) = solver.integrator
    options = merge(
        (checkpersist=true,), #default
        prob.options,
        solver.options,
    )
    (; M, M⁻¹) = build_mass_matrices(structure)
    
    q0 = traj.q[begin]
    s0 = traj.s[begin]
    λ0 = traj.λ[begin]
    q̇0 = traj.q̇[begin]
    p0 = traj.p[begin] .= M * q̇0
    T = eltype(q0)
    nq = length(q0)
    nλ = length(λ0)
    ns = length(s0)
    nu = get_num_of_actions(bot)
    F = zeros(T, nq)
    ∂F∂q = zeros(T, nq, nq)
    ∂F∂q̇ = zeros(T, nq, nq)
    ∂F∂u = zeros(T, nq, nu)
    ∂F∂s = zeros(T, nq, ns)
    ∂S∂q = zeros(T, ns, nq)
    ∂S∂s = zeros(T, ns, ns)
    A = zeros(T,nλ,nq)
    ∂Aᵀλ∂q = zeros(T, nq, nq)
    nx = nq + nλ + ns
    ny = 2nq
    qpₖ₊₁ = vcat(q0, p0)
    qpₖ = deepcopy(qpₖ₊₁)
    qp_rk_stages = [
        deepcopy(qpₖ₊₁)
        for _ in 1:tableau.s
    ]
    x = vcat(
        reduce(vcat, qp_rk_stages),
        λ0,
    )
    Res = zero(x)
    Jac = Res * transpose(Res)
    JK = zero(qpₖ₊₁) * transpose(qpₖ₊₁)

    consts = RungeKuttaConstants(nq, nλ, ns, nx, ny, dt, tableau)
    workspace = RungeKuttaJacobianWorkspace(
        M, M⁻¹, F, ∂F∂q, ∂F∂q̇, ∂F∂u, ∂F∂s, ∂S∂q, ∂S∂s,A,∂Aᵀλ∂q,
        x, qpₖ₊₁, qpₖ, qp_rk_stages, Res, Jac, JK
    )

    RungeKutta_Constant_Mass_Cache(solver, consts, workspace, options)
end

function interpolate!(workspace, consts, x)
    (; qpₖ₊₁, qpₖ, qp_rk_stages) = workspace
    (; tableau, ny, h) = consts
    (; a, b) = tableau
    nstage = tableau.s
    
    qpₖ₊₁ .= qpₖ
    for i = 1:nstage
        qp_rk_stages[i] .= qpₖ
    end

    for j = 1:nstage
        kⱼ = @view x[(ny*(j-1)+1):(ny*j)]
        for i = 1:nstage
            qp_rk_stages[i] .+= h * a[i, j] * kⱼ
        end
        qpₖ₊₁ .+= h * b[j] * kⱼ
    end
end

function guess_newton!(x, solver_state, consts, jacobian_workspace)
    (; state_k) = solver_state
    (; F) = jacobian_workspace
    (; ny, nq) = consts
    
    x .= 0.0
    for i = 1:consts.tableau.s
        x[(ny*(i-1)+   1):(ny*(i-1)+ nq)] .= state_k.q̇
        x[(ny*(i-1)+nq+1):(ny*(i-1)+2nq)] .= F
    end
end

function populate!(solver_state, x, structure, workspace, consts)
    (; state_kp1) = solver_state
    (; qpₖ₊₁, qpₖ, M⁻¹) = workspace
    (; nq, ny) = consts

    interpolate!(workspace, consts, x)
        
    state_kp1.q .= qpₖ₊₁[1:nq]
    state_kp1.p .= qpₖ₊₁[nq+1:ny]
    state_kp1.q̇ .= M⁻¹ * state_kp1.p
    state_kp1.t = solver_state.state_k.t + consts.h
    
    qpₖ .= qpₖ₊₁
end

function compute_rungekutta_residual!(Res, x, workspace, consts, solver_state, bot, env, policy)
    (; M, M⁻¹, A, qpₖ₊₁, qp_rk_stages) = workspace
    (; nq, nλ, ny, h, tableau) = consts
    (; tableau) = consts
    (; c) = tableau
    (; state_k, inst_state) = solver_state
    tₖ = state_k.t
    nstage = tableau.s
    
    interpolate!(workspace, consts, x)
    
    λₖ₊₁ = @view x[(nstage*ny+1):(nstage*ny+nλ)]
    
    # Internal stages
    for i = 1:nstage
        kᵢ = @view x[(ny*(i-1)+1):(ny*i)]
        
        y = qp_rk_stages[i]
        t = tₖ + c[i] * h
        
        inst_state.q .= y[1:nq]
        p = @view y[nq+1:2nq]
        inst_state.q̇ .= M⁻¹*p
        inst_state.F .= 0
        gen_force!(inst_state, bot, env.field, policy)
        F = inst_state.F
        cstr_jacobian!(A, bot.structure, inst_state)
        
        Res[(ny*(i-1)+1):(ny*(i-1)+nq)] .= kᵢ[1:nq] .- inst_state.q̇
        Res[(ny*(i-1)+nq+1):(ny*i)] .= kᵢ[nq+1:end] .- (F .- transpose(A) * λₖ₊₁)
    end
    
    # Constraint residual
    qₖ₊₁ = @view qpₖ₊₁[1:nq]
    inst_state.q .= qₖ₊₁
    ϕ = @view Res[(nstage*ny+1):(nstage*ny+nλ)] 
    #todo proper scaling 
    cstr_function!(ϕ, bot.structure, inst_state)
end

function compute_rungekutta_jacobian!(Jac, x, workspace, consts, solver_state, bot, env, policy)
    (;field) = env
    (; M, M⁻¹, ∂F∂q, ∂F∂q̇, ∂F∂u, ∂F∂s, A, ∂Aᵀλ∂q, qpₖ₊₁, qp_rk_stages, JK) = workspace
    (; nq, nλ, ny, h, tableau) = consts
    (; tableau) = consts
    (; a, b, c) = tableau
    (; state_k, inst_state) = solver_state
    tₖ = state_k.t
    nstage = tableau.s

    interpolate!(workspace, consts, x)
    
    Jac .= 0.0
    
    for i = 1:nstage
        is = ny*(i-1)
        for j = 1:ny
            Jac[is+j,is+j] = 1
        end
    end
    
    λₖ₊₁ = @view x[(nstage*ny+1):(nstage*ny+nλ)]
    inst_state.λ .= λₖ₊₁
    for j  = 1:nstage
        jsta = ny*(j-1) + 1
        jend = ny * j
        for i = 1:nstage
            ista = ny*(i-1)+1
            iend = ny * i
            
            ỹᵢ = qp_rk_stages[i]
            q̃ᵢ = @view ỹᵢ[1:nq]
            t = tₖ + c[i] * h
            #todo
            inst_state.t = t
            inst_state.q .= q̃ᵢ
            inst_state.q̇ .= zero(q̃ᵢ)
            
            ∂F∂q .= 0
            ∂F∂q̇ .= 0
            ∂F∂u .= 0
            gen_force_state_jacobian!(∂F∂q, ∂F∂q̇, ∂F∂u, bot, field, policy, inst_state, ∂F∂s)
            # q̃ᵢ, λₖ₊₁

            cstr_forces_jacobian!(∂Aᵀλ∂q, bot.structure, inst_state)
            
            JK .= 0
            JK[1:nq,(nq+1):ny] .= M⁻¹
            JK[(nq+1):ny, 1:nq] .= ∂F∂q .+ ∂Aᵀλ∂q
            
            Jac[ista:iend, jsta:jend] .-= h * a[i, j] * JK
        end
        
        qₖ₊₁ = @view qpₖ₊₁[1:nq]
        inst_state.q .= qₖ₊₁
        cstr_jacobian!(A, bot.structure, inst_state)
        A_qₖ₊₁ = A
        Jac[(nstage*ny+1):(nstage*ny+nλ), jsta:(ny*(j-1)+nq)] .= h * b[j] * A_qₖ₊₁
    end
    
    for i = 1:nstage
        ỹᵢ = qp_rk_stages[i]
        q̃ᵢ = @view ỹᵢ[1:nq]
        inst_state.q .= q̃ᵢ
        cstr_jacobian!(A, bot.structure, inst_state)
        Ãᵢ = A
        Jac[(ny*(i-1)+1+nq):(ny*(i-1)+ny), (nstage*ny+1):(nstage*ny+nλ)] .= transpose(Ãᵢ)
    end
    
    Jac[(nstage*ny+1):(nstage*ny+nλ),(nstage*ny+1):(nstage*ny+nλ)] .= 0.0
end

function solve!(sim::Simulator, cache::RungeKutta_Constant_Mass_Cache;
    dt, ftol=1e-14, verbose=false, maxiters=50,
    progress=true, exception=true)
    (; prob, controller, tspan, restart, totalstep) = sim
    (; bot,) = prob
    (; structure, traj) = bot
    (; env, policy) = prob
    
    (;
        jacobian_workspace,
        consts
    ) = cache
    (; x, Res, Jac, ) = jacobian_workspace
    (; h) = consts
    
    T = get_numbertype(bot)

    iteration_break = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    dg_step = ceil(Int, log10(totalstep)) + 1
    dg_dt = max(1, -floor(Int, log10(h)))
    wd_t = ceil(Int, log10(traj.t[end])) + dg_dt + 1 + 1
    progfmt = Printf.Format("Progress: %5.1f%%, step: %$(dg_step)u, time: %$(wd_t).$(dg_dt)f, iterations: %s \n")

    inst_state = deepcopy(traj[begin])
    for timestep = 1:totalstep
        #---------Step k Control-----------

        solver_state = RungeKuttaSolverState(traj[timestep], traj[timestep+1], inst_state, h)
        guess_newton!(x, solver_state, consts, jacobian_workspace)
        
        isconverged = false
        normRes = typemax(T)
        
        for iteration = 1:maxiters

            compute_rungekutta_residual!(Res, x, 
                jacobian_workspace, consts, solver_state, 
                bot, env, policy
            )
            
            normRes = norm(Res)
            iteration_break = iteration - 1
            if normRes < ftol
                isconverged = true
                break
            end

            compute_rungekutta_jacobian!(Jac, x, 
                jacobian_workspace, consts, solver_state, 
                bot, env, policy
            )
            
            x .+= Jac \ (-Res)

        end
        
        if !isconverged
            @warn "Newton max iterations $maxiters, Step=$timestep, normRes=$normRes"
            if exception
                error("Not Converged!")
            end
        end
        
        # Final update of qpₖ₊₁
        populate!(solver_state, x, structure, jacobian_workspace, consts)
        
        #---------Step k finisher-----------
        #---------Step k finisher-----------
        if verbose
            progstr = Printf.format(
                progfmt,
                floor(timestep / totalstep * 100; digits=1), timestep, traj.t[timestep+1], iteration_break
            )
            print(progstr)
        end
        next!(prog)
    end

    bot
end
