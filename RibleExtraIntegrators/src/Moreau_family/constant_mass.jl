struct Moreau_Constant_Mass_Cache{solT,cacheType}
    solver::solT
    cache::cacheType
end

struct MoreauConstants{T}
    h::T
    mass_norm::T
    nq::Int
    nλ::Int
    nx::Int
end

struct MoreauJacobianWorkspace{MT}
    ∂F∂q::MT
    ∂F∂q̇::MT
    ∂F∂u::MT
    ∂F∂s::MT
    A_buffer::MT
    Aq̇_buffer::MT
    ∂Aᵀλ_buffer::MT
end

mutable struct MoreauSolverState{StateT}
    state_k::StateT
    state_kp1::StateT
    state_θ::StateT
end

function update_solver_state!(state::MoreauSolverState, traj, timestep)
    state.state_k = traj[timestep]
    state.state_kp1 = traj[timestep+1]
    state.state_θ = deepcopy(state.state_k)
    state.state_θ.t = traj.t[timestep+1]
    state
end

function compute_moreau_residual!(
        res_out, xvec,
        solver_state::MoreauSolverState,
        consts::MoreauConstants,
        ws::MoreauJacobianWorkspace,
        structure, bot, field, policy,
        Ḿ, θ
    )
    (;h, mass_norm, nq, nλ) = consts
    state_k = solver_state.state_k
    state_kp1 = solver_state.state_kp1
    state_θ = solver_state.state_θ

    copyto!(state_kp1.q, 1, xvec, 1, nq)
    copyto!(state_kp1.λ, 1, xvec, nq+1, nλ)
    @. state_θ.q̇ = (state_kp1.q - state_k.q)/h
    @. state_kp1.q̇ = (1/θ) * state_θ.q̇ - (1/θ - 1) * state_k.q̇
    @. state_θ.q = (1-θ) * state_k.q + θ * state_kp1.q
    state_θ.λ .= state_kp1.λ

    gen_force!(state_θ, bot, field, policy)
    state_kp1.F .= state_θ.F
    copyto!(ws.A_buffer, cstr_jacobian(structure, state_kp1))

    res_out[   1:nq   ] .= h .* Ḿ * (state_kp1.q̇ .- state_k.q̇) .-
                            h^2 .* state_θ.F .+
                            mass_norm .* transpose(ws.A_buffer) * state_kp1.λ
    res_out[nq+1:nq+nλ] .= mass_norm .* h .* ws.A_buffer * state_kp1.q̇
    res_out
end

function compute_moreau_jacobian!(
        jac_out, xvec,
        solver_state::MoreauSolverState,
        consts::MoreauConstants,
        ws::MoreauJacobianWorkspace,
        structure, bot, field, policy,
        Ḿ, M̌, θ
    )
    (;h, mass_norm, nq, nλ) = consts
    state_k = solver_state.state_k
    state_kp1 = solver_state.state_kp1
    state_θ = solver_state.state_θ
    (;∂F∂q, ∂F∂q̇, ∂F∂u, ∂F∂s, A_buffer, Aq̇_buffer, ∂Aᵀλ_buffer) = ws

    copyto!(state_kp1.q, 1, xvec, 1, nq)
    copyto!(state_kp1.λ, 1, xvec, nq+1, nλ)
    @. state_θ.q̇ = (state_kp1.q - state_k.q)/h
    @. state_kp1.q̇ = (1/θ) * state_θ.q̇ - (1/θ - 1) * state_k.q̇
    @. state_θ.q = (1-θ) * state_k.q + θ * state_kp1.q
    state_θ.λ .= state_kp1.λ

    gen_force_state_jacobian!(∂F∂q, ∂F∂q̇, ∂F∂u, bot, field, policy, state_θ, ∂F∂s)
  
    copyto!(A_buffer, cstr_jacobian(structure, state_kp1))
    cstr_forces_jacobian!(∂Aᵀλ_buffer, structure, state_kp1)
    cstr_velocity_jacobian!(Aq̇_buffer, structure, state_kp1)

    jac_out[   1:nq ,   1:nq ] .=  1/θ .*M̌ .-
                                   (h^2) .*(θ .*∂F∂q .+ 1/h .*∂F∂q̇) .+
                                   mass_norm .* ∂Aᵀλ_buffer
    jac_out[   1:nq ,nq+1:end] .=  mass_norm .* transpose(A_buffer)
    jac_out[nq+1:end,   1:nq ] .=  mass_norm .* (h .* Aq̇_buffer .+ 1/θ .* A_buffer)
    jac_out[nq+1:end,nq+1:end] .=  0.0
    jac_out
end


function generate_cache(
        prob::DynamicsProblem,
        solver::DynamicsSolver{
            <:Moreau,
        },
        ::Val{true};
        dt,kargs...
    )
    (;bot) = prob
    (;traj,structure) = bot
    (;M,M⁻¹,M̌,M̌⁻¹,Ḿ,M̄)= build_mass_matrices(structure)
    q0 = traj.q[begin]
    λ0 = traj.λ[begin]
    q̇0 = traj.q̇[begin]
    p0 = traj.p[begin] .= Ḿ*q̇0
    T = eltype(q0)
    nq = length(q0)
    nλ = length(λ0)
    nu = get_num_of_actions(bot)
    ns = get_num_of_aux_var(structure)
    ∂F∂q = zeros(T,nq,nq)
    ∂F∂q̇ = zeros(T,nq,nq)
    ∂F∂u = zeros(T,nq,nu)
    ∂F∂s = zeros(T,nq,ns)
    nx = nq + nλ
    mass_norm = norm(Ḿ, Inf)
    consts = MoreauConstants(dt, mass_norm, nq, nλ, nx)
    A_buffer = zeros(T, nλ, nq)
    Aq̇_buffer = similar(A_buffer)
    ∂Aᵀλ_buffer = zeros(T, nq, nq)
    jac_workspace = MoreauJacobianWorkspace(
        ∂F∂q, ∂F∂q̇, ∂F∂u, ∂F∂s,
        A_buffer, Aq̇_buffer, ∂Aᵀλ_buffer
    )
    options = merge(
        (use_fd_jacobian=false,), #default
        prob.options,
        solver.options,
    )
    Moreau_Constant_Mass_Cache(solver,
        @eponymtuple(
            Ḿ,M̌,M̄,M̌⁻¹,
            jac_workspace,
            consts,
            options
        )
    )
end

function solve!(sim::Simulator,cache::Moreau_Constant_Mass_Cache;
        dt,ftol=1e-14,verbose=false,maxiters=50,
        progress=true,exception=true,
        use_fd_jacobian=nothing)
    (;prob,controller,tspan,restart,totalstep) = sim
    (;bot, policy, env) = prob
    (;traj) = bot
    (;
        Ḿ,M̌,M̄,M̌⁻¹,
        jac_workspace,
        consts,
        options
    ) = cache.cache
    structure = bot.structure
    use_fd_jacobian = something(use_fd_jacobian, get(options, :use_fd_jacobian, false))
    (;θ) = cache.solver.integrator
    (;h, mass_norm, nq, nλ, nx) = consts
    q0 = traj.q[begin]
    T = eltype(q0)
    newton_workspace = NewtonWorkspace(T, nx, 0, 0)
    ls_workspace = NewtonWorkspace(
        newton_workspace.xₖ, newton_workspace.Res, newton_workspace.Jac,
        newton_workspace.Δx, newton_workspace.x, newton_workspace.𝐰,
        newton_workspace.∂Γ∂x
    )

    iteration_break = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    dg_step = ceil(Int,log10(totalstep))+1
    dg_dt = max(1,-floor(Int,log10(h)))
    wd_t = ceil(Int,log10(traj.t[end]))+dg_dt+1+1
    progfmt = Printf.Format("Progress: %5.1f%%, step: %$(dg_step)u, time: %$(wd_t).$(dg_dt)f, iterations: %s \n")
    x = newton_workspace.x
    Res = newton_workspace.Res
    Jac = newton_workspace.Jac
    Δx = newton_workspace.Δx

    for timestep = 1:totalstep
        #---------Step k Control-----------
        
        #---------Step k Control-----------
        solver_state = MoreauSolverState(
            traj[timestep],
            traj[timestep+1],
            deepcopy(traj[timestep])
        )
        (;state_k,state_kp1,state_θ) = solver_state
        x[1:nq] .= state_k.q .+ h.*state_k.q̇
        x[nq+1:nq+nλ] .= 0

        residual! = (res_out, xvec) -> compute_moreau_residual!(
            res_out, xvec, solver_state, consts, jac_workspace, structure, bot, env.field, policy, Ḿ, θ
        )

        jacobian! = (jac_out, xvec) -> compute_moreau_jacobian!(
            jac_out, xvec, solver_state, consts, jac_workspace, structure, bot, env.field, policy, Ḿ, M̌, θ
        )

        call_res! = ws -> residual!(ws.Res, ws.x)
        call_jac! = ws -> begin
            if use_fd_jacobian
                FiniteDiff.finite_difference_jacobian!(ws.Jac, residual!, ws.x, Val{:central})
            else
                jacobian!(ws.Jac, ws.x)
            end
            nothing
        end

        isconverged = false
        normRes = typemax(T)
        for iteration = 1:maxiters
            residual!(Res, x)
            normRes = norm(Res)
            iteration_break = iteration-1
            if normRes < ftol
                isconverged = true
                break
            end                

            if use_fd_jacobian
                FiniteDiff.finite_difference_jacobian!(Jac, residual!, x, Val{:central})
            else
                jacobian!(Jac, x)
            end

            luJac = lu(Jac)
            Δx .= luJac \ (-Res)

            ϕ, dϕ, ϕdϕ = create_line_search_functions(
                newton_workspace, ls_workspace, call_res!, call_jac!
            )

            fx = 0.5 * normRes^2
            mul!(newton_workspace.JacΔx, Jac, Δx)
            dϕ_0 = dot(Res, newton_workspace.JacΔx)
            α, fx = Rible.LS.BackTracking(order=3)(ϕ, dϕ, ϕdϕ, 1.0, fx, dϕ_0)

            x .+= α .* Δx
        end

        if !isconverged
            if exception
                error("Not Converged! Step=$timestep")
            else
                # sim.convergence = false
                break
            end
        end
        copyto!(state_kp1.q, 1, x, 1, nq)
        copyto!(state_kp1.λ, 1, x, nq+1, nλ)
        @. state_θ.q̇ = (state_kp1.q - state_k.q)/h
        @. state_kp1.q̇ = (1/θ) * state_θ.q̇ - (1/θ - 1) * state_k.q̇
        # q̇ₖ .= invM̌*(pₖ.-M̄*q̃̇ₖ )
        #---------Step k finisher-----------
        #---------Step k finisher-----------
        if verbose
            progstr = Printf.format(
                progfmt,
                floor(timestep/totalstep*100;digits=1), timestep, state_k.t, iteration_break
            )
            print(progstr)
        end
        next!(prog)
    end

    bot
end
