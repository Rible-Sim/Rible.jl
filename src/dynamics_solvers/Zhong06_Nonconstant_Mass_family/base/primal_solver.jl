function generate_cache(
        prob::DynamicsProblem,
        solver::DynamicsSolver{
            <:Zhong06,
        },
        ::Val{false};
        dt, kargs...
    )
    @info "Zhong06_Nonconstant_Mass_Cache"
    (; bot, policy, env) = prob
    (; structure, traj) = bot
    Mₘ = assemble_M(structure)
    mass_norm = norm(Mₘ, Inf)
    consts = Zhong06Constants(bot, policy, structure, mass_norm, dt)
    jacobian_workspace = Zhong06JacobianWorkspace(bot)
    
    options = merge(
        (use_fd_jacobian=false,), #default
        prob.options,
        solver.options,
    )
    Zhong06_Nonconstant_Mass_Cache(
        bot, policy, env,
        jacobian_workspace,
        consts, options
    )
end

function pack_nonconstant_state!(solver_state::Zhong06SolverState, x, consts::Zhong06Constants)
    (; nq, nλ) = consts
    h = solver_state.dt
    n1 = nq
    n2 = n1 + nλ
    qₖ = solver_state.qₖ
    qₖ₊₁ = solver_state.qₖ₊₁
    λₘ = solver_state.λₘ
    qₘ = solver_state.qₘ
    q̇ₘ = solver_state.q̇ₘ

    qₖ₊₁ .= @view x[1:n1]
    λₘ .= @view x[n1+1:n2]

    @. qₘ = (qₖ₊₁ + qₖ) / 2
    @. q̇ₘ = (qₖ₊₁ - qₖ) / h
    solver_state.state_mid.t = (solver_state.tₖ + solver_state.tₖ₊₁) / 2
    nothing
end


function solve!(simulator::Simulator, solvercache::Zhong06_Nonconstant_Mass_Cache;
    dt, ftol=1e-14, verbose=false, maxiters=50,
    progress=true, exception=true,
    use_fd_jacobian=true)
    @info "Solving Zhong06_Nonconstant_Mass_Cache"
    (; prob, controller, tspan, restart, totalstep) = simulator
    (; bot) = prob
    (; structure, traj) = bot
    (; consts, jacobian_workspace) = solvercache
    (; Mₘ, M⁻¹ₘ, ∂Mₘhq̇ₘ∂qₘ, Aₖ, Aₖ₊₁) = jacobian_workspace
    q0 = traj.q[begin]
    q̇0 = traj.q̇[begin]
    T = eltype(q0)
    (; h, mass_norm, nq, nλ) = consts
    dt = h
    nx = nq + nλ
    assemble_M!(Mₘ, structure)
    p0 = traj.p[begin]
    mul!(p0, Mₘ, q̇0)
    step = 0
    newton_workspace = NewtonWorkspace(T, nx, 0, 0)
    ls_workspace = NewtonWorkspace(
        newton_workspace.xₖ, newton_workspace.Res, newton_workspace.Jac,
        newton_workspace.Δx, newton_workspace.x, newton_workspace.𝐰,
        newton_workspace.∂Γ∂x
    )
    xₖ = newton_workspace.x
    Res = newton_workspace.Res
    Jac = newton_workspace.Jac
    Δx = newton_workspace.Δx

    total_iterations = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    state_mid = deepcopy(traj[begin])
    for timestep = 1:totalstep
        #---------Step k Control-----------

        #---------Step k Control-----------
        solver_state = Zhong06SolverState(
            traj[timestep], traj[timestep+1], state_mid,
            h
        )
        guess_newton!(xₖ, solver_state, consts)
        call_res! = (ws) -> compute_nonconstant_mass_residual!(
            ws.Res, ws.x, solver_state, solvercache
        )
        call_jac! = (ws) -> begin
            if use_fd_jacobian
                fd_residual! = (res_out, xvec) -> compute_nonconstant_mass_residual!(
                    res_out, xvec, solver_state, solvercache
                )
                FiniteDiff.finite_difference_jacobian!(ws.Jac, fd_residual!, ws.x, Val{:central})
            else
                compute_nonconstant_mass_jacobian!(ws.Jac, ws.x, solver_state, solvercache)
            end
            nothing
        end
        isconverged = false
        normRes = typemax(T)
        iteration_break = 0
        # Newton iteration
        for iteration = 1:maxiters
            @debug "Iteration $iteration Begin"
            compute_nonconstant_mass_residual!(Res, xₖ, solver_state, solvercache)
            normRes = norm(Res)

            iteration_break = iteration - 1
            if normRes < ftol
                isconverged = true
                @debug "Iteration $iteration Break" normRes
                break
            end
            # @show timestep, iteration, norm(Res)
            if use_fd_jacobian
                fd_residual! = (res_out, xvec) -> compute_nonconstant_mass_residual!(
                    res_out, xvec, solver_state, solvercache
                )
                FiniteDiff.finite_difference_jacobian!(Jac, fd_residual!, xₖ, Val{:central})
            else
                compute_nonconstant_mass_jacobian!(Jac, xₖ, solver_state, solvercache)
            end
            # @show rank(Jac), size(Jac)
            luJac = lu(Jac)
            Δx .= luJac \ (-Res)
            ϕ, dϕ, ϕdϕ = create_line_search_functions(
                newton_workspace, ls_workspace, call_res!, call_jac!
            )

            fx = 0.5 * normRes^2
            mul!(newton_workspace.JacΔx, Jac, Δx)
            dϕ_0 = dot(Res, newton_workspace.JacΔx)
            #users should most likely use MoreThuente, HagerZhang or BackTracking
            α, fx = Rible.LS.BackTracking(order=3)(ϕ, dϕ, ϕdϕ, 1.0, fx, dϕ_0)
            @debug "Line Search" α, fx
            xₖ .+= α * Δx
            @debug "Iteration $iteration End  " normRes

        end
        if !isconverged
            @warn "Newton max iterations $maxiters, Step=$timestep, normRes=$normRes"
            if exception
                error("Not Converged!")
            else
                # sim.convergence = false
                ## break
            end
        end
        total_iterations += iteration_break
        populate!(solver_state,newton_workspace.x,structure,jacobian_workspace,consts)
    
        #---------Step k finisher-----------
        step += 1
        #---------Step k finisher-----------
        @debug "Timestep $timestep End  "
        next!(prog)
    end

    bot
end
