struct Zhong06_Nonconstant_Mass_Cache{CacheType}
    cache::CacheType
end

function generate_cache(
        simulator::Simulator{<:DynamicsProblem},
        solver::DynamicsSolver{Zhong06,Nothing},
        ::Val{false};
        dt,kargs...
    )
    (;bot) = simulator.prob
    (;structure) = bot
    Mₘ = assemble_M(structure)
    M⁻¹ₖ = assemble_M⁻¹(structure)
    ∂Mₘhq̇ₘ∂qₘ = assemble_∂Mq̇∂q(structure)
    M! = make_M!(structure)
    M⁻¹! = make_M⁻¹!(structure)
    M_and_Jac_M! = make_M_and_Jac_M!(structure)
    Φ = make_cstr_function(structure)
    A = make_cstr_jacobian(structure)
    F!(F,q,q̇,t) = generalized_force!(F,bot,q,q̇,t)
    Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t) = generalized_force_jacobian!(∂F∂q̌,∂F∂q̌̇,bot,q,q̇,t)
    cache = @eponymtuple(
        F!,Jac_F!,
        Mₘ,M⁻¹ₖ,
        ∂Mₘhq̇ₘ∂qₘ,
        M!,M⁻¹!,
        M_and_Jac_M!,
        Φ,A)
    Zhong06_Nonconstant_Mass_Cache(cache)
end

function solve!(simulator::Simulator,solvercache::Zhong06_Nonconstant_Mass_Cache;
                dt,ftol=1e-14,verbose=false,maxiters=50,
                progress=true,exception=true)
    (;prob,controller,tspan,restart,totalstep) = simulator
    (;bot,) = prob
    (;structure,traj) = bot
    (;
        F!,Jac_F!,
        Mₘ,M⁻¹ₖ,∂Mₘhq̇ₘ∂qₘ,
        M!,M⁻¹!,M_and_Jac_M!,
        Φ,A
    ) = solvercache.cache
    q0 = traj.q[begin]
    λ0 = traj.λ[begin]
    q̇0 = traj.q̇[begin]
    T = eltype(q0)
    nq = length(q0)
    nλ = length(λ0)
    ∂F∂q = zeros(T,nq,nq)
    ∂F∂q̇ = zeros(T,nq,nq)
    M!(Mₘ,q0)
    pᵏ⁻¹ = Mₘ*q̇0
    pᵏ   = zero(pᵏ⁻¹)
    step = 0
    initial_x = vcat(q0,λ0)
    initial_Res = zero(initial_x)
    mr = norm(Mₘ,Inf)
    mass_norm = mr

    function make_Res_stepk(qᵏ,λₘ,qᵏ⁻¹,pᵏ⁻¹,Mₘ,Fₘ,Aᵀ,tᵏ⁻¹)
        @inline @inbounds function inner_Res_stepk!(Res,x)
            h = dt
            qᵏ .= x[   1:nq   ]
            λₘ .= x[nq+1:nq+nλ]
            tₘ = tᵏ⁻¹+h/2
            qₘ = (qᵏ.+qᵏ⁻¹)./2
            q̇ₘ = (qᵏ.-qᵏ⁻¹)./h
            M!(Mₘ,qₘ)
            F!(Fₘ,qₘ,q̇ₘ,tₘ)
            Res[   1:nq   ] .= Mₘ*(qᵏ.-qᵏ⁻¹) .-
                               h.*pᵏ⁻¹ .-
                               (h^2)/2 .*Fₘ .+
                               mass_norm.*Aᵀ*λₘ
            Res[nq+1:nq+nλ] .= mass_norm.*Φ(qᵏ)
        end
    end

    function make_Jac_stepk(qᵏ,qᵏ⁻¹,Mₘ,∂Mₘhq̇ₘ∂qₘ,∂F∂q,∂F∂q̇,Aᵀ,tᵏ⁻¹)
        @inline @inbounds function inner_Jac_stepk!(Jac,x)
            h = dt
            qᵏ .= x[1:nq]
            tₘ = tᵏ⁻¹+h/2
            qₘ = (qᵏ.+qᵏ⁻¹)./2
            q̇ₘ = (qᵏ.-qᵏ⁻¹)./h
            M_and_Jac_M!(Mₘ,∂Mₘhq̇ₘ∂qₘ,qₘ,h.*q̇ₘ)
            Jac_F!(∂F∂q,∂F∂q̇,qₘ,q̇ₘ,tₘ)
            Jac[   1:nq ,   1:nq ] .=  Mₘ .+ 1/2 .*∂Mₘhq̇ₘ∂qₘ.-(h^2)/2 .*(1/2 .*∂F∂q.+1/h.*∂F∂q̇)
            Jac[   1:nq ,nq+1:end] .=  mass_norm.*Aᵀ
            Jac[nq+1:end,   1:nq ] .=  mass_norm.*A(qᵏ)
            Jac[nq+1:end,nq+1:end] .=  0.0
        end
    end

    total_iterations = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        # control!(sim,cache)
        #---------Step k Control-----------
        tᵏ⁻¹ = traj.t[timestep]
        tᵏ = traj.t[timestep+1]
        qᵏ⁻¹ = traj.q[timestep]
        qᵏ = traj.q[timestep+1]
        qᵏ = traj.q[timestep+1]
        q̇ᵏ = traj.q̇[timestep+1]
        λᵏ = traj.λ[timestep+1]
        F = traj.F[timestep+1]
        initial_x[   1:nq]    .= traj.q[timestep]
        initial_x[nq+1:nq+nλ] .= traj.λ[timestep]
        Aᵀ = transpose(A(qᵏ⁻¹))
        Res_stepk! = make_Res_stepk(qᵏ,λᵏ,qᵏ⁻¹,pᵏ⁻¹,Mₘ,F,Aᵀ,tᵏ⁻¹)
        if Jac_F! isa Missing
            dfk = OnceDifferentiable(Res_stepk!,initial_x,initial_Res)
        else

            Jac_stepk! = make_Jac_stepk(qᵏ,qᵏ⁻¹,Mₘ,∂Mₘhq̇ₘ∂qₘ,∂F∂q,∂F∂q̇,Aᵀ,tᵏ⁻¹)
            # Jac_ref = zeros(nq+nλ,nq+nλ)
            # FiniteDiff.finite_difference_jacobian!(Jac_ref,Res_stepk!,initial_x)
            # Jac_my = zeros(nq+nλ,nq+nλ)
            # Jac_stepk!(Jac_my,initial_x)
            # diff_Jac = Jac_my .- Jac_ref
            # diff_Jac[abs.(diff_Jac).<1e-5] .= 0.0
            # display(diff_Jac)
            # @show maximum(abs.(diff_Jac))
        # end
            dfk = OnceDifferentiable(Res_stepk!,Jac_stepk!,initial_x,initial_Res)
        end
        Res_stepk_result = nlsolve(dfk, initial_x; ftol, iterations=maxiters, method=:newton)

        if converged(Res_stepk_result) == false
            if exception
                error("Not Converged! Step=$timestep")
            else
                # sim.convergence = false
                break
            end
        end
        total_iterations += Res_stepk_result.iterations
        xᵏ = Res_stepk_result.zero
        qᵏ .= xᵏ[   1:nq   ]
        λᵏ .= xᵏ[nq+1:nq+nλ]     
        M!(Mₘ,(qᵏ.+qᵏ⁻¹)./2)
        Momentum_k!(pᵏ,pᵏ⁻¹,qᵏ,qᵏ⁻¹,λᵏ,Mₘ,A,Aᵀ,mass_norm,dt)           
        M⁻¹!(M⁻¹ₖ,qᵏ)
        q̇ᵏ .= M⁻¹ₖ*pᵏ
        #---------Step k finisher-----------
        step += 1
        pᵏ,pᵏ⁻¹ = pᵏ⁻¹,pᵏ
        #---------Step k finisher-----------
        if verbose
            dg_step = ceil(Int,log10(totalstep))+1
            dg_dt = max(1,-floor(Int,log10(dt)))
            wd_t = ceil(Int,log10(traj.t[end]))+dg_dt+1+1
            progfmt = Printf.Format("Progress: %5.1f%%, step: %$(dg_step)u, time: %$(wd_t).$(dg_dt)f, iterations: %s \n")
            progstr = Printf.format(progfmt,
                floor(timestep/totalstep*100;digits=1), timestep, tᵏ, Res_stepk_result.iterations
            )
            print(progstr)
        end
        next!(prog)
    end

    bot
end
