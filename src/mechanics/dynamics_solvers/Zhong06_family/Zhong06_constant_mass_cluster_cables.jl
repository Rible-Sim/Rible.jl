struct Zhong06_Constant_Mass_Cluster_Cables_Cache{solT,cacheType}
    solver::solT
    cache::cacheType
end

function generate_cache(
        simulator::Simulator{<:AbstractDynamicsProblem},
        solver::DynamicsSolver{Zhong06,Nothing,MonolithicApparatusSolver{SmoothedFischerBurmeister}},
        ::Val{true};
        dt,kargs...
    )
    (;prob) = simulator
    (;bot,policy) = prob
    (;traj,structure) = bot
    options = merge(
        (gravity=true,factor=1,checkpersist=true), #default
        prob.options,
        solver.options,
    )
    @show "Yes"
    (;M,M⁻¹,M̌,M̌⁻¹,Ḿ,M̄)= build_mass_matrices(structure)
    A = make_cstr_jacobian(structure)
    Φ = make_cstr_function(structure)
    ψ = build_ψ(structure)
    F!(F,q,s,q̇,t) = generalized_force!(F,bot,policy,q,q̇,t,s;gravity=options.gravity)
    Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t) = generalized_force_jacobian!(∂F∂q̌,∂F∂q̌̇,bot,policy,q,q̇,t)
    q̌0 = traj.q̌[begin]
    λ0 = traj.λ[begin]
    q̇0 = traj.q̇[begin]
    q̇0 = traj.q̇[begin]
    p̌0 = traj.p̌[begin] .= Ḿ*q̇0
    s0 = traj.s[begin]
    T = eltype(q̌0)
    nq̌ = length(q̌0)
    nλ = length(λ0)
    ns = length(s0)
    ∂F∂q̌ = zeros(T,nq̌,nq̌)
    ∂F∂q̌̇ = zeros(T,nq̌,nq̌)
    nx = nq̌ + nλ + ns
    xₖ = vcat(q̌0,λ0,s0)
    Res = zero(xₖ)
    Jac = Res*transpose(Res)
    Zhong06_Constant_Mass_Cluster_Cables_Cache(solver,
        @eponymtuple(
            F!,Jac_F!,
            Ḿ,M̌,M̄,M̌⁻¹,A,Φ,ψ,
            ∂F∂q̌,∂F∂q̌̇,
            nq̌,nλ,nx,
            xₖ,
            Res,
            Jac,
            options
        )
    )
end

function retrieve!(sim,cache::Zhong06_Constant_Mass_Cluster_Cables_Cache)

end

function solve!(sim::Simulator,cache::Zhong06_Constant_Mass_Cluster_Cables_Cache;
                dt,ftol=1e-14,verbose=false,maxiters=50,
                progress=true,exception=true)
    (;prob,controller,tspan,restart,totalstep) = sim
    (;bot,) = prob
    (;traj) = bot
    (;
        F!,Jac_F!,
        Ḿ,M̌,M̄,M̌⁻¹,A,Φ,ψ,
        ∂F∂q̌,∂F∂q̌̇,
        nq̌,nλ,nx,
        xₖ,
        Res,
        Jac
    ) = cache.cache
    mr = norm(Ḿ,Inf)
    mass_norm = mr
    h = dt
    T = get_numbertype(bot)
    function make_Res_stepk(qₖ,q̌ₖ,λₖ,sₖ,qₖ₋₁,p̌ₖ₋₁,F̌,Aᵀₖ₋₁,tₖ₋₁)
        @inline @inbounds function inner_Res_stepk!(Res,x)
            q̌ₖ .= x[   1:nq̌   ]
            λₖ .= x[nq̌+1:nq̌+nλ]
            sₖ .= x[nq̌+nλ+1:nx]
            qₘ = (qₖ.+qₖ₋₁)./2
            q̇ₘ = (qₖ.-qₖ₋₁)./h
            tₘ = tₖ₋₁+h/2
            F!(F̌,qₘ, sₖ, q̇ₘ, tₘ)
            Res[   1:nq̌   ] .= Ḿ*(qₖ.-qₖ₋₁) .-
                               h.*p̌ₖ₋₁ .-
                               (h^2)/2 .*F̌ .+
                               mass_norm.*Aᵀₖ₋₁*λₖ
            Res[nq̌+1:nq̌+nλ] .= mass_norm.*Φ(qₖ)
            Res[nq̌+nλ+1:nx] .= ψ(sₖ)
        end
    end

    function Jac_stepk!(Jac,x,qₖ,sₖ,q̌ₖ,qₖ₋₁,∂F∂q̌,∂F∂q̌̇,Aᵀₖ₋₁,tₖ₋₁)
        q̌ₖ .= x[1:nq̌]
        sₖ .= x[nq̌+nλ+1:end]
        qₘ = (qₖ.+qₖ₋₁)./2
        q̇ₘ = (qₖ.-qₖ₋₁)./h
        tₘ = tₖ₋₁+h/2
        Jac_F!(∂F∂q̌,∂F∂q̌̇,qₘ,q̇ₘ,tₘ)
        ζ = build_ζ(bot.structure)
        ∂ζ∂q = build_∂ζ∂q(bot.structure)
        ∂ζ∂s̄ = build_∂ζ∂s̄(bot.structure)
        n = length(ζ)
        ∂s̄∂s̄ = I(n)
        ∂F∂s̄ = build_∂Q̌∂s̄(bot.structure)
        κ₁ = 10; κ₂ = 10
        coes = diagm([((ζ[i]/κ₁)^2 + (κ₂*sₖ[i])^2)^(-1/2) for i in 1:n])
        coζ = coes*diagm([ζ[i]/κ₁^2 for i in 1:n]) - diagm([1/κ₁ for i in 1:n])
        cos̄ = coes*diagm([κ₂^2*sₖ[i] for i in 1:n]) - diagm([κ₂ for i in 1:n])

        Jac[   1:nq̌ ,   1:nq̌ ] .=  M̌.-(h^2)/2 .*(1/2 .*∂F∂q̌.+1/h.*∂F∂q̌̇)
        Jac[   1:nq̌ ,nq̌+1:nq̌+nλ] .=  mass_norm.*Aᵀₖ₋₁
        Jac[   1:nq̌ ,nq̌+nλ+1:end] .=  -1/2 * dt^2 .* ∂F∂s̄

        Jac[nq̌+1:nq̌+nλ,   1:nq̌ ] .=  mass_norm.*A(qₖ)
        Jac[nq̌+1:nq̌+nλ,   nq̌+1:nq̌+nλ ] .= 0.0 
        Jac[nq̌+1:nq̌+nλ,   nq̌+nλ+1:end ] .= 0.0 

        Jac[nq̌+nλ+1:end,1:nq̌] .= 1/2 * coζ*∂ζ∂q 
        Jac[nq̌+nλ+1:end,nq̌+1:nq̌+nλ] .= 0.0
        Jac[nq̌+nλ+1:end,nq̌+nλ+1:end] .= coζ*∂ζ∂s̄ + cos̄*∂s̄∂s̄
    end

    iteration_break = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    dg_step = ceil(Int,log10(totalstep))+1
    dg_dt = max(1,-floor(Int,log10(h)))
    wd_t = ceil(Int,log10(traj.t[end]))+dg_dt+1+1
    progfmt = Printf.Format("Progress: %5.1f%%, step: %$(dg_step)u, time: %$(wd_t).$(dg_dt)f, iterations: %s \n")
    
    for timestep = 1:totalstep
        #---------Step k Control-----------
        # control!(sim,cache)
        #---------Step k Control-----------
        tₖ₋₁ = traj.t[timestep]
        tₖ   = traj.t[timestep+1]
        qₖ₋₁ = traj.q[timestep]
        sₖ₋₁ = traj.s[timestep]
        q̌ₖ₋₁ = traj.q̌[timestep]
        q̌̇ₖ₋₁ = traj.q̌̇[timestep]
        p̌ₖ₋₁ = traj.p̌[timestep] 
        qₖ = traj.q[timestep+1]
        sₖ = traj.s[timestep+1]
        q̌ₖ = traj.q̌[timestep+1]
        q̌̇ₖ = traj.q̌̇[timestep+1]
        q̃̇ₖ = traj.q̃̇[timestep+1]
        λₖ = traj.λ[timestep+1]
        p̌ₖ = traj.p̌[timestep+1]
        F̌ = traj.F̌[timestep+1]
        xₖ[   1:nq̌]    .= q̌ₖ₋₁ .+ h.*q̌̇ₖ₋₁
        xₖ[nq̌+1:nq̌+nλ] .= 0
        xₖ[nq̌+nλ+1:nx] .= sₖ₋₁
        Aᵀₖ₋₁ = transpose(A(qₖ₋₁))
        Res_stepk! = make_Res_stepk(qₖ,q̌ₖ,λₖ,sₖ₋₁,qₖ₋₁,p̌ₖ₋₁,F̌,Aᵀₖ₋₁,tₖ₋₁)
        isconverged = false
        normRes = typemax(T)
        if false#Jac_F! isa Missing
            dfk = OnceDifferentiable(Res_stepk!,xₖ,Res)
            Res_stepk_result = nlsolve(dfk, xₖ; ftol, iterations=maxiters, method=:newton)
            isconverged = converged(Res_stepk_result)
            iteration_break = Res_stepk_result.iterations
            Res_stepk_result.iterations
            xₖ .= Res_stepk_result.zero
        else
            for iteration = 1:maxiters
                Res_stepk!(Res,xₖ)
                normRes = norm(Res)
                if normRes < ftol
                    isconverged = true
                    iteration_break = iteration-1
                    break
                end                
                # FiniteDiff.finite_difference_jacobian!(Jac,Res_stepk!,xₖ,Val{:central})
                Jac_stepk!(Jac,xₖ,qₖ,sₖ,q̌ₖ,qₖ₋₁,∂F∂q̌,∂F∂q̌̇,Aᵀₖ₋₁,tₖ₋₁)
                xₖ .+= Jac\(-Res)
            end
            # dfk = OnceDifferentiable(Res_stepk!,Jac_stepk!,xₖ,Res)
        end
        
        if !isconverged
            if exception
                error("Not Converged! Step=$timestep, normRes=$normRes")
            else
                # sim.convergence = false
                break
            end
        end
        q̌ₖ .= xₖ[   1:nq̌   ]
        λₖ .= xₖ[nq̌+1:nq̌+nλ]
        sₖ .= xₖ[nq̌+nλ+1:nx]
        Momentum_k!(p̌ₖ,p̌ₖ₋₁,qₖ,qₖ₋₁,λₖ,Ḿ,A,Aᵀₖ₋₁,mass_norm,h)
        q̌̇ₖ .= M̌⁻¹*(p̌ₖ.-M̄*q̃̇ₖ )
        #---------Step k finisher-----------
        #---------Step k finisher-----------
        if verbose
            progstr = Printf.format(
                progfmt,
                floor(timestep/totalstep*100;digits=1), timestep, tₖ, iteration_break
            )
            print(progstr)
        end
        next!(prog)
    end

    bot
end
