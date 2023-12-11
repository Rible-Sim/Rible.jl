struct Moreau_Constant_Mass_Cache{solT,cacheType}
    solver::solT
    cache::cacheType
end

function generate_cache(
        simulator::Simulator{<:DynamicsProblem},
        solver::DynamicsSolver{<:Moreau};
        dt,kargs...
    ) 

    generate_cache(
        simulator,
        solver,
        has_constant_mass_matrix(simulator.prob.bot);
        dt,kargs...
    )

end

function generate_cache(
        simulator::Simulator{<:DynamicsProblem},
        solver::DynamicsSolver{<:Moreau,Nothing},
        ::Val{true};
        dt,kargs...
    )
    (;bot) = simulator.prob
    (;traj,structure) = bot
    (;M,M⁻¹,M̌,M̌⁻¹,Ḿ,M̄)= build_mass_matrices(structure)
    A = make_cstr_jacobian(structure)
    Φ = make_cstr_function(structure)
    F!(F,q,q̇,t) = generalized_force!(F,bot,q,q̇,t)
    Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t) = generalized_force_jacobain!(∂F∂q̌,∂F∂q̌̇,bot,q,q̇,t)
    ∂Aᵀλ∂q(q::AbstractVector,λ) = cstr_forces_jacobian(structure,q,λ)
    ∂Aq̇∂q(q::AbstractVector,q̇) = cstr_velocity_jacobian(structure,q,q̇)
    q̌0 = traj.q̌[begin]
    λ0 = traj.λ[begin]
    q̇0 = traj.q̇[begin]
    q̇0 = traj.q̇[begin]
    p̌0 = traj.p̌[begin] .= Ḿ*q̇0
    T = eltype(q̌0)
    nq̌ = length(q̌0)
    nλ = length(λ0)
    ∂F∂q̌ = zeros(T,nq̌,nq̌)
    ∂F∂q̌̇ = zeros(T,nq̌,nq̌)
    nx = nq̌ + nλ
    initial_x = vcat(q̌0,λ0)
    initial_Res = zero(initial_x)
    initial_Jac = initial_Res*transpose(initial_Res)
    Moreau_Constant_Mass_Cache(solver,
        @eponymtuple(
            F!,Jac_F!,
            Ḿ,M̌,M̄,M̌⁻¹,A,Φ,
            ∂Aᵀλ∂q,∂Aq̇∂q,
            ∂F∂q̌,∂F∂q̌̇,
            nq̌,nλ,nx,
            initial_x,
            initial_Res,
            initial_Jac
        )
    )
end

function retrieve!(sim,cache::Moreau_Constant_Mass_Cache)

end

function solve!(sim::Simulator,cache::Moreau_Constant_Mass_Cache;
                dt,ftol=1e-14,verbose=false,maxiters=50,
                progress=true,exception=true)
    (;prob,controller,tspan,restart,totalstep) = sim
    (;bot,) = prob
    (;traj) = bot
    (;
        F!,Jac_F!,
        Ḿ,M̌,M̄,M̌⁻¹,A,Φ,
        ∂Aᵀλ∂q,∂Aq̇∂q,
        ∂F∂q̌,∂F∂q̌̇,
        nq̌,nλ,nx,
        initial_x,
        initial_Res,
        initial_Jac
    ) = cache.cache
    (;θ) = cache.solver.integrator
    mr = norm(Ḿ,Inf)
    mass_norm = mr
    h = dt
    function make_Res_stepk(qₖ₊₁,q̌ₖ₊₁,vₖ₊₁,λₖ₊₁,qₖ,vₖ,F̌,tₖ₊θ)
        function inner_Res_stepk!(Res,x)
            q̌ₖ₊₁ .= x[   1:nq̌   ]
            λₖ₊₁ .= x[nq̌+1:nq̌+nλ]
            vₖ₊θ = (qₖ₊₁.-qₖ)./h
            vₖ₊₁ .= (1/θ)*vₖ₊θ .- (1/θ-1)*vₖ
            qₖ₊θ = (1-θ)*qₖ.+θ*qₖ₊₁
            F!(F̌,qₖ₊θ,vₖ₊θ,tₖ₊θ)
            Aₖ₊₁ = A(qₖ₊₁)
            Res[   1:nq̌   ] .= h*Ḿ*(vₖ₊₁.-vₖ) .-
                               h^2*F̌ .+
                               mass_norm.*transpose(Aₖ₊₁)*λₖ₊₁
            Res[nq̌+1:nq̌+nλ] .= mass_norm.*h.*Aₖ₊₁*vₖ₊₁
            # Res[nq̌+1:nq̌+nλ] .= mass_norm.*Φ(qₖ₊₁)

        end
    end

    function Jac_stepk!(Jac,x,qₖ₊₁,q̌ₖ₊₁,vₖ₊₁,λₖ₊₁,qₖ,vₖ,∂F∂q̌,∂F∂q̌̇,tₖ₊θ)
        q̌ₖ₊₁ .= x[   1:nq̌   ]
        λₖ₊₁ .= x[nq̌+1:nq̌+nλ]
        vₖ₊θ = (qₖ₊₁.-qₖ)./h
        vₖ₊₁ .= (1/θ)*vₖ₊θ .- (1/θ-1)*vₖ
        qₖ₊θ = (1-θ)*qₖ.+θ*qₖ₊₁
        Jac_F!(∂F∂q̌,∂F∂q̌̇,qₖ₊θ,vₖ₊θ,tₖ₊θ)
        Aₖ₊₁ = A(qₖ₊₁)
        Jac[   1:nq̌ ,   1:nq̌ ] .=  1/θ .*M̌.-(h^2) .*(θ .*∂F∂q̌.+1/h.*∂F∂q̌̇) .+ mass_norm.*∂Aᵀλ∂q(qₖ₊₁,λₖ₊₁)
        Jac[   1:nq̌ ,nq̌+1:end] .=  mass_norm.*transpose(Aₖ₊₁)
        Jac[nq̌+1:end,   1:nq̌ ] .=  mass_norm.*(h.*∂Aq̇∂q(qₖ₊₁,vₖ₊₁) .+ 1/θ.*Aₖ₊₁)
        # Jac[nq̌+1:end,   1:nq̌ ] .=  mass_norm.*Aₖ₊₁
        Jac[nq̌+1:end,nq̌+1:end] .=  0.0
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
        tₖ = traj.t[timestep]
        tₖ₊θ   = traj.t[timestep+1]
        qₖ = traj.q[timestep]
        q̌ₖ = traj.q̌[timestep]
        vₖ = traj.q̇[timestep]
        v̌ₖ = traj.q̌̇[timestep]
        p̌ₖ = traj.p̌[timestep] 
        qₖ₊₁ = traj.q[timestep+1]
        q̌ₖ₊₁ = traj.q̌[timestep+1]
        vₖ₊₁ = traj.q̇[timestep+1]
        q̃̇ₖ₊₁ = traj.q̃̇[timestep+1]
        λₖ₊₁ = traj.λ[timestep+1]
        p̌ₖ₊₁ = traj.p̌[timestep+1]
        F̌ = traj.F̌[timestep+1]
        initial_x[   1:nq̌]    .= q̌ₖ .+ h.*v̌ₖ
        initial_x[nq̌+1:nq̌+nλ] .= 0
        Res_stepk! = make_Res_stepk(qₖ₊₁,q̌ₖ₊₁,vₖ₊₁,λₖ₊₁,qₖ,vₖ,F̌,tₖ₊θ)
        isconverged = false
        if false
            dfk = OnceDifferentiable(Res_stepk!,initial_x,initial_Res)
            Res_stepk_result = nlsolve(dfk, initial_x; ftol, iterations=maxiters, method=:newton)
            isconverged = converged(Res_stepk_result)
            iteration_break = Res_stepk_result.iterations
            Res_stepk_result.iterations
            xₖ₊₁ = Res_stepk_result.zero
        else
            for iteration = 1:maxiters
                Res_stepk!(initial_Res,initial_x)
                # @show initial_Res
                normRes = norm(initial_Res)
                if normRes < ftol
                    isconverged = true
                    iteration_break = iteration-1
                    break
                end                
                Jac_stepk!(initial_Jac,initial_x, qₖ₊₁,q̌ₖ₊₁,vₖ₊₁,λₖ₊₁,qₖ,vₖ,∂F∂q̌,∂F∂q̌̇,tₖ₊θ)
                # @show timestep, iteration
                # Jac_ref = zeros(nq̌+nλ,nq̌+nλ)
                # FiniteDiff.finite_difference_jacobian!(Jac_ref,Res_stepk!,initial_x)
                # diff_Jac = initial_Jac .- Jac_ref
                # diff_Jac[abs.(diff_Jac).<1e-5] .= 0.0
                # diff_rowidx = findall((x)-> x>1e-5,norm.(eachrow(diff_Jac)))
                # @show diff_rowidx #.- (nq̌+nλ)
                # @show diff_Jac[diff_rowidx,:]
                # @show maximum(abs.(diff_Jac))
                initial_x .+= initial_Jac\(-initial_Res)
            end
            xₖ₊₁ = initial_x
            # dfk = OnceDifferentiable(Res_stepk!,Jac_stepk!,initial_x,initial_Res)
        end
        
        if !isconverged
            if exception
                error("Not Converged! Step=$timestep")
            else
                # sim.convergence = false
                break
            end
        end
        q̌ₖ₊₁ .= xₖ₊₁[   1:nq̌   ]
        λₖ₊₁ .= xₖ₊₁[nq̌+1:nq̌+nλ]
        vₖ₊θ = (qₖ₊₁.-qₖ)./h
        vₖ₊₁ .= (1/θ)*vₖ₊θ .- (1/θ-1)*vₖ
        # q̌̇ₖ .= invM̌*(p̌ₖ.-M̄*q̃̇ₖ )
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
