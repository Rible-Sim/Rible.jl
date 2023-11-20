struct Zhong06QCache{CacheType}
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
    M = assemble_M(structure) |> Matrix
    ∂Mq̇∂q = assemble_∂Mq̇∂q(structure)
    M! = make_M!(structure)
    Jac_M! = make_Jac_M!(structure)
    Φ = make_cstr_function(structure)
    A = make_cstr_jacobian(structure)
    F!(F,q,q̇,t) = generalize_force!(F,bot,q,q̇,t)
    Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t) = generalize_force_jacobain!(∂F∂q̌,∂F∂q̌̇,bot,q,q̇,t)
    cache = @eponymtuple(F!,Jac_F!,M,∂Mq̇∂q,M!,Jac_M!,Φ,A)
    Zhong06QCache(cache)
end

function solve!(simulator::Simulator,solvercache::Zhong06QCache;
                dt,ftol=1e-14,verbose=false,maxiters=50,
                progress=true,exception=true)
    (;prob,controller,tspan,restart,totalstep) = simulator
    (;bot,) = prob
    (;structure,traj) = bot
    (;F!,Jac_F!,M,∂Mq̇∂q,M!,Jac_M!,Φ,A) = solvercache.cache
    q0 = traj.q[begin]
    λ0 = traj.λ[begin]
    q̇0 = traj.q̇[begin]
    T = eltype(q0)
    nq = length(q0)
    nλ = length(λ0)
    ∂F∂q = zeros(T,nq,nq)
    ∂F∂q̇ = zeros(T,nq,nq)
    M!(M,q0)
    pᵏ⁻¹ = M*q̇0
    pᵏ   = zero(pᵏ⁻¹)
    step = 0
    initial_x = vcat(q0,λ0)
    initial_Res = zero(initial_x)
    mr = norm(M,Inf)
    scaling = mr

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
                               scaling.*Aᵀ*λₘ
            Res[nq+1:nq+nλ] .= scaling.*Φ(qᵏ)
        end
    end

    function make_Jac_stepk(qᵏ,qᵏ⁻¹,Mₘ,∂Mq̇∂q,∂F∂q,∂F∂q̇,Aᵀ,tᵏ⁻¹)
        @inline @inbounds function inner_Jac_stepk!(Jac,x)
            h = dt
            qᵏ .= x[1:nq]
            tₘ = tᵏ⁻¹+h/2
            qₘ = (qᵏ.+qᵏ⁻¹)./2
            q̇ₘ = (qᵏ.-qᵏ⁻¹)./h
            M!(Mₘ,qₘ)
            Jac_M!(∂Mq̇∂q,qₘ,qᵏ)
            Jac_F!(∂F∂q,∂F∂q̇,qₘ,q̇ₘ,tₘ)
            Jac[   1:nq ,   1:nq ] .=  Mₘ .+ ∂Mq̇∂q.-(h^2)/2 .*(1/2 .*∂F∂q.+1/h.*∂F∂q̇)
            Jac[   1:nq ,nq+1:end] .=  scaling.*Aᵀ
            # Jac[nq+1:end,   1:nq ] .=  scaling.*transpose(Aᵀ)
            Jac[nq+1:end,   1:nq ] .=  scaling.*A(qᵏ)
            Jac[nq+1:end,nq+1:end] .=  0.0
        end
    end

    @inline @inbounds function Momentum_k!(pᵏ,pᵏ⁻¹,qᵏ,qᵏ⁻¹,λᵏ,Mₘ,A,Aᵀ,h)
        pᵏ .= -pᵏ⁻¹.+2/h.*Mₘ*(qᵏ.-qᵏ⁻¹) .-
            1/h.*scaling.*(transpose(A(qᵏ))-Aᵀ)*λᵏ
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
        Res_stepk! = make_Res_stepk(qᵏ,λᵏ,qᵏ⁻¹,pᵏ⁻¹,M,F,Aᵀ,tᵏ⁻¹)
        if Jac_F! isa Missing
            dfk = OnceDifferentiable(Res_stepk!,initial_x,initial_Res)
        else

            Jac_stepk! = make_Jac_stepk(qᵏ,qᵏ⁻¹,M,∂Mq̇∂q,∂F∂q,∂F∂q̇,Aᵀ,tᵏ⁻¹)
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
        M!(M,(qᵏ.+qᵏ⁻¹)./2)
        Momentum_k!(pᵏ,pᵏ⁻¹,qᵏ,qᵏ⁻¹,λᵏ,M,A,Aᵀ,dt)           
        M!(M,qᵏ)
        q̇ᵏ .= inv(Matrix(M))*(pᵏ)
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
