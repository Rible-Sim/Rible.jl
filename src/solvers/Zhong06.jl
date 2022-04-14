struct Zhong06 end

struct Zhong06Cache{MMT,AT,ΦT}
    mass_matrices::MMT
    A::AT
    Φ::ΦT
end

function generate_cache(::Zhong06,intor;dt,kargs...)
    (;prob,state) = intor
    (;bot,dynfuncs) = prob
    (;q,q̇) = state.now
    # F!,_ = dynfuncs
    mm = build_MassMatrices(bot)
    # (;M) = mm
    A = make_A(bot)
    Φ = make_Φ(bot)
    Zhong06Cache(mm,A,Φ)
end

function retrieve!(intor,cache::Zhong06Cache)

end

function solve!(intor::Integrator,cache::Zhong06Cache;
                dt,ftol=1e-14,verbose=false,iterations=50,
                progress=true,exception=true)
    (;prob,state,control!,tspan,restart,totalstep) = intor
    (;bot,dynfuncs) = prob
    (;traj) = bot
    # @unpack t,q,q̇,tprev,qprev,q̇prev = state
    (;F!,Jac_F!) = dynfuncs
    (;mass_matrices,A,Φ) = cache
    (;Ḿ,M̌,M̄,invM̌) = mass_matrices
    q̌0 = traj.q̌[begin]
    λ0 = traj.λ[begin]
    q̇0 = traj.q̇[begin]
    T = eltype(q̌0)
    nq̌ = length(q̌0)
    nλ = length(λ0)
    ∂F∂q̌ = zeros(T,nq̌,nq̌)
    ∂F∂q̌̇ = zeros(T,nq̌,nq̌)
    p̌ᵏ⁻¹ = Ḿ*q̇0
    p̌ᵏ   = zero(p̌ᵏ⁻¹)
    step = 0
    initial_x = vcat(q̌0,λ0)
    initial_Res = zero(initial_x)
    mr = norm(Ḿ,Inf)
    scaling = mr

    function make_Res_stepk(qᵏ,q̌ᵏ,λᵏ,qᵏ⁻¹,p̌ᵏ⁻¹,F̌,Aᵀ,tᵏ⁻¹)
        @inline @inbounds function inner_Res_stepk!(Res,x)
            h = dt
            q̌ᵏ .= x[   1:nq̌   ]
            λᵏ .= x[nq̌+1:nq̌+nλ]
            F!(F̌,(qᵏ.+qᵏ⁻¹)./2,(qᵏ.-qᵏ⁻¹)./h,tᵏ⁻¹+h/2)
            Res[   1:nq̌   ] .= Ḿ*(qᵏ.-qᵏ⁻¹) .-
                               h.*p̌ᵏ⁻¹ .-
                               (h^2)/2 .*F̌ .-
                               scaling.*Aᵀ*λᵏ
            Res[nq̌+1:nq̌+nλ] .= scaling.*Φ(qᵏ)
        end
    end

    function make_Jac_stepk(qᵏ,q̌ᵏ,qᵏ⁻¹,∂F∂q̌,∂F∂q̌̇,Aᵀ,tᵏ⁻¹)
        @inline @inbounds function inner_Jac_stepk!(Jac,x)
            h = dt
            q̌ᵏ .= x[1:nq̌]
            Jac_F!(∂F∂q̌,∂F∂q̌̇,(qᵏ.+qᵏ⁻¹)./2,(qᵏ.-qᵏ⁻¹)./h,tᵏ⁻¹+h/2)
            Jac[   1:nq̌ ,   1:nq̌ ] .=  M̌.-(h^2)/2 .*(1/2 .*∂F∂q̌.+1/h.*∂F∂q̌̇)
            Jac[   1:nq̌ ,nq̌+1:end] .= -scaling.*Aᵀ
            Jac[nq̌+1:end,   1:nq̌ ] .=  scaling.*A(qᵏ)
            Jac[nq̌+1:end,nq̌+1:end] .=  0.0
        end
    end

    @inline @inbounds function Momentum_k!(p̌ᵏ,p̌ᵏ⁻¹,qᵏ,qᵏ⁻¹,λᵏ,M,A,Aᵀ,h)
        p̌ᵏ .= -p̌ᵏ⁻¹.+2/h.*Ḿ*(qᵏ.-qᵏ⁻¹) .+
            1/h.*scaling.*(transpose(A(qᵏ))-Aᵀ)*λᵏ
    end

    total_iterations = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        # control!(intor,cache)
        #---------Step k Control-----------
        tᵏ⁻¹ = traj.t[timestep]
        tᵏ = traj.t[timestep+1]
        qᵏ⁻¹ = traj.q[timestep]
        qᵏ = traj.q[timestep+1]
        q̌ᵏ = traj.q̌[timestep+1]
        q̌̇ᵏ = traj.q̌̇[timestep+1]
        q̃̇ᵏ = traj.q̃̇[timestep+1]
        λᵏ = traj.λ[timestep+1]
        F̌ = traj.F̌[timestep+1]
        initial_x[   1:nq̌]    .= traj.q̌[timestep]
        initial_x[nq̌+1:nq̌+nλ] .= traj.λ[timestep]
        Aᵀ = transpose(A(qᵏ⁻¹))
        Res_stepk! = make_Res_stepk(qᵏ,q̌ᵏ,λᵏ,qᵏ⁻¹,p̌ᵏ⁻¹,F̌,Aᵀ,tᵏ⁻¹)
        if Jac_F! isa Nothing
            dfk = OnceDifferentiable(Res_stepk!,initial_x,initial_Res)
        else

            Jac_stepk! = make_Jac_stepk(qᵏ,q̌ᵏ,qᵏ⁻¹,∂F∂q̌,∂F∂q̌̇,Aᵀ,tᵏ⁻¹)
            # Jac_ref = zeros(nq̌+nλ,nq̌+nλ)
            # FiniteDiff.finite_difference_jacobian!(Jac_ref,Res_stepk!,initial_x)
            # Jac_my = zeros(nq̌+nλ,nq̌+nλ)
            # Jac_stepk!(Jac_my,initial_x)
            # diff_Jac = Jac_my .- Jac_ref
            # diff_Jac[abs.(diff_Jac).<1e-5] .= 0.0
            # display(diff_Jac)
            # @show maximum(abs.(diff_Jac))
        # end
            dfk = OnceDifferentiable(Res_stepk!,Jac_stepk!,initial_x,initial_Res)
        end
        Res_stepk_result = nlsolve(dfk, initial_x; ftol, iterations, method=:newton)

        if converged(Res_stepk_result) == false
            if exception
                error("Not Converged! Step=$timestep")
            else
                # intor.convergence = false
                break
            end
        end
        total_iterations += Res_stepk_result.iterations
        xᵏ = Res_stepk_result.zero
        q̌ᵏ .= xᵏ[   1:nq̌   ]
        λᵏ .= xᵏ[nq̌+1:nq̌+nλ]
        Momentum_k!(p̌ᵏ,p̌ᵏ⁻¹,qᵏ,qᵏ⁻¹,λᵏ,Ḿ,A,Aᵀ,dt)
        q̌̇ᵏ .= invM̌*(p̌ᵏ.-M̄*q̃̇ᵏ )
        #---------Step k finisher-----------
        step += 1
        p̌ᵏ,p̌ᵏ⁻¹ = p̌ᵏ⁻¹,p̌ᵏ
        state.prv = traj[timestep]
        state.now = traj[timestep+1]
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

    return intor,cache
end
