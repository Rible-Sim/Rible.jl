struct Alpha{T}
    αm::T
    αf::T
    γ::T
    β::T
end

function Alpha(ρ∞)
    αm = (2ρ∞-1)/(ρ∞+1)
    αf = ρ∞/(ρ∞+1)
    γ = 1/2 + αf - αm
    β = 1/4*(γ+1/2)^2
    Alpha(αm,αf,γ,β)
end

struct AlphaCache{solverT,MMT,funcsT,T}
    solver::solverT
    mass_matrices::MMT
    funcs::funcsT
    β′::T
    γ′::T
end

function generate_cache(solver::Alpha,intor;dt,kargs...)
    (;prob,state) = intor
    (;bot,dynfuncs) = prob
    (;q,q̇) = state.now
    (;αm,αf,γ,β) = solver
    # F!,_ = dynfuncs
    mm = build_MassMatrices(bot)
    # (;M) = mm
    A = make_A(bot)
    Φ = make_Φ(bot)
    ∂Ǎᵀλ∂q̌ = (q,λ)->∂Aᵀλ∂q̌(bot.tg,λ)
    funcs = @eponymtuple(A,Φ,∂Ǎᵀλ∂q̌)
    h = dt
    β′ = (1-αm)/((h^2)*β*(1-αf))
    γ′ = γ/(h*β)
    AlphaCache(solver,mm,funcs,β′,γ′)
end

function retrieve!(intor,cache::AlphaCache)

end

function solve!(intor::Integrator,cache::AlphaCache;
                dt,ftol=1e-14,verbose=false,iterations=50,
                progress=true,exception=true)
    (;prob,state,control!,tspan,restart,totalstep) = intor
    (;bot,dynfuncs) = prob
    (;traj) = bot
    # @unpack t,q,q̇,tprev,qprev,q̇prev = state
    (;F!,Jac_F!) = dynfuncs
    (;solver,mass_matrices,funcs,β′,γ′) = cache
    (;A,Φ,∂Ǎᵀλ∂q̌) = funcs
    (;αm,αf,γ,β) = solver
    (;Ḿ,M̌,M̄,invM̌) = mass_matrices
    q̌0 = traj.q̌[begin]
    λ0 = traj.λ[begin]
    q0 = traj.q[begin]
    q̇0 = traj.q̇[begin]
    q̈0 = traj.q̈[begin]
    q̌̈0 = traj.q̌̈[begin]
    F̌0 = traj.F̌[begin]
    t0 = traj.t[begin]
    T = eltype(q̌0)
    nq̌ = length(q̌0)
    nλ = length(λ0)
    ∂F∂q̌ = zeros(T,nq̌,nq̌)
    ∂F∂q̌̇ = zeros(T,nq̌,nq̌)
    F!(F̌0,q0,q̇0,t0)
    q̌̈0 .= invM̌*F̌0
    aᵏ⁻¹ = copy(q̈0)
    aᵏ = copy(aᵏ⁻¹)
    step = 0
    initial_x = vcat(q̌0,λ0)
    initial_Res = zero(initial_x)
    mr = norm(Ḿ,Inf)
    scaling = (β*dt^2)/mr

    function make_Res_stepk(qᵏ,q̇ᵏ,q̈ᵏ,aᵏ,q̌ᵏ,λᵏ,qᵏ⁻¹,q̇ᵏ⁻¹,q̈ᵏ⁻¹,aᵏ⁻¹,F̌,tᵏ)
        @inline @inbounds function inner_Res_stepk!(Res,x)
            h = dt
            q̌ᵏ .= x[   1:nq̌   ]
            λᵏ .= x[nq̌+1:nq̌+nλ]
            aᵏ .= 1/(h^2)/β.*(qᵏ.-qᵏ⁻¹.-h.*q̇ᵏ⁻¹.-(h^2)*(0.5-β).*aᵏ⁻¹)
            q̇ᵏ .= q̇ᵏ⁻¹ .+ h.*((1-γ).*aᵏ⁻¹ .+ γ.*aᵏ)
            q̈ᵏ .= 1/(1-αf).*((1-αm).*aᵏ .+ αm.*aᵏ⁻¹ .- αf.*q̈ᵏ⁻¹)
            F!(F̌,qᵏ,q̇ᵏ,tᵏ)
            Res[   1:nq̌   ] .= scaling.*Ḿ*q̈ᵏ .-
                               scaling.*F̌ .+
                               transpose(A(qᵏ))*λᵏ
            Res[nq̌+1:nq̌+nλ] .= Φ(qᵏ)
        end
    end

    function make_Jac_stepk(qᵏ,q̇ᵏ,q̈ᵏ,aᵏ,q̌ᵏ,λᵏ,qᵏ⁻¹,q̇ᵏ⁻¹,q̈ᵏ⁻¹,aᵏ⁻¹,F̌,tᵏ)
        @inline @inbounds function inner_Jac_stepk!(Jac,x)
            h = dt
            q̌ᵏ .= x[1:nq̌]
            λᵏ .= x[nq̌+1:nq̌+nλ]
            aᵏ .= 1/(h^2)/β.*(qᵏ.-qᵏ⁻¹.-h.*q̇ᵏ⁻¹.-(h^2)*(0.5-β).*aᵏ⁻¹)
            q̇ᵏ .= q̇ᵏ⁻¹ .+ h.*((1-γ).*aᵏ⁻¹ .+ γ.*aᵏ)
            q̈ᵏ .= 1/(1-αf).*((1-αm).*aᵏ .+ αm.*aᵏ⁻¹ .- αf.*q̈ᵏ⁻¹)
            Jac_F!(∂F∂q̌,∂F∂q̌̇,qᵏ,q̇ᵏ,tᵏ)
            Jac[   1:nq̌ ,   1:nq̌ ] .= scaling.*(β′.*M̌ .- ∂F∂q̌) .+ ∂Ǎᵀλ∂q̌(qᵏ,λᵏ)
            Jac[   1:nq̌ ,nq̌+1:end] .= transpose(A(qᵏ))
            Jac[nq̌+1:end,   1:nq̌ ] .= transpose(Jac[   1:nq̌ ,nq̌+1:end])
            Jac[nq̌+1:end,nq̌+1:end] .= 0.0
        end
    end

    total_iterations = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        # control!(intor,cache)
        #---------Step k Control-----------
        tᵏ⁻¹ = traj.t[timestep]
        qᵏ⁻¹ = traj.q[timestep]
        q̌ᵏ⁻¹ = traj.q̌[timestep]
        q̇ᵏ⁻¹ = traj.q̇[timestep]
        q̌̇ᵏ⁻¹ = traj.q̌̇[timestep]
        q̈ᵏ⁻¹ = traj.q̈[timestep]
        q̌̈ᵏ⁻¹ = traj.q̌̈[timestep]
        tᵏ = traj.t[timestep+1]
        qᵏ = traj.q[timestep+1]
        q̌ᵏ = traj.q̌[timestep+1]
        q̇ᵏ = traj.q̇[timestep+1]
        q̌̇ᵏ = traj.q̌̇[timestep+1]
        q̈ᵏ = traj.q̈[timestep+1]
        # q̌̈ᵏ = traj.q̌̈[timestep+1]
        λᵏ = traj.λ[timestep+1]
        F̌ = traj.F̌[timestep+1]
        q̌ᵏ .= q̌ᵏ⁻¹ .+ dt.*q̌̇ᵏ⁻¹ .+ (dt^2) .*(0.5-β).*q̌̈ᵏ⁻¹
        # q̌̇ᵏ .= q̌̇ᵏ⁻¹ .+ dt.*(1-γ).*q̌̈ᵏ⁻¹
        initial_x[   1:nq̌]    .= q̌ᵏ
        initial_x[nq̌+1:nq̌+nλ] .= 0.0
        Res_stepk! = make_Res_stepk(qᵏ,q̇ᵏ,q̈ᵏ,aᵏ,q̌ᵏ,λᵏ,qᵏ⁻¹,q̇ᵏ⁻¹,q̈ᵏ⁻¹,aᵏ⁻¹,F̌,tᵏ)
        if Jac_F! isa Nothing
            dfk = OnceDifferentiable(Res_stepk!,initial_x,initial_Res)
        else
            Jac_stepk! = make_Jac_stepk(qᵏ,q̇ᵏ,q̈ᵏ,aᵏ,q̌ᵏ,λᵏ,qᵏ⁻¹,q̇ᵏ⁻¹,q̈ᵏ⁻¹,aᵏ⁻¹,F̌,tᵏ)
            dfk = OnceDifferentiable(Res_stepk!,Jac_stepk!,initial_x,initial_Res)
            # Jac_ref = zeros(nq̌+nλ,nq̌+nλ)
            # FiniteDiff.finite_difference_jacobian!(Jac_ref,Res_stepk!,initial_x)
            # Jac_my = zeros(nq̌+nλ,nq̌+nλ)
            # Jac_stepk!(Jac_my,initial_x)
            # diff_Jac = Jac_my .- Jac_ref
            # diff_Jac[abs.(diff_Jac).<1e-5] .= 0.0
            # display(diff_Jac)
            # @show maximum(abs.(diff_Jac))
        end
        Res_stepk_result = nlsolve(dfk, initial_x; ftol, iterations, method=:newton)

        if converged(Res_stepk_result) == false
            if exception
                error("Not Converged! Step=$timestep")
            else
                @warn "Not Converged! Step=$timestep"
                break
            end
        end
        total_iterations += Res_stepk_result.iterations
        xᵏ = Res_stepk_result.zero
        q̌ᵏ .= xᵏ[   1:nq̌   ]
        λᵏ .= xᵏ[nq̌+1:nq̌+nλ]
        #---------Step k finisher-----------
        step += 1
        aᵏ,aᵏ⁻¹ = aᵏ⁻¹,aᵏ
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
