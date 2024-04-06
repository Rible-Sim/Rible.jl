struct Adjoint_Zhong06_Nonconstant_Mass_Cache{CacheType}
    cache::CacheType
end

function generate_cache(
        simulator::Simulator{<:AdjointDynamicsProblem},
        solver::DiscreteAdjointDynamicsSolver{Zhong06},
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
    ∂Aᵀλ∂q(q,λ) = cstr_forces_jacobian(structure,q,λ)
    adjoint_traj = StructArray(
        [
            ComponentArray(
                a = zero(state.q),
                b = zero(state.p),
                μ = zero(state.λ)
            )
            for state in bot.traj
        ]
    )
    cache = @eponymtuple(
        F!,Jac_F!,
        Mₘ,M⁻¹ₖ,
        ∂Mₘhq̇ₘ∂qₘ,
        M!,M⁻¹!,
        M_and_Jac_M!,
        Φ,A,∂Aᵀλ∂q,
        adjoint_traj
    )
    Adjoint_Zhong06_Nonconstant_Mass_Cache(cache)
end

function solve!(simulator::Simulator,solvercache::Adjoint_Zhong06_Nonconstant_Mass_Cache;
                dt,ftol=1e-14,verbose=false,maxiters=50,
                progress=true,exception=true)
    (;prob,controller,tspan,restart,totalstep) = simulator
    (;bot,) = prob
    (;structure,traj) = bot
    (;
        F!,Jac_F!,
        Mₘ,∂Mₘhq̇ₘ∂qₘ,
        M!,M_and_Jac_M!,
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
    step = 0
    mr = norm(Mₘ,Inf)
    mass_norm = mr

    function Jac!(Jacₖ₊₁,Jacₖ,qₖ₊₁,qₖ,λₘ,Mₘ,∂Mₘhq̇ₘ∂qₘ,∂F∂q,∂F∂q̇,Aₖ₊₁,Aₖ,tₖ)
        h = dt
        tₘ = tₖ+h/2
        qₘ = (qₖ₊₁.+qₖ)./2
        q̇ₘ = (qₖ₊₁.-qₖ)./h
        M_and_Jac_M!(Mₘ,∂Mₘhq̇ₘ∂qₘ,qₘ,h.*q̇ₘ)
        Jac_F!(∂F∂q,∂F∂q̇,qₘ,q̇ₘ,tₘ)
        Jacₖ .= 0.0
        Jacₖ[   1:nq ,   1:nq ] .=  1/2 .*∂Mₘhq̇ₘ∂qₘ .- Mₘ .-(h^2)/2 .*(1/2 .*∂F∂q .-1/h.*∂F∂q̇) .- mass_norm.*∂Aᵀλ∂q(qₖ,λₘ)
        Jacₖ[nq+1:2nq,   1:nq ] .= -1/2 .*∂Mₘhq̇ₘ∂qₘ .+ Mₘ .-(h^2)/2 .*(1/2 .*∂F∂q .-1/h.*∂F∂q̇)
        Jacₖ[   1:nq ,nq+1:2nq] .= -I(nq)
        
        Jacₖ₊₁ .= 0.0
        Jacₖ₊₁[    1:nq ,    1:nq ] .=  1/2 .*∂Mₘhq̇ₘ∂qₘ .+ Mₘ .-(h^2)/2 .*(1/2 .*∂F∂q .+1/h.*∂F∂q̇)
        Jacₖ₊₁[ nq+1:2nq,    1:nq ] .= -1/2 .*∂Mₘhq̇ₘ∂qₘ .- Mₘ .-(h^2)/2 .*(1/2 .*∂F∂q .+1/h.*∂F∂q̇) .- mass_norm.*∂Aᵀλ∂q(qₖ₊₁,λₘ)
        Jacₖ₊₁[2nq+1:end,    1:nq ] .=  mass_norm.*Aₖ₊₁
        Jacₖ₊₁[ nq+1:2nq, nq+1:2nq] .=  I(nq)
        Jacₖ₊₁[    1:nq ,2nq+1:end] .= -mass_norm.*transpose(Aₖ)
        Jacₖ₊₁[ nq+1:2nq,2nq+1:end] .= -mass_norm.*transpose(Aₖ₊₁)
    end

    Jac!(Jacₖ₊₁,Jacₖ,qₖ₊₁,qₖ,λₘ,Mₘ,∂Mₘhq̇ₘ∂qₘ,∂F∂q,∂F∂q̇,Aₖ₊₁,Aₖ,tₖ)
    yN = adjoint_traj[end]
    ## yN .= transpose(Jacₖ₊₁)\(-1/2*∂ϕ∂xᵀ-∂S∂xᵀ)
    yN .= transpose(Jacₖ₊₁)\(-∂S∂xᵀ)

    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = totalstep-1:-1:1
        #---------Step k Control-----------
        # control!(sim,cache)
        #---------Step k Control-----------
        tₖ₊₁ = traj.t[timestep+1]
        tₖ = traj.t[timestep]
        tₘ = 1/2*(tₖ₊₁ + tₖ)
        yₖ₊₁ = adjoint_traj[timestep+1]
        yₖ = adjoint_traj[timestep]
        aₖ₊₁ = yₖ₊₁.a
        bₖ₊₁ = yₖ₊₁.b
        μₖ₊₁ = yₖ₊₁.μ
        aₖ  = yₖ.a
        bₖ  = yₖ.b
        μₖ  = yₖ.μ
        ## aₖ₊₁ == bₖ
        Aₖ₊₁ .= A(qₖ₊₁)
        Aₖ .= A(qₖ)
        Jac!(Jacᵏ⁺¹ₖ,Jacᵏₖ,qₖ₊₁,qₖ,λₘ,Mₘ,∂Mₘhq̇ₘ∂qₘ,∂F∂q,∂F∂q̇,Aₖ₊₁,Aₖ,tₘ)
        ## yₖ .= transpose(Jacᵏₖ)\(-1/2*(∂ϕ∂xᵀ+∂ϕ∂xᵀ)-transpose(Jacᵏ⁺¹ₖ)*yₖ₊₁)
        yₖ .= transpose(Jacᵏₖ)\(-transpose(Jacᵏ⁺¹ₖ)*yₖ₊₁)
        #---------Step k finisher-----------
        step += 1
        #---------Step k finisher-----------
        if verbose
            dg_step = ceil(Int,log10(totalstep))+1
            dg_dt = max(1,-floor(Int,log10(dt)))
            wd_t = ceil(Int,log10(traj.t[end]))+dg_dt+1+1
            progfmt = Printf.Format("Progress: %5.1f%%, step: %$(dg_step)u, time: %$(wd_t).$(dg_dt)f, iterations: %s \n")
            progstr = Printf.format(progfmt,
                floor(timestep/totalstep*100;digits=1), timestep, tₖ, Res_stepk_result.iterations
            )
            print(progstr)
        end
        next!(prog)
    end

    bot
end
