struct Adjoint_Zhong06_Constant_Mass_Cache{CacheType}
    cache::CacheType
end

function generate_cache(
        simulator::Simulator{<:AbstractDynamicsProblem},
        solver::DiscreteAdjointDynamicsSolver{Zhong06},
        ::Val{true};
        dt,kargs...
    )
    (;prob) = simulator
    (;bot,policy) = prob
    (;structure) = bot
    options = merge(
        (gravity=true,factor=1,checkpersist=true), #default
        prob.options,
        solver.forward_solver.options,
    )
    Mₘ = assemble_M(structure)
    M⁻¹ₖ = assemble_M⁻¹(structure)
    ∂Mₘhq̇ₘ∂qₘ = assemble_∂Mq̇∂q(structure)
    M! = make_M!(structure)
    M⁻¹! = make_M⁻¹!(structure)
    M_and_Jac_M! = make_M_and_Jac_M!(structure)
    Φ = make_cstr_function(structure)
    A = make_cstr_jacobian(structure)
    F!(F,q,q̇,t) = generalized_force!(F,bot,policy,q,q̇,t;gravity=options.gravity)
    Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t) = generalized_force_jacobian!(∂F∂q̌,∂F∂q̌̇,bot,policy,q,q̇,t)
    ∂Aᵀλ∂q(q,λ) = cstr_forces_jacobian(structure,q,λ)
    adjoint_traj = VectorOfArray(
        [
            ComponentArray(
                a = zero(state.q),
                b = zero(state.p),
                μ = zero(state.λ)
            )
            for state in bot.traj
        ]
    )
    Jac_ϕ!(∂ϕ∂qᵀ,∂ϕ∂q̇ᵀ,∂ϕ∂pᵀ,q,q̇,t) = cost_jacobian!(∂ϕ∂qᵀ,∂ϕ∂q̇ᵀ,∂ϕ∂pᵀ,bot,q,q̇,t;gravity=options.gravity)
    cache = @eponymtuple(
        F!,Jac_F!,
        Mₘ,M⁻¹ₖ,
        ∂Mₘhq̇ₘ∂qₘ,
        M!,M⁻¹!,
        M_and_Jac_M!,
        Φ,A,∂Aᵀλ∂q,
        adjoint_traj,
        Jac_ϕ!
    )
    Adjoint_Zhong06_Constant_Mass_Cache(cache)
end

function solve!(simulator::Simulator,solvercache::Adjoint_Zhong06_Constant_Mass_Cache;
                dt,ftol=1e-14,verbose=false,maxiters=50,
                progress=true,exception=true)
    (;prob,controller,tspan,restart,totalstep) = simulator
    (;bot,) = prob
    (;structure,traj) = bot
    (;
        F!,Jac_F!,
        Mₘ,∂Mₘhq̇ₘ∂qₘ,
        M!,M_and_Jac_M!,
        Φ,A,∂Aᵀλ∂q,
        adjoint_traj,Jac_ϕ!
    ) = solvercache.cache
    T = get_numbertype(structure)
    nq = get_num_of_free_coords(structure)
    nλ = get_num_of_cstr(structure)
    nx = ny = 2nq+nλ
    ∂S∂xᵀ = zeros(T,nx)
    ∂S∂qᵀ = zeros(T,nq)
    ∂S∂q̇ᵀ = zeros(T,nq)
    ∂S∂pᵀ = zeros(T,nq)
    ∂F∂q = zeros(T,nq,nq)
    ∂F∂q̇ = zeros(T,nq,nq)
    Jacᵏₖ = Jacᵏ⁺¹ₖ₊₁  = zeros(T,ny,ny)
    Jacᵏₖ₋₁ = Jacᵏ⁺¹ₖ = zeros(T,ny,ny)
    Jacᵏ⁺¹ₖ_backup = zeros(T,ny,ny)
    q0 = traj.q[begin]
    M!(Mₘ,q0)
    step = 0
    mr = norm(Mₘ,Inf)
    mass_norm = mr

    function Jac!(Jacᵏ⁺¹ₖ₊₁,Jacᵏ⁺¹ₖ,qₖ₊₁,qₖ,λₘ,Mₘ,∂Mₘhq̇ₘ∂qₘ,∂F∂q,∂F∂q̇,Aₖ₊₁,Aₖ,tₘ)
        h = dt
        qₘ = (qₖ₊₁.+qₖ)./2
        q̇ₘ = (qₖ₊₁.-qₖ)./h
        M_and_Jac_M!(Mₘ,∂Mₘhq̇ₘ∂qₘ,qₘ,h.*q̇ₘ)
        Jac_F!(∂F∂q,∂F∂q̇,qₘ,q̇ₘ,tₘ)
        Jacᵏ⁺¹ₖ .= 0.0
        Jacᵏ⁺¹ₖ[   1:nq ,   1:nq ]     .=  1/2 .*∂Mₘhq̇ₘ∂qₘ .- Mₘ .-(h^2)/2 .*(1/2 .*∂F∂q .-1/h.*∂F∂q̇) .- mass_norm.*∂Aᵀλ∂q(qₖ,λₘ)
        Jacᵏ⁺¹ₖ[nq+1:2nq,   1:nq ]     .= -1/2 .*∂Mₘhq̇ₘ∂qₘ .+ Mₘ .-(h^2)/2 .*(1/2 .*∂F∂q .-1/h.*∂F∂q̇)
        Jacᵏ⁺¹ₖ[   1:nq ,nq+1:2nq]     .= -h*I(nq)
        
        Jacᵏ⁺¹ₖ₊₁ .= 0.0
        Jacᵏ⁺¹ₖ₊₁[    1:nq ,    1:nq ] .=  1/2 .*∂Mₘhq̇ₘ∂qₘ .+ Mₘ .-(h^2)/2 .*(1/2 .*∂F∂q .+1/h.*∂F∂q̇)
        Jacᵏ⁺¹ₖ₊₁[ nq+1:2nq,    1:nq ] .= -1/2 .*∂Mₘhq̇ₘ∂qₘ .- Mₘ .-(h^2)/2 .*(1/2 .*∂F∂q .+1/h.*∂F∂q̇) .- mass_norm.*∂Aᵀλ∂q(qₖ₊₁,λₘ)
        Jacᵏ⁺¹ₖ₊₁[2nq+1:end,    1:nq ] .=  mass_norm.*Aₖ₊₁
        Jacᵏ⁺¹ₖ₊₁[ nq+1:2nq, nq+1:2nq] .=  h*I(nq)
        Jacᵏ⁺¹ₖ₊₁[    1:nq ,2nq+1:end] .= -mass_norm.*transpose(Aₖ)
        Jacᵏ⁺¹ₖ₊₁[ nq+1:2nq,2nq+1:end] .= -mass_norm.*transpose(Aₖ₊₁)
    end
    qₖ₊₁ = traj.q[totalstep+1]
    q̇ₖ₊₁ = traj.q̇[totalstep+1]
    qₖ = traj.q[totalstep]
    λₘ = traj.λ[totalstep+1]
    tₖ₊₁ = traj.t[totalstep+1]
    tₖ = traj.t[totalstep]
    tₘ = 1/2*(tₖ₊₁ + tₖ)
    Aₖ₊₁ = A(qₖ₊₁)
    Aₖ = A(qₖ)
    Jac!(Jacᵏ⁺¹ₖ₊₁,Jacᵏ⁺¹ₖ,qₖ₊₁,qₖ,λₘ,Mₘ,∂Mₘhq̇ₘ∂qₘ,∂F∂q,∂F∂q̇,Aₖ₊₁,Aₖ,tₘ)
    yN = adjoint_traj[totalstep+1]
    ## yN .= transpose(Jacᵏ⁺¹ₖ₊₁)\(-1/2*∂ϕ∂xᵀ-∂S∂xᵀ)
    Jac_ϕ!(∂S∂qᵀ,∂S∂q̇ᵀ,∂S∂pᵀ,qₖ₊₁,q̇ₖ₊₁,tₖ₊₁)
    ∂S∂xᵀ[   1:nq] .= ∂S∂qᵀ
    ∂S∂xᵀ[nq+1:2nq] .= ∂S∂pᵀ
    yN .= transpose(Jacᵏ⁺¹ₖ₊₁)\(-∂S∂xᵀ)
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = totalstep:-1:2
        #---------Step k Control-----------
        # control!(sim,cache)
        #---------Step k Control-----------
        Jacᵏ⁺¹ₖ_backup .= Jacᵏ⁺¹ₖ
        tₖ = traj.t[timestep]
        tₖ₋₁ = traj.t[timestep-1]
        tₘ = 1/2*(tₖ + tₖ₋₁)
        yₖ₊₁ = adjoint_traj[timestep+1]
        yₖ = adjoint_traj[timestep]
        qₖ = traj.q[timestep]
        qₖ₋₁ = traj.q[timestep-1]
        λₘ = traj.λ[timestep]
        aₖ₊₁ = yₖ₊₁.a
        bₖ₊₁ = yₖ₊₁.b
        μₖ₊₁ = yₖ₊₁.μ
        aₖ  = yₖ.a
        bₖ  = yₖ.b
        μₖ  = yₖ.μ
        Aₖ = A(qₖ)
        Aₖ₋₁ = A(qₖ₋₁)
        Jac!(Jacᵏₖ,Jacᵏₖ₋₁,qₖ,qₖ₋₁,λₘ,Mₘ,∂Mₘhq̇ₘ∂qₘ,∂F∂q,∂F∂q̇,Aₖ,Aₖ₋₁,tₘ)
        ## yₖ .= transpose(Jacᵏₖ)\(-1/2*(∂ϕ∂xᵀ+∂ϕ∂xᵀ)-transpose(Jacᵏ⁺¹ₖ)*yₖ₊₁)
        yₖ .= transpose(Jacᵏₖ)\(-transpose(Jacᵏ⁺¹ₖ_backup)*yₖ₊₁)
        ## @show bₖ .- aₖ₊₁ |> norm
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

struct Adjoint_Sensitivity_Zhong06_Constant_Mass_Cache{CacheType}
    cache::CacheType
end

function generate_cache(
        simulator::Simulator{<:DynamicsSensitivityProblem},
        solver::AdjointDynamicsSensitivitySolver{
            <:DynamicsSolver{<:Zhong06},
            <:DiscreteAdjointDynamicsSolver{<:Zhong06}
        },
        ::Val{true};
        dt,kargs...
    )
    (;prob) = simulator
    (;bot,policy) = prob
    (;structure,hub) = bot
    options = merge(
        (gravity=true,factor=1,checkpersist=true), #default
        prob.options,
        solver.forward_solver.options,
    )
    forward_cache = solve!(simulator,solver.forward_solver;dt,kargs...)
    adjoint_cache = solve!(simulator,solver.adjoint_solver;dt,kargs...)
    ∂J∂θᵀ = [
        zero(policy.nt.ps)[:]
        for _ in  1:simulator.totalstep
    ]
    (;num_of_free_coords) = structure.connectivity.indexed
    T = get_numbertype(bot)
    nu = get_num_of_actions(bot)
    nθ = get_num_of_params(policy)
    nq = get_num_of_free_coords(structure)
    nλ = get_num_of_cstr(structure)
    nx = ny = 2nq+nλ
    ∂Fₘ∂θ = zeros(T,nq,nθ)
    Jacᵏ⁺¹ₘ = zeros(T,nx,nθ)
    function ∂f∂θ!(Jacᵏ⁺¹ₘ,∂Fₘ∂θ,q,q̇,t,h)
        generalized_force_jacobian!(∂Fₘ∂θ,structure,policy,q,q̇,t;gravity=options.gravity)
        Jacᵏ⁺¹ₘ[    1: nq,1:nθ] .= -h^2/2*∂Fₘ∂θ
        Jacᵏ⁺¹ₘ[ nq+1:2nq,1:nθ] .= -h^2/2*∂Fₘ∂θ
        Jacᵏ⁺¹ₘ[2nq+1:end,1:nθ] .= 0.0
    end
    cache = @eponymtuple(
        forward_cache,
        adjoint_cache,
        ∂J∂θᵀ,∂f∂θ!,∂Fₘ∂θ,Jacᵏ⁺¹ₘ
    )
    Adjoint_Sensitivity_Zhong06_Constant_Mass_Cache(cache)
end

function solve!(simulator::Simulator,solvercache::Adjoint_Sensitivity_Zhong06_Constant_Mass_Cache;
                dt,ftol=1e-14,verbose=false,maxiters=50,
                progress=true,exception=true)
    (;prob,controller,tspan,restart,totalstep) = simulator
    (;bot,) = prob
    (;structure,traj,control_traj) = bot
    (;
        forward_cache,
        adjoint_cache,
        ∂J∂θᵀ,∂f∂θ!,∂Fₘ∂θ,Jacᵏ⁺¹ₘ
    ) = solvercache.cache
    (;adjoint_traj) = adjoint_cache.cache
    h = dt
    T = get_numbertype(structure)
    nq = get_num_of_free_coords(structure)
    nλ = get_num_of_cstr(structure)
    nx = ny = 2nq+nλ
    ∂S∂xᵀ = zeros(T,nx)
    ∂S∂qᵀ = zeros(T,nq)
    ∂S∂q̇ᵀ = zeros(T,nq)
    ## q0 = traj.q[begin]
    ## M!(Mₘ,q0)
    ## step = 0
    ## mr = norm(Mₘ,Inf)
    ## mass_norm = mr
    for timestep = 1:totalstep
        qₖ₊₁ = traj.q[timestep+1]
        qₖ = traj.q[timestep]
        uₘ = control_traj.u[timestep+1]
        qₘ = (qₖ₊₁.+qₖ)./2
        q̇ₘ = (qₖ₊₁.-qₖ)./h
        tₖ₊₁ = traj.t[totalstep+1]
        tₖ = traj.t[totalstep]
        tₘ = 1/2*(tₖ₊₁ + tₖ)
        ∂f∂θ!(Jacᵏ⁺¹ₘ,∂Fₘ∂θ,qₘ,q̇ₘ,tₘ,h)
        ∂J∂θᵀ[timestep] = transpose(Jacᵏ⁺¹ₘ)*adjoint_traj[timestep+1]#+η*∂ϕ∂u
    end
end
