
struct Zhong06_CCP_Nonconstant_Mass_Cache{CacheType}
    cache::CacheType
end

function generate_cache(
        simulator::Simulator{<:DynamicsProblem{
            RobotType,
            EnvType,
            RestitutionFrictionCombined{NewtonRestitution,CoulombFriction}
        }},
        solver::DynamicsSolver{
            Zhong06,
            <:InnerLayerContactSolver
        },
        ::Val{false};
        dt,kargs...
    )   where {RobotType,EnvType}
    (;bot,env) = simulator.prob
    (;structure) = bot
    F!(F,q,q̇,t) = generalized_force!(F,bot,q,q̇,t)
    Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t) = generalized_force_jacobain!(∂F∂q̌,∂F∂q̌̇,bot,q,q̇,t)
    Mₘ = assemble_M(structure) 
    M⁻¹ₖ = assemble_M⁻¹(structure)
    ∂Mₘhq̇ₘ∂qₘ = assemble_∂Mq̇∂q(structure)
    ∂M⁻¹ₖpₖ∂qₖ = assemble_∂M⁻¹p∂q(structure)
    M! = make_M!(structure)
    M⁻¹! = make_M⁻¹!(structure)
    M_and_Jac_M! = make_M_and_Jac_M!(structure)
    M⁻¹_and_Jac_M⁻¹! = make_M⁻¹_and_Jac_M⁻¹!(structure)
    Φ = make_cstr_function(structure)
    A = make_cstr_jacobian(structure)

    nq = size(Mₘ,2)
    T = get_numbertype(structure)
    Fₘ = zeros(T,nq)
    ∂Fₘ∂qₘ = zeros(T,nq,nq)
    ∂Fₘ∂q̇ₘ = zeros(T,nq,nq)
    Ψ(q,q̇) = Vector{T}()
    ∂Ψ∂q(q,q̇) = Matrix{T}(undef,0,nq)
    B(q) = Matrix{T}(undef,0,nq)
    # ∂𝐌𝐚∂𝐪(q,a) = zeros(T,nq,nq)
    ∂Aᵀλ∂q(q,λ) = cstr_forces_jacobian(structure,q,λ)
    # ∂𝚽𝐪𝐯∂𝒒(q,v) = RB.∂Aq̇∂q(st,v)
    ∂Bᵀμ∂q(q,μ) = zeros(T,nq,nq)
    (;
        contacts_bits,
        persistent_bits,
        μs_sys,
        es_sys,
        gaps_sys
    ) = prepare_contacts(bot,env)
    cache = @eponymtuple(
        F!,Jac_F!,
        Mₘ,M⁻¹ₖ,
        ∂Mₘhq̇ₘ∂qₘ,
        ∂M⁻¹ₖpₖ∂qₖ,
        Fₘ,∂Fₘ∂qₘ,∂Fₘ∂q̇ₘ,
        M!,M⁻¹!,
        M_and_Jac_M!,
        M⁻¹_and_Jac_M⁻¹!,
        Φ,A,Ψ,B,
        ∂Ψ∂q,
        ∂Aᵀλ∂q,
        ∂Bᵀμ∂q,
        contacts_bits,
        persistent_bits,
        μs_sys,
        es_sys,
        gaps_sys
    )
    Zhong06_CCP_Nonconstant_Mass_Cache(cache)
end

function make_step_k(
        solver_cache::Zhong06_CCP_Nonconstant_Mass_Cache,
        nq,nλ,na,
        qₖ₋₁,vₖ₋₁,pₖ₋₁,tₖ₋₁,
        pₖ,vₖ,
        structure,
        contact_cache,
        h,mass_norm,
    )
    (;
        F!,Jac_F!,
        Mₘ,M⁻¹ₖ,
        ∂Mₘhq̇ₘ∂qₘ,
        ∂M⁻¹ₖpₖ∂qₖ,
        Fₘ,∂Fₘ∂qₘ,∂Fₘ∂q̇ₘ,
        M!,M⁻¹!,
        M_and_Jac_M!,
        M⁻¹_and_Jac_M⁻¹!,
        Φ,A,
        ∂Aᵀλ∂q,
    ) = solver_cache.cache
    # T = eltype(qₖ₋₁)
    n1 = nq
    n2 = nq+nλ
    nΛ = 3na
    nx = n2
    function ns_stepk!(
            𝐫𝐞𝐬,𝐉,
            𝐁,𝐛,𝐜ᵀ,𝐍,𝐫,
            x,Λₘ,
            Dₖ₋₁,ŕₖ₋₁,
            timestep,iteration
        )
        # @show timestep, iteration, na, persistent_idx
        qₖ = @view x[   1:n1]
        λₘ = @view x[n1+1:n2]
        qₘ = (qₖ.+qₖ₋₁)./2
        q̇ₘ = (qₖ.-qₖ₋₁)./h
        vₘ = q̇ₘ
        tₘ = tₖ₋₁+h/2
        M_and_Jac_M!(Mₘ,∂Mₘhq̇ₘ∂qₘ,qₘ,h.*q̇ₘ)
        F!(Fₘ,qₘ,q̇ₘ,tₘ)
        Jac_F!(∂Fₘ∂qₘ,∂Fₘ∂q̇ₘ,qₘ,q̇ₘ,tₘ)
        Aₖ₋₁ = A(qₖ₋₁)
        Aₖ   = A(qₖ)

        𝐫𝐞𝐬[   1:n1] .= h.*Mₘ*vₘ .- 
                        h.*pₖ₋₁ .-
                        (h^2)/2 .*Fₘ .-
                        mass_norm.*transpose(Aₖ₋₁)*λₘ 
        𝐫𝐞𝐬[n1+1:n2] .= mass_norm.*Φ(qₖ)
        
        𝐉 .= 0.0
        𝐉[   1:n1,   1:n1] .=  Mₘ .+ 1/2 .*∂Mₘhq̇ₘ∂qₘ .-h^2/2 .*(1/2 .*∂Fₘ∂qₘ .+ 1/h.*∂Fₘ∂q̇ₘ)
        𝐉[   1:n1,n1+1:n2] .= -mass_norm.*transpose(Aₖ₋₁)
        𝐉[n1+1:n2,   1:n1] .=  mass_norm.*Aₖ
        
        if na != 0
            (;
                H,
                restitution_coefficients,
                persistent_idx
            ) = contact_cache.cache
            𝐫𝐞𝐬[   1:n1] .-= mass_norm*h .*transpose(Dₖ₋₁)*H*Λₘ 
            # get_directions_and_positions!(Dₖ,Dper, Dimp, ∂Dₖvₖ∂qₖ, ∂DᵀₖHΛₘ∂qₖ,ŕₖ,qₖ, vₖ, H*Λₘ,bodyid2act_idx)
            get_frictional_directions_and_positions!(structure, contact_cache, qₖ, vₖ, H*Λₘ)
            Dₖ = contact_cache.cache.D
            ŕₖ = contact_cache.cache.ŕ
            ∂Dₖvₖ∂qₖ = contact_cache.cache.∂Dq̇∂q
            ∂DᵀₖHΛₘ∂qₖ = contact_cache.cache.∂DᵀΛ∂q
            pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,Λₘ,Dₖ₋₁,Dₖ,H,mass_norm,h)
            M⁻¹_and_Jac_M⁻¹!(M⁻¹ₖ,∂M⁻¹ₖpₖ∂qₖ,qₖ,pₖ)
            vₖ .= M⁻¹ₖ*pₖ
            ∂Aᵀₖλₘ∂qₖ = ∂Aᵀλ∂q(qₖ,λₘ)
            ∂pₖ∂qₖ = 2/h.*Mₘ + 
                     1/h.*∂Mₘhq̇ₘ∂qₘ .+
                     mass_norm/h.*∂Aᵀₖλₘ∂qₖ .+ 
                     mass_norm.*∂DᵀₖHΛₘ∂qₖ
            ∂pₖ∂λₘ = mass_norm/h.*transpose(Aₖ-Aₖ₋₁)
            ∂vₖ∂qₖ = M⁻¹ₖ*∂pₖ∂qₖ .+ ∂M⁻¹ₖpₖ∂qₖ
            ∂vₖ∂λₘ = M⁻¹ₖ*∂pₖ∂λₘ
            𝐁 .= 0
            𝐁[  1:n1,   1:nΛ] .= mass_norm.*h .*transpose(Dₖ₋₁)*H
            v́ₖ = Dₖ*vₖ
            ∂v́ₖ∂qₖ = Dₖ*∂vₖ∂qₖ .+ ∂Dₖvₖ∂qₖ 
            ∂v́ₘ∂qₖ = Dₖ./h 
            𝐜ᵀ .= 0
            v́ₖ₋₁ = Dₖ₋₁*vₖ₋₁
            v́ₘ = (ŕₖ .- ŕₖ₋₁)./h
            v́⁺ = copy(v́ₖ)
            for i = 1:na
                is = 3(i-1)
                vⁱₖ₋₁ = @view v́ₖ₋₁[is+1:is+3]
                # vₜⁱₖ₋₁ = norm(vⁱₖ₋₁[2:3])
                if i in persistent_idx
                    vⁱ⁺   = @view v́ₘ[is+1:is+3]
                    vₙⁱₖ₋₁ = zero(vⁱₖ₋₁[1])
                else
                    vⁱ⁺   = @view v́ₖ[is+1:is+3]
                    vₙⁱₖ₋₁ = vⁱₖ₋₁[1]
                end
                v́⁺[is+1:is+3] = vⁱ⁺
                vₜⁱ⁺   = norm(vⁱ⁺[2:3])
                # vₙⁱ   = vⁱ⁺[1]
                # @show timestep,iteration, vₙⁱₖ₋₁, vₙⁱ, vₜⁱₖ₋₁, vₜⁱ, Λₘ
                v́ₜⁱ = vₜⁱ⁺ + restitution_coefficients[i]*min(vₙⁱₖ₋₁,zero(vₙⁱₖ₋₁))
                𝐛[is+1:is+3] .= [v́ₜⁱ,0,0]
                
                Dⁱₖ = @view Dₖ[is+1:is+3,:]                
                if i in persistent_idx
                    𝐜ᵀ[is+1:is+3,   1:n1] .= ∂v́ₘ∂qₖ[is+1:is+3,:]                     
                    𝐜ᵀ[is+1     ,   1:n1] .+= 1/(norm(v́ₘ[is+2:is+3])+1e-14)*(v́ₘ[is+2]*∂v́ₘ∂qₖ[is+2,:] .+ v́ₘ[is+3]*∂v́ₘ∂qₖ[is+3,:])
                    𝐜ᵀ[is+1:is+3,n1+1:n2] .= 0
                else
                    𝐜ᵀ[is+1:is+3,   1:n1] .= ∂v́ₖ∂qₖ[is+1:is+3,:]
                    𝐜ᵀ[is+1     ,   1:n1] .+= 1/(norm(v́ₖ[is+2:is+3])+1e-14)*(v́ₖ[is+2]*∂v́ₖ∂qₖ[is+2,:] .+ v́ₖ[is+3]*∂v́ₖ∂qₖ[is+3,:])
                    𝐜ᵀ[is+1:is+3,n1+1:n2] .= Dⁱₖ*∂vₖ∂λₘ
                end
            end

            𝐜ᵀinv𝐉 = 𝐜ᵀ*inv(𝐉)
            𝐍 .= 𝐜ᵀinv𝐉*𝐁
            # debug
            # @show norm(D*vₖ + 𝐛), norm(𝐫𝐞𝐬)
            # @show Λₘ, D*vₖ, 𝐛
            # @show v́ₖ
            # @show Λₘ[1:3]⋅(v́ₖ + 𝐛)[1:3]
            𝐫 .= (v́⁺ + 𝐛) - 𝐜ᵀinv𝐉*(𝐫𝐞𝐬 + 𝐁*Λₘ)
        end

    end
    ns_stepk!
end

function solve!(sim::Simulator,solver_cache::Zhong06_CCP_Nonconstant_Mass_Cache;
        dt,
        ftol=1e-14,xtol=ftol,maxiters=50,
        verbose=false, verbose_contact=false,
        progress=true,
        exception=true,
    )
    (;prob,totalstep) = sim
    (;bot,env,) = prob
    (;structure,traj,contacts_traj) = bot
    (;Mₘ,M⁻¹ₖ,M!,M⁻¹!,A,contacts_bits) = solver_cache.cache
    q0 = traj.q[begin]
    λ0 = traj.λ[begin]
    q̇0 = traj.q̇[begin]
    activate_contacts!(structure,env,solver_cache,q0)
    M!(Mₘ,q0)
    pₖ₋₁ = Mₘ*q̇0
    pₖ   = deepcopy(pₖ₋₁)
    qₖ₋½ = deepcopy(q0)
    T = eltype(q0)
    nq = length(q0)
    nλ = length(λ0)
    nx = nq + nλ
    Δx = zeros(T,nx)
    x = zero(Δx)
    Res = zero(Δx)
    Jac = zeros(T,nx,nx)
    mr = norm(Mₘ,Inf)
    mass_norm = mr
    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Time Step k Control-----------
        # control!(sim,cache)
        #---------Time Step k Control-----------
        cₖ₋₁ = contacts_traj[timestep]
        cₖ = contacts_traj[timestep+1]
        qₖ₋₁ = traj.q[timestep]
        q̇ₖ₋₁ = traj.q̇[timestep]
        tₖ₋₁ = traj.t[timestep]
        qₖ   = traj.q[timestep+1]
        q̇ₖ   = traj.q̇[timestep+1]
        λₘ   = traj.λ[timestep+1]
        qₖ₋½ .= qₖ₋₁ .+ dt./2 .*q̇ₖ₋₁
        qₖ .= qₖ₋₁ .+ dt .*q̇ₖ₋₁
        q̇ₖ .= q̇ₖ₋₁
        # na,bodyid2act_idx,persistent_idx,contacts_bits,
        # H,restitution_coefficients,Dₖ₋₁, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕₖ₋₁, 
        # L = prepare_contacts!(qₖ₋½)
        contact_cache = activate_frictional_contacts!(structure,env,solver_cache,qₖ₋½)
        (;na,) = contact_cache.cache
        isconverged = false
        normRes = typemax(T)
        iteration_break = 0
        nΛ = 3na
        Λₘ = zeros(T,nΛ)
        Λʳₖ = copy(Λₘ)
        ΔΛₖ = copy(Λₘ)
        𝐁 = zeros(T,nx,nΛ)
        𝐛 = zeros(T,nΛ)
        𝐜ᵀ = zeros(T,nΛ,nx)
        𝐍 = zeros(T,nΛ,nΛ)
        𝐫 = zeros(T,nΛ)
        # get_directions_and_positions!(Dₖ₋₁, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕₖ₋₁, qₖ₋₁, q̇ₖ₋₁, Λₘ, bodyid2act_idx,)
        get_frictional_directions_and_positions!(structure, contact_cache, qₖ₋₁, q̇ₖ₋₁, Λₘ)
        (;
            H
        ) = contact_cache.cache
        Dₖ₋₁ = deepcopy(contact_cache.cache.D)
        ŕₖ₋₁ = deepcopy(contact_cache.cache.ŕ)
        ns_stepk! = make_step_k(
            solver_cache,
            nq,nλ,na,
            qₖ₋₁,q̇ₖ₋₁,pₖ₋₁,tₖ₋₁,
            pₖ,q̇ₖ,
            structure,
            contact_cache,
            dt,mass_norm,
        )
        restart_count = 0
        Λ_guess = 1.0
        while restart_count < 10
            Λₘ .= repeat([Λ_guess,0,0],na)
            x[      1:nq]          .= qₖ
            x[   nq+1:nq+nλ]       .= 0.0
            Λʳₖ .= Λₘ
            Nmax = 50
            for iteration = 1:maxiters
                # @show iteration,D,ηs,restitution_coefficients,gaps
                ns_stepk!(
                    Res,Jac,
                    𝐁,𝐛,𝐜ᵀ,𝐍,𝐫,
                    x,Λₘ,
                    Dₖ₋₁,ŕₖ₋₁,
                    timestep,iteration
                )
                normRes = norm(Res)
                if  normRes < ftol
                    isconverged = true
                    iteration_break = iteration-1
                    break
                elseif normRes > 1e10
                    # force restart
                    iteration_break = iteration-1
                    isconverged = false
                    break
                elseif iteration == maxiters
                    iteration_break = iteration-1
                    isconverged = false
                end
                if na == 0
                    Δx .= Jac\(-Res)
                    x .+= Δx
                else # na!=0
                    # @show timestep,iteration,normRes,Λₘ
                    # Λₘini = repeat([Λ_guess,0,0],na)
                    get_distribution_law!(structure,contact_cache,x[1:nq])
                    (;L) = contact_cache.cache
                    Λₘini = deepcopy(Λₘ)
                    Λₘini[begin+1:3:end] .= 0.0
                    Λₘini[begin+2:3:end] .= 0.0
                    # yini = deepcopy(Λₘini)
                    𝐍 .+= L
                    yini = 𝐍*Λₘ + 𝐫
                    yini .= abs.(yini)
                    yini[begin+1:3:end] .= 0.0
                    yini[begin+2:3:end] .= 0.0
                    IPM!(Λₘ,na,nΛ,Λₘini,yini,𝐍,𝐫;ftol,Nmax)
                    ΔΛₖ .= Λₘ - Λʳₖ
                    minusResΛ = -Res + 𝐁*(ΔΛₖ)
                    normRes = norm(minusResΛ)
                    Δx .= Jac\minusResΛ
                    Λʳₖ .= Λₘ
                    x .+= Δx
                    # @show timestep, iteration, normRes, norm(Δx), norm(ΔΛₖ),persistent_idx
                end
            end
            if isconverged
                break
            end
            restart_count += 1
            Λ_guess /= 10
            # @warn "restarting step: $timestep, count: $restart_count, Λ_guess = $Λ_guess"
        end
        qₖ .= x[      1:nq]
        λₘ .= x[   nq+1:nq+nλ]
        qₖ₋½ .= (qₖ.+qₖ₋₁)./2
        M!(Mₘ,qₖ₋½)
        # get_directions_and_positions!(Dₖ, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕₖ, qₖ, q̇ₖ, Λₘ, bodyid2act_idx)
        get_frictional_directions_and_positions!(structure, contact_cache, qₖ, q̇ₖ, Λₘ,)
        Dₖ = contact_cache.cache.D
        pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,Λₘ,Dₖ₋₁,Dₖ,H,mass_norm,dt)
        M⁻¹!(M⁻¹ₖ,qₖ)
        q̇ₖ .= M⁻¹ₖ*pₖ
        if na != 0
            update_contacts!(cₖ[contacts_bits],cₖ₋₁[contacts_bits],Dₖ*q̇ₖ,Λₘ./(mass_norm*dt))
        end
        if !isconverged
            @warn "Newton max iterations $maxiters, at timestep=$timestep, normRes=$(normRes)"
            if exception
                @error "Not converged!"
                break
            else
                # sim.convergence = false
                # break
            end
        end

        #---------Time Step k finisher-----------
        pₖ₋₁, pₖ = pₖ, pₖ₋₁
        if verbose || (na>0 && verbose_contact)
            dg_step = ceil(Int,log10(totalstep))+1
            dg_dt = max(1,-floor(Int,log10(dt)))
            wd_t = ceil(Int,log10(traj.t[end]))+dg_dt+1+1
            progfmt = Printf.Format("Prog.: %5.1f%%, step: %$(dg_step)u, time: %$(wd_t).$(dg_dt)f, iters: %s, contacts: %s \n")
            progstr = Printf.format(progfmt,
                floor(timestep/totalstep*100;digits=1), timestep, traj.t[timestep], iteration_break, na
            )
            print(progstr)
        end
        next!(prog)
    end
    bot
end
