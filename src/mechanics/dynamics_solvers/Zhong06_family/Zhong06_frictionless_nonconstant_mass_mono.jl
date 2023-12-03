
struct Zhong06_Frictionless_Nonconstant_Mass_Mono_Cache{CacheType}
    cache::CacheType
end

function generate_cache(
        simulator::Simulator{<:DynamicsProblem{
            RobotType,
            EnvType,
            RestitutionFrictionCombined{NewtonRestitution,Frictionless}
        }},
        solver::DynamicsSolver{
            Zhong06,
            <:MonolithicContactSolver,
        },
        ::Val{false};
        dt,kargs...
    )   where {RobotType,EnvType}
    (;bot,env) = simulator.prob
    (;structure) = bot
    F!(F,q,q̇,t) = generalize_force!(F,bot,q,q̇,t)
    Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t) = generalize_force_jacobain!(∂F∂q̌,∂F∂q̌̇,bot,q,q̇,t)
    Mₘ = assemble_M(structure)
    ∂Mₘhq̇ₘ∂qₘ = assemble_∂Mq̇∂q(structure)
    M⁻¹ₖ = assemble_M⁻¹(structure)
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
    Zhong06_Frictionless_Nonconstant_Mass_Mono_Cache(cache)
end

function make_step_k(
        solver_cache::Zhong06_Frictionless_Nonconstant_Mass_Mono_Cache,
        nq,nλ,na,
        qₖ₋₁,vₖ₋₁,pₖ₋₁,tₖ₋₁,
        pₖ,vₖ,
        structure,
        contact_cache,
        h,scaling,
    )
    (;
        F!,Jac_F!,
        Mₘ,M⁻¹ₖ,
        ∂Mₘhq̇ₘ∂qₘ,
        ∂M⁻¹ₖpₖ∂qₖ,
        Fₘ,∂Fₘ∂qₘ,∂Fₘ∂q̇ₘ,
        M_and_Jac_M!,
        M⁻¹_and_Jac_M⁻¹!,
        Φ,A,
        ∂Aᵀλ∂q,
    ) = solver_cache.cache

    # T = eltype(qₖ₋₁)
    n1 = nq
    n2 = nq+nλ
    nx = n2+2na
    
    function ns_stepk!(𝐫𝐞𝐬,𝐉,𝐰,
            x,Λₘ,y,∂y∂x,
            Λ_split,y_split,
            Dₖ₋₁,ŕₖ₋₁,
            timestep,iteration)
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
                        scaling.*transpose(Aₖ₋₁)*λₘ
        𝐫𝐞𝐬[n1+1:n2] .= scaling.*Φ(qₖ)
        
        𝐉 .= 0.0
        𝐉[   1:n1,   1:n1] .=  Mₘ .+ 1/2 .*∂Mₘhq̇ₘ∂qₘ .-h^2/2 .*(1/2 .*∂Fₘ∂qₘ .+ 1/h.*∂Fₘ∂q̇ₘ)
        𝐉[   1:n1,n1+1:n2] .= -scaling.*transpose(Aₖ₋₁)
        𝐉[n1+1:n2,   1:n1] .=  scaling.*Aₖ
        
        if na != 0
            (;
                H,
                persistent_idx,
                restitution_coefficients                
            ) = contact_cache.cache
            𝐫𝐞𝐬[   1:n1] .-= scaling*h .*transpose(Dₖ₋₁)*H*Λₘ 
            get_directions_and_positions!(structure, contact_cache, qₖ, vₖ, H*Λₘ)
            Dₖ = contact_cache.cache.D
            ŕₖ = contact_cache.cache.ŕ
            ∂Dₖvₖ∂qₖ = contact_cache.cache.∂Dq̇∂q
            ∂DᵀₖHΛₘ∂qₖ = contact_cache.cache.∂DᵀΛ∂q
            pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,Λₘ,Dₖ₋₁,Dₖ,H,scaling,h)
            M⁻¹_and_Jac_M⁻¹!(M⁻¹ₖ,∂M⁻¹ₖpₖ∂qₖ,qₖ,pₖ)
            vₖ .= M⁻¹ₖ*pₖ
            ∂Aᵀₖλₘ∂qₖ = ∂Aᵀλ∂q(qₖ,λₘ)
            ∂pₖ∂qₖ = 2/h.*Mₘ + 
                     1/h.*∂Mₘhq̇ₘ∂qₘ .+
                    scaling/h.*∂Aᵀₖλₘ∂qₖ .+ 
                    scaling.*∂DᵀₖHΛₘ∂qₖ
            ∂pₖ∂λₘ = scaling/h.*transpose(Aₖ-Aₖ₋₁)
            ∂vₖ∂qₖ = M⁻¹ₖ*∂pₖ∂qₖ .+ ∂M⁻¹ₖpₖ∂qₖ
            ∂vₖ∂λₘ = M⁻¹ₖ*∂pₖ∂λₘ
            v́ₖ = Dₖ*vₖ
            ∂v́ₖ∂qₖ = Dₖ*∂vₖ∂qₖ .+ ∂Dₖvₖ∂qₖ 
            ∂v́ₘ∂qₖ = Dₖ./h 
            ∂y∂x .= 0
            v́ₖ₋₁ = Dₖ₋₁*vₖ₋₁
            v́ₘ = (ŕₖ .- ŕₖ₋₁)./h
            v́⁺ = copy(v́ₖ)
            for i = 1:na
                vⁱₖ₋₁ = v́ₖ₋₁[i]
                # vₜⁱₖ₋₁ = norm(vⁱₖ₋₁[2:3])
                if i in persistent_idx
                    vⁱ⁺   = v́ₘ[i]
                    vₙⁱₖ₋₁ = zero(vⁱₖ₋₁)
                else
                    vⁱ⁺   = v́ₖ[i]
                    vₙⁱₖ₋₁ = vⁱₖ₋₁
                end
                v́⁺[i] = vⁱ⁺
                # vₙⁱ   = vⁱ⁺
                # @show timestep,iteration, vₙⁱₖ₋₁, vₙⁱ, vₜⁱₖ₋₁, vₜⁱ, Λₘ
                𝐰[i] = restitution_coefficients[i]*min(vₙⁱₖ₋₁,zero(vₙⁱₖ₋₁))
                
                Dⁱₖ = @view Dₖ[[i],:]                
                if i in persistent_idx
                    ∂y∂x[i,   1:n1] .= ∂v́ₘ∂qₖ[i,:]                     
                    ∂y∂x[i,n1+1:n2] .= 0
                else
                    ∂y∂x[i,   1:n1] .= ∂v́ₖ∂qₖ[i,:]
                    # @show size(Dₖ), size(Dⁱₖ), size(∂vₖ∂λₘ)
                    ∂y∂x[[i],n1+1:n2] .= Dⁱₖ*∂vₖ∂λₘ
                end
            end
            𝐫𝐞𝐬[(n2   +1):(n2+ na)] .= (h.*(v́⁺ .+ 𝐰) .- h.*y)
            𝐫𝐞𝐬[n2+na+1:n2+2na]     .= reduce(vcat,Λ_split⊙y_split)

            𝐉[      1:n1    , n2+   1:n2+ na] .=  -scaling*h .*transpose(Dₖ₋₁)*H
            𝐉[n2+1:n2+ na,      1:n2    ]     .=  h.*∂y∂x
            𝐉[n2+1:n2+ na,    n2+na+1:n2+2na] .= -h.*I(na)
            𝐉[n2+na+1:n2+2na, n2+   1:n2+ na] .=  BlockDiagonal(mat.(y_split))
            𝐉[n2+na+1:n2+2na, n2+na+1:n2+2na] .=  BlockDiagonal(mat.(Λ_split))
            # @show cond(𝐉)
        end
    end
    ns_stepk!
end

function solve!(sim::Simulator,solver_cache::Zhong06_Frictionless_Nonconstant_Mass_Mono_Cache;
        dt,
        ftol=1e-14,xtol=ftol,maxiters=50,
        verbose=false, verbose_contact=false,
        progress=true,
        exception=true,
    )
    (;prob,totalstep) = sim
    (;bot,env) = prob
    (;structure,traj,contacts_traj) = bot
    (;Mₘ,M⁻¹ₖ,M!,M⁻¹!,A,contacts_bits) = solver_cache.cache
    q0 = traj.q[begin]
    λ0 = traj.λ[begin]
    q̇0 = traj.q̇[begin]
    M!(Mₘ,q0)
    pₖ₋₁ = Mₘ*q̇0
    pₖ   = deepcopy(pₖ₋₁)
    qₖ₋½ = deepcopy(q0)
    T = eltype(q0)
    nq = length(q0)
    nλ = length(λ0)
    activate_contacts!(structure,env,solver_cache,q0)
    mr = norm(Mₘ,Inf)
    scaling = mr
    iteration = 0
    α_record = ones(T,maxiters)
    condition_number = typemax(T)
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
        contact_cache = activate_contacts!(structure,env,solver_cache,qₖ₋½)
        (;na,) = contact_cache.cache
        isconverged = false
        normRes = typemax(T)
        iteration_break = 0
        n1 = nq
        n2 = n1 + nλ
        nx = n2 + 2na
        Δx = zeros(T,nx)
        x = zero(Δx)
        Res = zero(Δx)
        Jac = zeros(T,nx,nx)
        Λₘ = @view x[(n2+1):n2+na]
        y  = @view x[n2+na+1:n2+2na]
        𝐰 = zeros(T,na)
        ∂y∂x = zeros(T,na,n2)
        𝐞 = ones(T,na)
        Λ_split = split_by_lengths(Λₘ,1)
        y_split = split_by_lengths(y,1)
        Λp = zero(Λₘ)
        yp = zero(y)
        Δxp = zeros(T,nx)
        ΔΛp = @view Δxp[(n2+1):n2+na]
        Δyp = @view Δxp[n2+na+1:n2+2na]
        ΔΛp_split = split_by_lengths(ΔΛp,1)
        Δyp_split = split_by_lengths(Δyp,1)
        Δxc = zeros(T,nx)
        ΔΛc = @view Δxc[(n2+1):n2+na]
        Δyc = @view Δxc[n2+na+1:n2+2na]
        ΔΛc_split = split_by_lengths(ΔΛc,1)
        Δyc_split = split_by_lengths(Δyc,1)
        get_directions_and_positions!(structure, contact_cache, qₖ₋₁, q̇ₖ₋₁, Λₘ)
        (;
            H
        ) = contact_cache.cache
        Dₖ₋₁ = deepcopy(contact_cache.cache.D)
        ŕₖ₋₁ = deepcopy(contact_cache.cache.ŕ)
        ns_stepk! = make_step_k(
            solver_cache,
            nq,nλ,na,qₖ₋₁,q̇ₖ₋₁,pₖ₋₁,tₖ₋₁,pₖ,q̇ₖ,
            structure,contact_cache,
            dt,scaling,
        )
        restart_count = 0
        Λ_guess = 0.1
        while restart_count < 10
            Λₘ .= repeat([Λ_guess],na)
            y .= Λₘ
            x[      1:nq]          .= qₖ
            x[   nq+1:nq+nλ]       .= 0.0
            Nmax = 50
            α_record = fill(typemax(T),maxiters)
            for iteration = 1:maxiters
                # @show iteration,D,ηs,restitution_coefficients,gaps
                ns_stepk!(
                    Res,Jac,
                    𝐰,x,Λₘ,y,∂y∂x,
                    Λ_split,y_split,
                    Dₖ₋₁,ŕₖ₋₁,
                    timestep,iteration
                )
                condition_number = cond(Jac)
                if na == 0
                    normRes = norm(Res)
                    if normRes < ftol
                        isconverged = true
                        iteration_break = iteration-1
                        break
                    end
                    Δx .= Jac\(-Res)
                    x .+= Δx
                else # na!=0
                    # @show timestep,iteration,normRes,Λₘ
                    # Λₘini = repeat([Λ_guess,0,0],na)
                    μ = transpose(y)*Λₘ/na
                    lu𝐉 = lu(Jac)
                    Δxp .= lu𝐉\(-Res)
                    αp_Λ = find_nonnegative_step_length(Λ_split,ΔΛp_split)
                    αp_y = find_nonnegative_step_length(y_split,Δyp_split)
                    αpmax = min(αp_Λ,αp_y)
                    # αpmax = find_nonnegative_step_length(z_split,W_blocks,Δyp_split,ΔΛp_split,J)
                    αp = min(one(αpmax),0.99αpmax)
                    Λp .= Λₘ .+ αp.*ΔΛp
                    yp .= y .+ αp.*Δyp
                    μp = transpose(yp)*Λp/na
                    σ = (μp/μ)^3
                    if σ == NaN || μ == 0
                        break
                    end
                    τ = σ*μp
                    Res_c = -τ.*𝐞.+(Δyp.*ΔΛp)
                    Res[n2+na+1:n2+2na] .+= Res_c
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
                    Δxc .= lu𝐉\(-Res)
                    # η = exp(-0.1μ) + 0.9
                    α_Λ = find_nonnegative_step_length(Λ_split,ΔΛc_split)
                    # @show Λ_split,ΔΛc_split
                    α_y = find_nonnegative_step_length(y_split,Δyc_split)
                    αmax = min(α_Λ,α_y)
                    α = min(1,0.999αmax)
                    α_record[iteration] = α
                    # Λₘ .+= α.*ΔΛc
                    # y .+= α.*Δyc

                    # @show Λₘ, y
                    x .+= α.*Δxc
                    μ = transpose(y)*Λₘ/na
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
        get_directions_and_positions!(structure, contact_cache, qₖ, q̇ₖ, Λₘ,)
        Dₖ = contact_cache.cache.D
        pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,Λₘ,Dₖ₋₁,Dₖ,H,scaling,dt)
        M⁻¹!(M⁻¹ₖ,qₖ)
        q̇ₖ .= M⁻¹ₖ*pₖ
        if na != 0
            update_contacts!(cₖ[contacts_bits],cₖ₋₁[contacts_bits],Dₖ*q̇ₖ,Λₘ./(scaling*dt))
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
        record!(
            sim.solver_history,
            (
                residual=normRes,
                iteration=iteration_break,
                walltime = 1.0,
                num_of_contacts = na,
                stepsizes = α_record,
                condition_number = condition_number
            ),
            timestep
        )
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
