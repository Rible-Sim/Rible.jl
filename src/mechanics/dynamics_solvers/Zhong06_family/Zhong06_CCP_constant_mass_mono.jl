
struct Zhong06_CCP_Constant_Mass_Mono_Cache{CacheType}
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
            <:MonolithicContactSolver
        },
        ::Val{true};
        dt,kargs...
    )   where {RobotType,EnvType}
    (;prob) = simulator
    (;bot,env) = prob
    (;structure) = bot
    F!(F,q,q̇,t) = generalized_force!(F,bot,q,q̇,t;gravity=true)
    Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t) = generalized_force_jacobain!(∂F∂q̌,∂F∂q̌̇,bot,q,q̇,t)
    
    M = Matrix(assemble_M(structure))
    Φ = make_cstr_function(bot)
    A = make_cstr_jacobian(bot)

    nq = size(M,2)
    T = get_numbertype(bot)
    Ψ(q,q̇) = Vector{T}()
    ∂Ψ∂q(q,q̇) = Matrix{T}(undef,0,nq)
    B(q) = Matrix{T}(undef,0,nq)

    # ∂𝐌𝐚∂𝐪(q,a) = zeros(T,nq,nq)
    ∂Aᵀλ∂q(q::AbstractVector,λ) = cstr_forces_jacobian(structure,q,λ)
    # ∂𝚽𝐪𝐯∂𝒒(q,v) = RB.∂Aq̇∂q(structure,v)
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
        M,Φ,A,Ψ,B,∂Ψ∂q,∂Aᵀλ∂q,∂Bᵀμ∂q,
        contacts_bits,
        persistent_bits,
        μs_sys,
        es_sys,
        gaps_sys
    )
    Zhong06_CCP_Constant_Mass_Mono_Cache(cache)
end

function make_step_k(
        solver_cache::Zhong06_CCP_Constant_Mass_Mono_Cache,
        nq,nλ,na,
        qₖ₋₁,vₖ₋₁,pₖ₋₁,tₖ₋₁,
        pₖ,vₖ,
        invM,
        h,mass_norm)
    (;F!,Jac_F!,M,Φ,A,∂Aᵀλ∂q) = solver_cache.cache

    n1 = nq
    n2 = nq+nλ
    nΛ = 3na
    nx = n2
    function ns_stepk!(
            𝐫𝐞𝐬,𝐉,
            Fₘ,∂F∂q,∂F∂q̇,
            𝐰,x,Λₖ,y,∂y∂x,
            Λ_split,y_split,
            structure,
            contact_cache,
            timestep,iteration
        )
        # @show timestep, iteration, na
        qₖ = @view x[   1:n1]
        λₘ = @view x[n1+1:n2]
        qₘ = (qₖ.+qₖ₋₁)./2
        q̇ₘ = (qₖ.-qₖ₋₁)./h
        vₘ = q̇ₘ
        tₘ = tₖ₋₁+h/2
        F!(Fₘ,qₘ,q̇ₘ,tₘ)
        Jac_F!(∂F∂q,∂F∂q̇,qₘ,q̇ₘ,tₘ)

        Aₖ₋₁ = A(qₖ₋₁)
        Aₖ   = A(qₖ)

        𝐫𝐞𝐬[   1:n1] .= -h.*pₖ₋₁ .+ M*(qₖ.-qₖ₋₁) .-
                        mass_norm.*transpose(Aₖ₋₁)*λₘ .-
                        (h^2)/2 .*Fₘ
        𝐫𝐞𝐬[n1+1:n2] .= -mass_norm.*Φ(qₖ)
        
        𝐉 .= 0.0
        𝐉[   1:n1,   1:n1] .=  M .-h^2/2 .*(1/2 .*∂F∂q .+ 1/h.*∂F∂q̇)
        𝐉[   1:n1,n1+1:n2] .= -mass_norm.*transpose(Aₖ₋₁)

        𝐉[n1+1:n2,   1:n1] .=  -mass_norm.*Aₖ

        if na != 0
            (;
                H,
                restitution_coefficients,
                D,
            ) = contact_cache.cache
            Dₘ = contact_cache.cache.Dper
            Dₖ = contact_cache.cache.Dimp
            𝐫𝐞𝐬[   1:n1]  .-= h.*mass_norm.*transpose(D)*H*Λₖ 

            pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,M,A,mass_norm,h)
            vₖ .= invM*pₖ        
            ∂vₘ∂qₖ = 1/h*I
            ∂vₖ∂qₖ = 2/h*I  + mass_norm/(h).*invM*(∂Aᵀλ∂q(qₖ,λₘ))
            ∂vₖ∂λₘ = mass_norm.*invM*transpose(Aₖ-Aₖ₋₁)/(h)
            
            v́⁺ = Dₘ*vₘ .+ Dₖ*vₖ
            ∂v́⁺∂qₖ = Dₘ*∂vₘ∂qₖ .+ Dₖ*∂vₖ∂qₖ
            ∂y∂x .= 0
            v́ₖ₋₁ = Dₖ*vₖ₋₁
            for i = 1:na
                is = 3(i-1)
                vⁱₖ₋₁ = @view v́ₖ₋₁[is+1:is+3]
                vⁱ⁺   = @view v́⁺[is+1:is+3]
                vₜⁱₖ₋₁ = norm(vⁱₖ₋₁[2:3])
                vₜⁱ⁺   = norm(vⁱ⁺[2:3])
                vₙⁱₖ₋₁ = vⁱₖ₋₁[1]
                vₙⁱ   = vⁱ⁺[1]
                v́ₜⁱ = vₜⁱ⁺ + restitution_coefficients[i]*min(vₙⁱₖ₋₁,0)
                𝐰[is+1:is+3] .= [v́ₜⁱ,0,0]
                Dⁱₘ = @view Dₘ[is+1:is+3,:]
                Dⁱₖ = @view Dₖ[is+1:is+3,:]
                ∂y∂x[is+1     ,   1:n1] .= 1/(norm(v́⁺[is+2:is+3])+1e-14)*(v́⁺[is+2]*∂v́⁺∂qₖ[is+2,:] .+ v́⁺[is+3]*∂v́⁺∂qₖ[is+3,:])
                ∂y∂x[is+1:is+3,   1:n1] .+= ∂v́⁺∂qₖ[is+1:is+3,:]
                ∂y∂x[is+1:is+3,n1+1:n2] .= Dⁱₖ*∂vₖ∂λₘ
            end
            𝐫𝐞𝐬[(n2   +1):(n2+ nΛ)] .= (h.*(v́⁺ .+ 𝐰) .- h.*y)
            𝐫𝐞𝐬[n2+nΛ+1:n2+2nΛ]     .= reduce(vcat,Λ_split⊙y_split)
            𝐉[      1:n1    , n2+   1:n2+ nΛ] .=  -mass_norm*h .*transpose(D)*H
            𝐉[n2+1:n2+ nΛ,      1:n2    ]     .=  h.*∂y∂x
            𝐉[n2+1:n2+ nΛ,    n2+nΛ+1:n2+2nΛ] .= -h.*I(nΛ)
            𝐉[n2+nΛ+1:n2+2nΛ, n2+   1:n2+ nΛ] .=  BlockDiagonal(mat.(y_split))
            𝐉[n2+nΛ+1:n2+2nΛ, n2+nΛ+1:n2+2nΛ] .=  BlockDiagonal(mat.(Λ_split))
            if timestep == 289
                @show y_split[1], Λ_split[1]
                @show v́⁺, 𝐰, y
                @show qr(𝐉).R |> diag
            end
        end
        # debug
        # @show norm(D*vₖ + 𝐛), norm(𝐫𝐞𝐬)
        # @show Λₖ, D*vₖ, 𝐛
        # @show Λₖ[1:3]⋅(D*vₖ + 𝐛)[1:3]

    end
    ns_stepk!
end

function solve!(sim::Simulator,solver_cache::Zhong06_CCP_Constant_Mass_Mono_Cache;
                dt,
                ftol=1e-14,xtol=ftol,
                verbose=false,verbose_contact=false,
                maxiters=50,
                max_restart=3,
                progress=true,exception=true)
    (;prob,controller,tspan,restart,totalstep) = sim
    (;bot,env) = prob
    (;structure,traj,contacts_traj) = bot
    (;M,A,contacts_bits) = solver_cache.cache
    q0 = traj.q[begin]
    λ0 = traj.λ[begin]
    q̇0 = traj.q̇[begin]
    activate_contacts!(structure,env,solver_cache,q0)
    invM = inv(M)
    pₖ₋₁ = M*q̇0
    pₖ   = zero(pₖ₋₁)
    T = eltype(q0)
    nq = length(q0)
    nλ = length(λ0)
    F = zeros(T,nq)
    ∂F∂q = zeros(T,nq,nq)
    ∂F∂q̇ = zeros(T,nq,nq)
    nx = nq + nλ
    Δx = zeros(T,nx)
    x = zero(Δx)
    Res = zero(Δx)
    Jac = zeros(T,nx,nx)
    mr = norm(M,Inf)
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
        # pₖ₋₁ = traj.p[timestep]
        # λₖ₋₁ = traj.λ[timestep]
        tₖ₋₁ = traj.t[timestep]
        qₖ   = traj.q[timestep+1]
        q̇ₖ   = traj.q̇[timestep+1]
        # pₖ   = traj.p[timestep+1]
        λₘ   = traj.λ[timestep+1]
        pₖ₋₁ = M*q̇ₖ₋₁
        qˣ = qₖ₋₁ .+ dt./2 .*q̇ₖ₋₁
        qₖ .= qₖ₋₁ .+ dt .*q̇ₖ₋₁
        q̇ₖ .= q̇ₖ₋₁
        contact_cache = activate_frictional_contacts!(structure,env,solver_cache,qˣ;checkpersist=true)
        (;na) = contact_cache.cache
        (;L) = contact_cache.cache
        isconverged = false
        normRes = typemax(T)
        iteration_break = 0
        isconverged = false
        nΛ = 3na
        n1 = nq
        n2 = n1 + nλ
        nx = n2 + 2nΛ
        Δx = zeros(T,nx)
        x = zero(Δx)
        Res = zero(Δx)
        Jac = zeros(T,nx,nx)
        Λₖ = @view x[(n2+1):n2+nΛ]
        y  = @view x[n2+nΛ+1:n2+2nΛ]
        𝐰 = zeros(T,nΛ)
        ∂y∂x = zeros(T,nΛ,n2)
        J = Diagonal(SVector(one(T),-one(T),-one(T)))
        𝐞_split = [SVector(one(T),zero(T),zero(T)) for i = 1:na]
        Λ_split = split_by_lengths(Λₖ,3)
        y_split = split_by_lengths(y,3)
        Λp = zero(Λₖ)
        yp = zero(y)
        Δxp = zeros(T,nx)
        ΔΛp = @view Δxp[(n2+1):n2+nΛ]
        Δyp = @view Δxp[n2+nΛ+1:n2+2nΛ]
        ΔΛp_split = split_by_lengths(ΔΛp,3)
        Δyp_split = split_by_lengths(Δyp,3)
        Δxc = zeros(T,nx)
        ΔΛc = @view Δxc[(n2+1):n2+nΛ]
        Δyc = @view Δxc[n2+nΛ+1:n2+2nΛ]
        ΔΛc_split = split_by_lengths(ΔΛc,3)
        Δyc_split = split_by_lengths(Δyc,3)
        get_frictional_directions_and_positions!(structure, contact_cache, qₖ₋₁, q̇ₖ₋₁, Λₖ)
        ns_stepk! = make_step_k(
            solver_cache,
            nq,nλ,na,
            qₖ₋₁,q̇ₖ₋₁,pₖ₋₁,tₖ₋₁,
            pₖ,q̇ₖ,
            invM,
            dt,mass_norm
        )
        restart_count = 0
        Λ_guess = 1.0
        while restart_count < max_restart
            Λₖ .= repeat([Λ_guess,0,0],na)
            y .= Λₖ
            x[      1:nq]          .= qₖ
            x[   nq+1:nq+nλ]       .= 0.0
            Nmax = 50
            for iteration = 1:maxiters
                ns_stepk!(
                    Res,Jac,
                    F,∂F∂q,∂F∂q̇,
                    𝐰,x,Λₖ,y,∂y∂x,
                    Λ_split,y_split,
                    structure,
                    contact_cache,
                    timestep,iteration
                )
                normRes = norm(Res,Inf)
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
                lu𝐉 = lu(Jac)
                if na == 0
                    Δx .= lu𝐉\(-Res)
                    x .+= Δx
                else # na!=0
                    get_distribution_law!(structure,contact_cache,x[1:nq])
                    if iteration < 2
                        Nmax = 50
                    else
                        Nmax = 50
                    end
                    μ = transpose(y)*Λₖ/nΛ
                    Δxp .= lu𝐉\(-Res)
                    αp_Λ = find_cone_step_length(Λ_split,ΔΛp_split,J)
                    αp_y = find_cone_step_length(y_split,Δyp_split,J)
                    αpmax = min(αp_Λ,αp_y)
                    # αpmax = find_cone_step_length(z_split,W_blocks,Δyp_split,ΔΛp_split,J)
                    αp = min(one(αpmax),0.99αpmax)
                    Λp .= Λₖ .+ αp.*ΔΛp
                    yp .= y .+ αp.*Δyp
                    μp = transpose(yp)*Λp/nΛ
                    σ = (μp/μ)^3
                    if σ == NaN || μ == 0
                        break
                    end
                    τ = σ*μp
                    𝐫𝐞𝐬_c_split = -τ.*𝐞_split.+((Δyp_split)⊙(ΔΛp_split))
                    Res[n2+nΛ+1:n2+2nΛ] .+= reduce(vcat,𝐫𝐞𝐬_c_split)
                    Δxc .= lu𝐉\(-Res)
                    # η = exp(-0.1μ) + 0.9
                    α_Λ = find_cone_step_length(Λ_split,ΔΛc_split,J)
                    # @show Λ_split,ΔΛc_split
                    α_y = find_cone_step_length(y_split,Δyc_split,J)
                    αmax = min(α_Λ,α_y)
                    α = min(1,0.99αmax)
                    # α_record[iteration] = α
                    x .+= α.*Δxc
                    μ = transpose(y)*Λₖ/nΛ
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
        pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,M,A,mass_norm,dt)
        q̇ₖ .= invM*pₖ
        Dₘ = contact_cache.cache.Dper
        Dₖ = contact_cache.cache.Dimp
        if na != 0
            update_contacts!(cₖ[contacts_bits],cₖ₋₁[contacts_bits],Dₘ*(qₖ.-qₖ₋₁).+Dₖ*q̇ₖ,2*Λₖ./(mass_norm*dt))
        end

        if !isconverged
            @warn "Newton max iterations $maxiters, at timestep=$timestep, normRes=$(normRes), restart_count=$(restart_count), num_active_contacts=$(na)"
            if exception
                @error "Not converged!"
                break
            else
                # sim.convergence = false
                # break
            end
        end

        #---------Time Step k finisher-----------
        if verbose || (na > 0 && verbose_contact)
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
