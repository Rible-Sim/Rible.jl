
struct Zhong06_CCP_Constant_Mass_Cache{cacheType}
    cache::cacheType
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
        ::Val{true};
        dt,kargs...
    )   where {RobotType,EnvType}
    (;prob) = simulator
    (;bot,env) = prob
    (;structure) = bot
    options = merge(
        (gravity=true,factor=1,checkpersist=true), #default
        prob.options,
        solver.options,
    )
    F!(F,q,q̇,t) = generalized_force!(F,bot,q,q̇,t;gravity=options.gravity)
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
        gaps_sys,
        options,
    )
    Zhong06_CCP_Constant_Mass_Cache(cache)
end

function make_step_k(
        solver_cache::Zhong06_CCP_Constant_Mass_Cache,
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
            𝐁,𝐛,𝐜ᵀ,𝐍,𝐫,
            x,Λₖ,
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

        lu𝐉 = lu(𝐉)

        if na != 0
            (;
                H,
                restitution_coefficients,
                D,
            ) = contact_cache.cache
            Dₘ = contact_cache.cache.Dper
            Dₖ = contact_cache.cache.Dimp
            𝐁 .= 0
            𝐁[   1:n1,1:nΛ] .= h.*mass_norm.*transpose(D)*H
            𝐫𝐞𝐬  .-= 𝐁*Λₖ 

            pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,M,A,mass_norm,h)
            vₖ .= invM*pₖ        
            ∂vₘ∂qₖ = 1/h*I
            ∂vₖ∂qₖ = 2/h*I  + mass_norm/(h).*invM*(∂Aᵀλ∂q(qₖ,λₘ))
            ∂vₖ∂λₘ = mass_norm.*invM*transpose(Aₖ-Aₖ₋₁)/(h)
            
            v́⁺ = Dₘ*vₘ .+ Dₖ*vₖ
            ∂v́⁺∂qₖ = Dₘ*∂vₘ∂qₖ .+ Dₖ*∂vₖ∂qₖ
            𝐜ᵀ .= 0
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
                𝐛[is+1:is+3] .= [v́ₜⁱ,0,0]
                Dⁱₘ = @view Dₘ[is+1:is+3,:]
                Dⁱₖ = @view Dₖ[is+1:is+3,:]
                𝐜ᵀ[is+1     ,   1:n1] .= 1/(norm(v́⁺[is+2:is+3])+1e-14)*(v́⁺[is+2]*∂v́⁺∂qₖ[is+2,:] .+ v́⁺[is+3]*∂v́⁺∂qₖ[is+3,:])
                𝐜ᵀ[is+1:is+3,   1:n1] .+= ∂v́⁺∂qₖ[is+1:is+3,:]
                𝐜ᵀ[is+1:is+3,n1+1:n2] .= Dⁱₖ*∂vₖ∂λₘ
            end
            # 𝐜ᵀinv𝐉 = 𝐜ᵀ*inv(𝐉)
            𝐍 .= 𝐜ᵀ*(lu𝐉\𝐁)
            𝐫 .= (v́⁺ + 𝐛) .-𝐜ᵀ*(lu𝐉\(𝐫𝐞𝐬 + 𝐁*Λₖ))
        end
        lu𝐉
        # debug
        # @show norm(D*vₖ + 𝐛), norm(𝐫𝐞𝐬)
        # @show Λₖ, D*vₖ, 𝐛
        # @show Λₖ[1:3]⋅(D*vₖ + 𝐛)[1:3]

    end
    ns_stepk!
end

function solve!(sim::Simulator,solver_cache::Zhong06_CCP_Constant_Mass_Cache;
                dt,
                ftol=1e-14,xtol=ftol,
                verbose=false,verbose_contact=false,
                maxiters=50,
                max_restart=3,
                progress=true,exception=true)
    (;prob,controller,tspan,restart,totalstep) = sim
    (;bot,env) = prob
    (;structure,traj,contacts_traj) = bot
    (;M,A,contacts_bits,options) = solver_cache.cache
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
        contact_cache = activate_frictional_contacts!(structure,env,solver_cache,qˣ;checkpersist=options.checkpersist)
        (;na) = contact_cache.cache
        (;L) = contact_cache.cache
        isconverged = false
        normRes = typemax(T)
        iteration_break = 0
        isconverged = false
        nΛ = 3na
        Λₖ = zeros(T,nΛ)
        Λʳₖ = copy(Λₖ)
        ΔΛₖ = copy(Λₖ)
        𝐁 = zeros(T,nx,nΛ)
        𝐛 = zeros(T,nΛ)
        𝐜ᵀ = zeros(T,nΛ,nx)
        𝐍 = zeros(T,nΛ,nΛ)
        𝐫 = zeros(T,nΛ)
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
        Λ_guess = 0.1
        while restart_count < max_restart
            Λₖ .= repeat([Λ_guess,0,0],na)
            x[      1:nq]          .= qₖ
            x[   nq+1:nq+nλ]       .= 0.0
            Λʳₖ .= Λₖ
            Nmax = 50
            for iteration = 1:maxiters
                luJac = ns_stepk!(
                    Res,Jac,
                    F,∂F∂q,∂F∂q̇,
                    𝐁,𝐛,𝐜ᵀ,𝐍,𝐫,
                    x,Λₖ,
                    structure,
                    contact_cache,
                    timestep,iteration
                )
                normRes = norm(Res)
                if na == 0
                    if  normRes < ftol
                        isconverged = true
                        iteration_break = iteration-1
                        break
                    end
                    Δx .= luJac\(-Res)
                    x .+= Δx
                else # na!=0
                    get_distribution_law!(structure,contact_cache,x[1:nq])
                    if iteration < 2
                        Nmax = 50
                    else
                        Nmax = 50
                    end
                    # Λₖini = repeat([Λ_guess,0,0],na)
                    Λₖini = deepcopy(Λₖ)
                    Λₖini[begin+1:3:end] .= 0.0
                    Λₖini[begin+2:3:end] .= 0.0
                    if false
                        # @show timestep, iteration
                        # @show 𝐍
                        @show iteration
                        # @show qr(L).R |> diag
                        @show :befor, size(𝐍), rank(𝐍), cond(𝐍)
                    end
                    𝐍 .+= L
                    yₖini = 𝐍*Λₖ + 𝐫
                    yₖini .= abs.(yₖini)
                    yₖini[begin+1:3:end] .= 0.0
                    yₖini[begin+2:3:end] .= 0.0
                    IPM!(Λₖ,na,nΛ,Λₖini,yₖini,𝐍,𝐫;ftol,Nmax)
                    ΔΛₖ .= Λₖ - Λʳₖ
                    minusResΛ = -Res + 𝐁*(ΔΛₖ)
                    normRes = norm(minusResΛ)
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
                    if false
                        @show :after, size(𝐍), rank(𝐍), cond(𝐍)
                        @show minusResΛ
                        @show Λₖ, Λʳₖ
                        @show cond(Jac)
                        # @show yₖini
                        # display(L)
                        # display(nullspace(L))
                        # @show Λₖ
                        # display(L*Λₖ)
                    end
                    Δx .= luJac\minusResΛ
                    Λʳₖ .= Λₖ
                    x .+= Δx
                    # @show timestep, iteration, normRes, norm(Δx), norm(ΔΛₖ)
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
