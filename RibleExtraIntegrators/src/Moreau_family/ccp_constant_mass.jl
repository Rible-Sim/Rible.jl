
struct Moreau_CCP_Constant_Mass_Cache{CacheType}
    cache::CacheType
end

function generate_cache(
        simulator::Simulator{<:DynamicsProblem{
            RobotType,
            policyType,
            EnvType,
            RestitutionFrictionCombined{NewtonRestitution,CoulombFriction}
        }},
        solver::DynamicsSolver{
            <:Moreau,
            <:AbstractBodySolver,
            <:AbstractApparatusSolver,
            <:InnerLayerContactSolver
        },
        ::Val{true};
        dt,kargs...
    )   where {RobotType,policyType,EnvType}
    (;prob) = simulator
    (;bot,policy,env) = prob
    (;structure) = bot
    (;field) = env
    options = merge(
        (checkpersist=true,), #default
        prob.options,
        solver.options,
    )
    F!(F,q,q̇,t) = gen_force!(F,bot,field,policy,q,q̇,t;)
    Jac_F!(∂F∂q̌, ∂F∂q̌̇, q, q̇, t) = gen_force_state_jacobian!(∂F∂q̌, ∂F∂q̌̇, bot, field, policy, q, q̇, t)
    
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
    ∂Aq̇∂q(q::AbstractVector,q̇) = cstr_velocity_jacobian(structure,q,q̇)
    ∂Bᵀμ∂q(q,μ) = zeros(T,nq,nq)
    (;
        contacts_bits,
        persistent_bits,
        μs_sys,
        es_sys,
        gaps_sys
    ) = prepare_contacts(bot,env)
    
    cache = @eponymtuple(
        solver,
        F!,Jac_F!,
        M,Φ,A,∂Aᵀλ∂q,∂Aq̇∂q,
        Ψ,B,∂Ψ∂q,∂Bᵀμ∂q,
        contacts_bits,
        persistent_bits,
        μs_sys,
        es_sys,
        gaps_sys
    )
    Moreau_CCP_Constant_Mass_Cache(cache)
end

function make_step_k(
        solver_cache::Moreau_CCP_Constant_Mass_Cache,
        nq,nλ,na,
        qₖ,vₖ,tₖ₊θ,
        vₖ₊₁,
        invM,
        h,mass_norm)
    (;F!,Jac_F!,M,A,Φ,∂Aᵀλ∂q,∂Aq̇∂q,solver) = solver_cache.cache
    (;θ) = solver.integrator

    n1 = nq
    n2 = nq+nλ
    nΛ = 3na
    nx = n2
    function ns_stepk!(
            𝐫𝐞𝐬,𝐉,
            F,∂F∂q,∂F∂q̇,
            𝐁,𝐛,𝐜ᵀ,𝐍,𝐫,
            x,Λₖ₊₁,
            structure,
            contact_cache,
            timestep,iteration
        )
        # @show timestep, iteration, na
        qₖ₊₁ = @view x[   1:n1]
        λₖ₊₁ = @view x[n1+1:n2]

        vₖ₊θ = (qₖ₊₁.-qₖ)./h
        vₖ₊₁ .= (1/θ)*vₖ₊θ .- (1/θ-1)*vₖ
        qₖ₊θ = (1-θ)*qₖ.+θ*qₖ₊₁

        F!(F,qₖ₊θ,vₖ₊θ,tₖ₊θ)
        Jac_F!(∂F∂q,∂F∂q̇,qₖ₊θ,vₖ₊θ,tₖ₊θ)

        Aₖ₊₁ = A(qₖ₊₁)
        

        𝐫𝐞𝐬[   1:n1] .= h.*M*(vₖ₊₁.-vₖ) .-
                        mass_norm.*transpose(Aₖ₊₁)*λₖ₊₁ .-
                        (h^2) .*F

        𝐫𝐞𝐬[n1+1:n2] .= -mass_norm.*h.*Aₖ₊₁*vₖ₊₁
        
        𝐉 .= 0.0
        𝐉[   1:n1,   1:n1] .=  1/θ .*M .-h^2 .*(θ .*∂F∂q .+ 1/h.*∂F∂q̇) .- mass_norm.*∂Aᵀλ∂q(qₖ₊₁,λₖ₊₁)
        𝐉[   1:n1,n1+1:n2] .= -mass_norm.*transpose(Aₖ₊₁)

        𝐉[n1+1:n2,   1:n1] .= -mass_norm.*(h.*∂Aq̇∂q(qₖ₊₁,vₖ₊₁) .+ 1/θ.*Aₖ₊₁)

        lu𝐉 = lu(𝐉)

        if na != 0
            (;
                H,
                restitution_coefficients,
            ) = contact_cache.cache
            Dₖ₊₁ = contact_cache.cache.Dimp
            𝐫𝐞𝐬[   1:n1] .-= h.*mass_norm.*transpose(Dₖ₊₁)*H*Λₖ₊₁ 

            𝐁 .= 0
            𝐁[   1:n1,1:nΛ] .= h.*mass_norm.*transpose(Dₖ₊₁)*H

            ∂vₖ₊₁∂q₊₁ = 1/(θ*h)*I
            v́⁺ = Dₖ₊₁*vₖ₊₁
            ∂v́⁺∂qₖ₊₁ = Dₖ₊₁*∂vₖ₊₁∂q₊₁
            𝐜ᵀ .= 0
            v́ₖ = Dₖ₊₁*vₖ
            for i = 1:na
                is = 3(i-1)
                vⁱₖ = @view v́ₖ[is+1:is+3]
                vⁱ⁺   = @view v́⁺[is+1:is+3]
                vₜⁱ⁺   = norm(vⁱ⁺[2:3])
                v́ₖₙⁱ = vⁱₖ[1]
                v́ₜⁱ = vₜⁱ⁺ + restitution_coefficients[i]*min(v́ₖₙⁱ,0)
                𝐛[is+1:is+3] .= [v́ₜⁱ,0,0]
                𝐜ᵀ[is+1     ,   1:n1] .= 1/(norm(v́⁺[is+2:is+3])+1e-14)*(v́⁺[is+2]*∂v́⁺∂qₖ₊₁[is+2,:] .+ v́⁺[is+3]*∂v́⁺∂qₖ₊₁[is+3,:])
                𝐜ᵀ[is+1:is+3,   1:n1] .+= ∂v́⁺∂qₖ₊₁[is+1:is+3,:]
                𝐜ᵀ[is+1:is+3,n1+1:n2] .= 0.0
            end
            # 𝐜ᵀinv𝐉 = 𝐜ᵀ*inv(𝐉)
            𝐍 .= 𝐜ᵀ*(lu𝐉\𝐁)
            𝐫 .= (v́⁺ + 𝐛) .-𝐜ᵀ*(lu𝐉\(𝐫𝐞𝐬 + 𝐁*Λₖ₊₁))
        end
        lu𝐉
        # debug
        # @show norm(D*vₖ + 𝐛), norm(𝐫𝐞𝐬)
        # @show Λₖ, D*vₖ, 𝐛
        # @show Λₖ[1:3]⋅(D*vₖ + 𝐛)[1:3]

    end
    ns_stepk!
end

function solve!(sim::Simulator,solver_cache::Moreau_CCP_Constant_Mass_Cache;
                dt,
                ftol=1e-14,xtol=ftol,
                verbose=false,verbose_contact=false,
                maxiters=50,max_restart=1,
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
    pₖ = M*q̇0
    pₖ   = zero(pₖ)
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
        
        #---------Time Step k Control-----------
        cₖ = contacts_traj[timestep]
        cₖ₊₁ = contacts_traj[timestep+1]
        qₖ = traj.q[timestep]
        q̇ₖ = traj.q̇[timestep]
        # pₖ = traj.p[timestep]
        # λₖ = traj.λ[timestep]
        tₖ = traj.t[timestep]
        qₖ₊₁ = traj.q[timestep+1]
        q̇ₖ₊₁ = traj.q̇[timestep+1]
        # pₖ₊₁ = traj.p[timestep+1]
        λₖ₊₁ = traj.λ[timestep+1]
        qˣ = qₖ .+ dt./2 .*q̇ₖ
        qₖ₊₁ .= qₖ .+ dt .*q̇ₖ₊₁
        q̇ₖ₊₁ .= q̇ₖ
        contact_cache = activate_frictional_contacts!(structure,env,solver_cache,qˣ;checkpersist=false)
        (;na,L) = contact_cache.cache
        isconverged = false
        normRes = typemax(T)
        iteration_break = 0
        isconverged = false
        nΛ = 3na
        Λₖ₊₁ = zeros(T,nΛ)
        Λʳₖ₊₁ = copy(Λₖ₊₁)
        ΔΛₖ₊₁ = copy(Λₖ₊₁)
        𝐁 = zeros(T,nx,nΛ)
        𝐛 = zeros(T,nΛ)
        𝐜ᵀ = zeros(T,nΛ,nx)
        𝐍 = zeros(T,nΛ,nΛ)
        𝐫 = zeros(T,nΛ)
        get_frictional_directions_and_positions!(structure, contact_cache, qₖ₊₁, q̇ₖ₊₁, Λₖ₊₁)
        ns_stepk! = make_step_k(
            solver_cache,
            nq,nλ,na,
            qₖ,q̇ₖ,tₖ,
            q̇ₖ₊₁,
            invM,
            dt,mass_norm
        )
        restart_count = 0
        Λ_guess = 1.0
        while restart_count < max_restart
            Λₖ₊₁ .= repeat([Λ_guess,0,0],na)
            x[      1:nq]          .= qₖ₊₁
            x[   nq+1:nq+nλ]       .= 0.0
            Λʳₖ₊₁ .= Λₖ₊₁
            Nmax = 50
            for iteration = 1:maxiters
                luJac = ns_stepk!(
                    Res,Jac,
                    F,∂F∂q,∂F∂q̇,
                    𝐁,𝐛,𝐜ᵀ,𝐍,𝐫,
                    x,Λₖ₊₁,
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
                    Λₖ₊₁ini = deepcopy(Λₖ₊₁)
                    Λₖ₊₁ini[begin+1:3:end] .= 0.0
                    Λₖ₊₁ini[begin+2:3:end] .= 0.0
                    if false 
                        # @show timestep, iteration
                        # @show norm(𝐍),norm(L)
                        @show norm(L*Λₖ)
                        # @show qr(L).R |> diag
                        # @show :befor, size(𝐍), rank(𝐍), cond(𝐍)
                    end
                    𝐍 .+= L
                    yₖ₊₁ini = 𝐍*Λₖ₊₁ + 𝐫
                    if false 
                        # @show :after, size(𝐍), rank(𝐍), cond(𝐍)
                        # @show yₖini
                    end
                    yₖ₊₁ini .= abs.(yₖ₊₁ini)
                    yₖ₊₁ini[begin+1:3:end] .= 0.0
                    yₖ₊₁ini[begin+2:3:end] .= 0.0
                    IPM!(Λₖ₊₁,na,nΛ,Λₖ₊₁ini,yₖ₊₁ini,𝐍,𝐫;ftol,Nmax)                    
                    ΔΛₖ₊₁ .= Λₖ₊₁ - Λʳₖ₊₁
                    minusResΛ = -Res + 𝐁*(ΔΛₖ₊₁)
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
                    Δx .= luJac\minusResΛ
                    Λʳₖ₊₁ .= Λₖ₊₁
                    x .+= Δx
                end
            end
            if isconverged
                break
            end
            restart_count += 1
            Λ_guess =  max(Λ_guess/10,maximum(abs.(Λₖ₊₁[begin:3:end])))
        end
        qₖ₊₁ .= x[      1:nq]
        λₖ₊₁ .= x[   nq+1:nq+nλ]
        Dₖ₊₁ = contact_cache.cache.Dimp
        if na != 0
            update_contacts!(cₖ₊₁[contacts_bits],cₖ[contacts_bits],Dₖ₊₁*q̇ₖ₊₁,2*Λₖ₊₁./(mass_norm*dt))
        end

        if !isconverged
            @warn "Newton max iterations $maxiters, at timestep=$timestep, normRes=$(normRes), restart_count=$(restart_count)"
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
