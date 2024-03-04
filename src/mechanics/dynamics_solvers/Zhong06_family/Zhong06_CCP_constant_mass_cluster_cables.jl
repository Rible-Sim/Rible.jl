struct Zhong06_CCP_Constant_Mass_Cluster_Cables_Cache{solT,cacheType}
    solver::solT
    cache::cacheType
end

function generate_cache(
        simulator::Simulator{<:DynamicsProblem{
            RobotType,
            policyType,
            EnvType,
            RestitutionFrictionCombined{NewtonRestitution,CoulombFriction},
            EulerEytelwein,
        }},
        solver::DynamicsSolver{
            Zhong06,
            <:InnerLayerContactSolver,
            <:MonolithicApparatusSolver
        },
        ::Val{true};
        dt,kargs...
    ) where {RobotType,policyType,EnvType}
    (;prob) = simulator
    (;bot,policy) = prob
    (;traj,structure) = bot
    options = merge(
        (gravity=true,factor=1,checkpersist=true), #default
        prob.options,
        solver.options,
    )
    @show "YesYesYes"
    (;M,M⁻¹)= build_mass_matrices(structure)
    A = make_cstr_jacobian(structure)
    Φ = make_cstr_function(structure)
    ψ = build_ψ(structure)
    F!(F,q,s,q̇,t) = generalized_force!(F,bot,policy,q,q̇,t,s;gravity=options.gravity)
    Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t) = generalized_force_jacobian!(∂F∂q,∂F∂q̇,bot,policy,q,q̇,t)
    ∂Aᵀλ∂q(q::AbstractVector,λ) = cstr_forces_jacobian(structure,q,λ)
    q0 = traj.q[begin]
    λ0 = traj.λ[begin]
    q̇0 = traj.q̇[begin]
    q̇0 = traj.q̇[begin]
    p0 = traj.p[begin] .= M*q̇0
    s0 = traj.s[begin]
    T = eltype(q0)
    nq = length(q0)
    nλ = length(λ0)
    ns = length(s0)
    ∂F∂q = zeros(T,nq,nq)
    ∂F∂q̇ = zeros(T,nq,nq)
    nx = nq + nλ + ns
    xₖ = vcat(q0,λ0,s0)
    Res = zero(xₖ)
    Jac = Res*transpose(Res)
    Zhong06_CCP_Constant_Mass_Cluster_Cables_Cache(solver,
        @eponymtuple(
            F!,Jac_F!,∂Aᵀλ∂q,
            M,M⁻¹,A,Φ,ψ,
            ∂F∂q,∂F∂q̇,
            nq,nλ,ns,nx,
            xₖ,
            Res,
            Jac,
            options
        )
    )
end

function make_step_k(
        solver_cache::Zhong06_CCP_Constant_Mass_Cluster_Cables_Cache,
        nq,nλ,ns,na,
        qₖ₋₁,vₖ₋₁,pₖ₋₁,tₖ₋₁,
        pₖ,vₖ,
        M⁻¹,
        h,mass_norm)
    (;F!,Jac_F!,M,Φ,A,ψ,∂Aᵀλ∂q) = solver_cache.cache

    n1 = nq
    n2 = nq+nλ
    n3 = nq+nλ+ns
    nΛ = 3na
    nx = n3

    function ns_stepk!(
            𝐫𝐞𝐬,𝐉,
            F,∂F∂q,∂F∂q̇,
            𝐁,𝐛,𝐜ᵀ,𝐲,
            x,Λₖ,
            structure,
            contact_cache,
            timestep,iteration,doin=true
        )
        # @show timestep, iteration, na
        qₖ = @view x[   1:nq   ]
        λₖ = @view x[nq+1:nq+nλ]
        sₖ = @view x[nq+nλ+1:nx]
        qₘ = (qₖ.+qₖ₋₁)./2
        q̇ₘ = (qₖ.-qₖ₋₁)./h
        tₘ = tₖ₋₁+h/2

        vₘ = q̇ₘ
    
        F!(F,qₘ, sₖ, q̇ₘ, tₘ)

        Jac_F!(∂F∂q,∂F∂q̇,qₘ,q̇ₘ,tₘ)

        Aₖ₋₁ = A(qₖ₋₁)
        Aᵀₖ₋₁ = transpose(Aₖ₋₁)
        Aₖ   = A(qₖ)

        𝐫𝐞𝐬[   1:nq   ] .= M*(qₖ.-qₖ₋₁) .-
                           h.*pₖ₋₁ .-
                           (h^2)/2 .*F .+
                           mass_norm.*Aᵀₖ₋₁*λₖ
        𝐫𝐞𝐬[nq+1:nq+nλ] .= mass_norm.*Φ(qₖ)
        𝐫𝐞𝐬[nq+nλ+1:nx] .= ψ(sₖ)


        Jac_F!(∂F∂q,∂F∂q̇,qₘ,q̇ₘ,tₘ)
        ζ = build_ζ(structure)
        ∂ζ∂q = build_∂ζ∂q(structure)
        ∂ζ∂s̄ = build_∂ζ∂s̄(structure)
        n = length(ζ)
        ∂s̄∂s̄ = I(n)
        ∂F∂s̄ = build_∂Q̌∂s̄(structure)
        κ₁ = 10; κ₂ = 10
        coes = diagm([((ζ[i]/κ₁)^2 + (κ₂*sₖ[i])^2)^(-1/2) for i in 1:n])
        coζ = coes*diagm([ζ[i]/κ₁^2 for i in 1:n]) - diagm([1/κ₁ for i in 1:n])
        cos̄ = coes*diagm([κ₂^2*sₖ[i] for i in 1:n]) - diagm([κ₂ for i in 1:n])

        𝐉 .= 0.0
        𝐉[   1:nq ,   1:nq ] .=  M.-(h^2)/2 .*(1/2 .*∂F∂q.+1/h.*∂F∂q̇)
        𝐉[   1:nq ,nq+1:nq+nλ] .=  mass_norm.*Aᵀₖ₋₁
        𝐉[   1:nq ,nq+nλ+1:end] .=  -1/2 * h^2 .* ∂F∂s̄

        𝐉[nq+1:nq+nλ,   1:nq ] .=  mass_norm.*Aₖ
        𝐉[nq+1:nq+nλ,   nq+1:nq+nλ ] .= 0.0 
        𝐉[nq+1:nq+nλ,   nq+nλ+1:end ] .= 0.0 

        𝐉[nq+nλ+1:end,1:nq] .= 1/2 * coζ*∂ζ∂q 
        𝐉[nq+nλ+1:end,nq+1:nq+nλ] .= 0.0
        𝐉[nq+nλ+1:end,nq+nλ+1:end] .= coζ*∂ζ∂s̄ + cos̄*∂s̄∂s̄

        lu𝐉 = lu(𝐉)

        if (na != 0)
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

            if doin
                pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,M,A,mass_norm,h)
                vₖ .= M⁻¹*pₖ        
                ∂vₘ∂qₖ = 1/h*I
                ∂vₖ∂qₖ = 2/h*I  + mass_norm/(h).*M⁻¹*(∂Aᵀλ∂q(qₖ,λₘ))
                ∂vₖ∂λₘ = mass_norm.*M⁻¹*transpose(Aₖ-Aₖ₋₁)/(h)
                
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
                𝐲 .= (v́⁺ + 𝐛)
            end
        end
        lu𝐉
        # debug
        # @show norm(D*vₖ + 𝐛), norm(𝐫𝐞𝐬)
        # @show Λₖ, D*vₖ, 𝐛
        # @show Λₖ[1:3]⋅(D*vₖ + 𝐛)[1:3]

    end
    ns_stepk!
end

function solve!(sim::Simulator,solver_cache::Zhong06_CCP_Constant_Mass_Cluster_Cables_Cache;
                dt,ftol=1e-14,verbose=false,maxiters=50,
                progress=true,exception=true)
    (;prob,controller,tspan,restart,totalstep) = sim
    (;bot,) = prob
    (;structure,traj) = bot
    (;
        M,M⁻¹,A,Φ,ψ,
        ∂F∂q,∂F∂q̇,
        nq,nλ,ns,nx,
        xₖ,
        Res,
        Jac
    ) = solver_cache.cache
    mr = norm(M,Inf)
    mass_norm = mr
    h = dt
    T = get_numbertype(bot)

    iteration_break = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    dg_step = ceil(Int,log10(totalstep))+1
    dg_dt = max(1,-floor(Int,log10(h)))
    wd_t = ceil(Int,log10(traj.t[end]))+dg_dt+1+1
    progfmt = Printf.Format("Progress: %5.1f%%, step: %$(dg_step)u, time: %$(wd_t).$(dg_dt)f, iterations: %s \n")
    
    for timestep = 1:totalstep
        #---------Step k Control-----------
        # control!(sim,cache)
        #---------Step k Control-----------
        tₖ₋₁ = traj.t[timestep]
        tₖ   = traj.t[timestep+1]
        qₖ₋₁ = traj.q[timestep]
        sₖ₋₁ = traj.s[timestep]
        qₖ₋₁ = traj.q[timestep]
        q̇ₖ₋₁ = traj.q̇[timestep]
        pₖ₋₁ = traj.p[timestep] 
        qₖ = traj.q[timestep+1]
        sₖ = traj.s[timestep+1]
        qₖ = traj.q[timestep+1]
        q̇ₖ = traj.q̇[timestep+1]
        λₖ = traj.λ[timestep+1]
        pₖ = traj.p[timestep+1]
        F = traj.F[timestep+1]
        xₖ[   1:nq]    .= qₖ₋₁ .+ h.*q̇ₖ₋₁
        xₖ[nq+1:nq+nλ] .= 0
        xₖ[nq+nλ+1:nx] .= sₖ₋₁
        

        na = 0 #INFO 激活的接触点数目
        nΛ = 3na
        Λₖ = zeros(T,nΛ)
        Λʳₖ = copy(Λₖ)
        ΔΛₖ = copy(Λₖ)
        𝐁 = zeros(T,nx,nΛ)
        𝐁t = zeros(T,nx,nΛ)
        𝐛 = zeros(T,nΛ)
        𝐜ᵀ = zeros(T,nΛ,nx)
        𝐍 = zeros(T,nΛ,nΛ)
        𝐲 = zeros(T,nΛ)

        ## Res_stepk! = make_Res_stepk(qₖ,qₖ,sₖ₋₁,qₖ₋₁,pₖ₋₁,F,tₖ₋₁)
        #INFO 已将Res，Jac合成一个函数
        ns_stepk! = make_step_k(
            solver_cache,
            nq,nλ,ns,na,
            qₖ₋₁,q̇ₖ₋₁,pₖ₋₁,tₖ₋₁,
            pₖ,q̇ₖ,
            M⁻¹,
            dt,mass_norm
        )
        isconverged = false
        normRes = typemax(T)
        
        for iteration = 1:maxiters
            ## Res_stepk!(Res,xₖ)
            ## Jac_stepk!(Jac,xₖ,qₖ,sₖ,qₖ,qₖ₋₁,∂F∂q,∂F∂q̇,tₖ₋₁)
            contact_cache = nothing
            luJac = ns_stepk!(
                Res,Jac,
                F,∂F∂q,∂F∂q̇,
                𝐁,𝐛,𝐜ᵀ,𝐲,
                xₖ,Λʳₖ,
                structure,
                contact_cache,
                timestep,iteration
            )

            normRes = norm(Res)
            if normRes < ftol
                isconverged = true
                iteration_break = iteration-1
                break
            end                
            # FiniteDiff.finite_difference_jacobian!(Jac,Res_stepk!,xₖ,Val{:central})
            xₖ .+= luJac\(-Res)
        end
        
        if !isconverged
            if exception
                error("Not Converged! Step=$timestep, normRes=$normRes")
            else
                # sim.convergence = false
                break
            end
        end
        qₖ .= xₖ[   1:nq   ]
        λₖ .= xₖ[nq+1:nq+nλ]
        sₖ .= xₖ[nq+nλ+1:nx]
        pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₖ,M,A,mass_norm,h)
        #---------Step k finisher-----------
        #---------Step k finisher-----------
        if verbose
            progstr = Printf.format(
                progfmt,
                floor(timestep/totalstep*100;digits=1), timestep, tₖ, iteration_break
            )
            print(progstr)
        end
        next!(prog)
    end

    bot
end
