struct ZhongQCCP <: AbstractSolver end

struct ZhongQCCPCache{CacheType}
    cache::CacheType
end

function generate_cache(::ZhongQCCP,intor;dt,kargs...)
    (;tg) = intor.prob.bot
    M = build_M(tg) 
    M⁻¹ = build_M⁻¹(tg) 
    ∂Mq̇∂q = build_∂Mq̇∂q(tg)
    ∂M⁻¹p∂q = build_∂M⁻¹p∂q(tg)
    M! = make_M!(tg)
    M⁻¹! = make_M⁻¹!(tg)
    Jac_M! = make_Jac_M!(tg)
    Jac_M⁻¹! = make_Jac_M⁻¹!(tg)
    Φ = make_Φ(tg)
    A = make_A(tg)

    nq = size(M,2)
    T = get_numbertype(tg)
    ∂F∂q = zeros(T,nq,nq)
    ∂F∂q̇ = zeros(T,nq,nq)
    Ψ(q,q̇) = Vector{T}()
    ∂Ψ∂q(q,q̇) = Matrix{T}(undef,0,nq)
    B(q) = Matrix{T}(undef,0,nq)

    # ∂𝐌𝐚∂𝐪(q,a) = zeros(T,nq,nq)
    ∂Aᵀλ∂q(q,λ) = ∂Aᵀλ∂q̌(tg,λ)
    # ∂𝚽𝐪𝐯∂𝒒(q,v) = TR.∂Aq̇∂q(tg,v)
    ∂Bᵀμ∂q(q,μ) = zeros(T,nq,nq)
    cache = @eponymtuple(
        M,M⁻¹,∂Mq̇∂q,∂M⁻¹p∂q,
        M!,Jac_M!,M⁻¹!,Jac_M⁻¹!,
        Φ,A,Ψ,B,∂Ψ∂q,∂Aᵀλ∂q,∂Bᵀμ∂q,∂F∂q,∂F∂q̇)
    ZhongQCCPCache(cache)
end

function Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,Λₘ,Dₖ₋₁,Dₖ,H,scaling,h)
    pₖ = -pₖ₋₁ .+ 
        2/h.*Mₘ*(qₖ.-qₖ₋₁) .+ 
        scaling/h.*(transpose(A(qₖ))-transpose(A(qₖ₋₁)))*λₘ .+
        scaling.*(transpose(Dₖ)-transpose(Dₖ₋₁))*H*Λₘ
end

function make_zhongccp_ns_stepk(
        nq,nλ,na,qₖ₋₁,vₖ₋₁,pₖ₋₁,tₖ₋₁,pₖ,vₖ,
        F!,Jac_F!,get_directions_and_positions!,
        cache,h,scaling,persistent_indices,mem2act_idx
    )
    (;M!,Jac_M!,M⁻¹!,Jac_M⁻¹!,Φ,A,∂Aᵀλ∂q) = cache
    T = eltype(qₖ₋₁)
    Fₘ = zeros(T,nq)
    ∂Fₘ∂qₘ = cache.∂F∂q
    ∂Fₘ∂q̇ₘ = cache.∂F∂q̇
    ∂Mₘqₖ∂qₘ = cache.∂Mq̇∂q
    ∂Mₘq̇ₘ∂qₘ = zero(∂Mₘqₖ∂qₘ)
    ∂M⁻¹ₖpₖ∂qₖ = cache.∂M⁻¹p∂q
    Mₘ = cache.M
    M⁻¹ₘ = cache.M⁻¹
    M⁻¹ₖ = deepcopy(M⁻¹ₘ)
    n1 = nq
    n2 = nq+nλ
    nΛ = 3na
    nx = n2
    function ns_stepk!(𝐫𝐞𝐬,𝐉,𝐁,𝐛,𝐜ᵀ,𝐍,𝐫,
            x,Λₘ,Dₖ₋₁,ŕₖ₋₁,
            Dₖ,Dper, Dimp, ∂Dₖvₖ∂qₖ, ∂DᵀₖHΛₘ∂qₖ, ŕₖ,H,
            es,timestep,iteration)
        # @show timestep, iteration, na, persistent_indices
        qₖ = @view x[   1:n1]
        λₘ = @view x[n1+1:n2]
        qₘ = (qₖ.+qₖ₋₁)./2
        q̇ₘ = (qₖ.-qₖ₋₁)./h
        vₘ = q̇ₘ
        tₘ = tₖ₋₁+h/2
        M!(Mₘ,qₘ)
        Jac_M!(∂Mₘqₖ∂qₘ,qₘ,qₖ)
        F!(Fₘ,qₘ,q̇ₘ,tₘ)
        Jac_F!(∂Fₘ∂qₘ,∂Fₘ∂q̇ₘ,qₘ,q̇ₘ,tₘ)
        Aₖ₋₁ = A(qₖ₋₁)
        Aₖ   = A(qₖ)

        𝐫𝐞𝐬[   1:n1] .= h.*Mₘ*vₘ .- 
                        h.*pₖ₋₁ .-
                        (h^2)/2 .*Fₘ .-
                        scaling.*transpose(Aₖ₋₁)*λₘ .-
                        scaling*h .*transpose(Dₖ₋₁)*H*Λₘ 
        𝐫𝐞𝐬[n1+1:n2] .= scaling.*Φ(qₖ)
        
        𝐉 .= 0.0
        𝐉[   1:n1,   1:n1] .=  Mₘ .+ 1/2 .*∂Mₘqₖ∂qₘ .-h^2/2 .*(1/2 .*∂Fₘ∂qₘ .+ 1/h.*∂Fₘ∂q̇ₘ)
        𝐉[   1:n1,n1+1:n2] .= -scaling.*transpose(Aₖ₋₁)
        𝐉[n1+1:n2,   1:n1] .=  scaling.*Aₖ
        
        if na != 0
            pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,Λₘ,Dₖ₋₁,Dₖ,H,scaling,h)
            M⁻¹!(M⁻¹ₖ,qₖ) 
            vₖ .= M⁻¹ₖ*pₖ
            M⁻¹!(M⁻¹ₘ,qₘ)
            Jac_M!(∂Mₘq̇ₘ∂qₘ,qₘ,q̇ₘ)
            Jac_M⁻¹!(∂M⁻¹ₖpₖ∂qₖ,qₖ,pₖ)
            ∂Aᵀₖλₘ∂qₖ = ∂Aᵀλ∂q(qₖ,λₘ)
            get_directions_and_positions!(Dₖ,Dper, Dimp, ∂Dₖvₖ∂qₖ, ∂DᵀₖHΛₘ∂qₖ,ŕₖ,qₖ, vₖ, H*Λₘ,mem2act_idx)
            ∂pₖ∂qₖ = 2/h.*Mₘ + 
                    ∂Mₘq̇ₘ∂qₘ .+
                    scaling/(h).*∂Aᵀₖλₘ∂qₖ .+ 
                    scaling.*∂DᵀₖHΛₘ∂qₖ
            ∂vₖ∂qₖ = M⁻¹ₖ*∂pₖ∂qₖ .+ ∂M⁻¹ₖpₖ∂qₖ
            ∂vₖ∂λₘ = scaling/h.*M⁻¹ₘ*transpose(Aₖ-Aₖ₋₁)
            𝐁 .= 0
            𝐁[  1:n1,   1:nΛ] .= scaling.*h .*transpose(Dₖ₋₁)*H
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
                if i in persistent_indices
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
                v́ₜⁱ = vₜⁱ⁺ + es[i]*min(vₙⁱₖ₋₁,zero(vₙⁱₖ₋₁))
                𝐛[is+1:is+3] .= [v́ₜⁱ,0,0]
                
                Dⁱₖ = @view Dₖ[is+1:is+3,:]                
                if i in persistent_indices
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

function solve!(intor::Integrator,solvercache::ZhongQCCPCache;
        dt,
        ftol=1e-14,xtol=ftol,maxiters=50,
        verbose=false, verbose_contact=false,
        progress=true,
        exception=true,
    )
    (;prob,totalstep) = intor
    (;bot,dynfuncs) = prob
    (;traj,contacts_traj) = bot
    (;F!, Jac_F!,
        prepare_contacts!,
        get_directions_and_positions!,
        get_distribution_law!
    ) = dynfuncs
    (;cache) = solvercache
    (;M,M⁻¹,M!,M⁻¹!,A) = cache
    q0 = traj.q[begin]
    λ0 = traj.λ[begin]
    q̇0 = traj.q̇[begin]
    M!(M,q0)
    pₖ₋₁ = M*q̇0
    pₖ   = deepcopy(pₖ₋₁)
    qₖ₋½ = deepcopy(q0)
    T = eltype(q0)
    nq = length(q0)
    nλ = length(λ0)
    prepare_contacts!(q0)
    nx = nq + nλ
    Δx = zeros(T,nx)
    x = zero(Δx)
    Res = zero(Δx)
    Jac = zeros(T,nx,nx)
    mr = norm(M,Inf)
    scaling = mr
    @show mr
    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Time Step k Control-----------
        # control!(intor,cache)
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
        na,mem2act_idx,persistent_indices,contacts_bits,
        H,es,Dₖ₋₁, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕₖ₋₁, 
        L = prepare_contacts!(qₖ₋½)
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
        get_directions_and_positions!(Dₖ₋₁, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕₖ₋₁, qₖ₋₁, q̇ₖ₋₁, Λₘ, mem2act_idx,)
        Dₖ = deepcopy(Dₖ₋₁)
        ŕₖ = deepcopy(ŕₖ₋₁)
        ns_stepk! = make_zhongccp_ns_stepk(
            nq,nλ,na,qₖ₋₁,q̇ₖ₋₁,pₖ₋₁,tₖ₋₁,pₖ,q̇ₖ,
            F!,Jac_F!,
            get_directions_and_positions!,
            cache,dt,scaling,persistent_indices,mem2act_idx
        )
        restart_count = 0
        Λ_guess = 0.1
        while restart_count < 10
            Λₘ .= repeat([Λ_guess,0,0],na)
            x[      1:nq]          .= qₖ
            x[   nq+1:nq+nλ]       .= 0.0
            Λʳₖ .= Λₘ
            Nmax = 50
            for iteration = 1:maxiters
                # @show iteration,D,ηs,es,gaps
                ns_stepk!(Res,Jac,
                    𝐁,𝐛,𝐜ᵀ,𝐍,𝐫,x,Λₘ,
                    Dₖ₋₁,ŕₖ₋₁,Dₖ,Dper,Dimp,∂Dq̇∂q,∂DᵀΛ∂q,ŕₖ,H,
                    es,timestep,iteration
                )
                if na == 0
                    normRes = norm(Res)
                    if normRes < ftol
                        isconverged = true
                        iteration_break = iteration-1
                        break
                    end
                    Δx .= -Jac\Res
                    x .+= Δx
                else # na!=0
                    # @show timestep,iteration,normRes,Λₘ
                    # Λₘini = repeat([Λ_guess,0,0],na)
                    Λₘini = deepcopy(Λₘ)
                    Λₘini[begin+1:3:end] .= 0.0
                    Λₘini[begin+2:3:end] .= 0.0
                    # yini = deepcopy(Λₘini)
                    yini = 𝐍*Λₘ + 𝐫
                    yini .= abs.(yini)
                    yini[begin+1:3:end] .= 0.0
                    yini[begin+2:3:end] .= 0.0
                    IPM!(Λₘ,na,nΛ,Λₘini,yini,𝐍,𝐫;ftol=1e-14,Nmax)
                    ΔΛₖ .= Λₘ - Λʳₖ
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
                    Δx .= Jac\minusResΛ
                    Λʳₖ .= Λₘ
                    x .+= Δx
                    # @show timestep, iteration, normRes, norm(Δx), norm(ΔΛₖ),persistent_indices
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
        M!(M,qₖ₋½)
        get_directions_and_positions!(Dₖ, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕₖ, qₖ, q̇ₖ, Λₘ, mem2act_idx)
        pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,M,A,Λₘ,Dₖ₋₁,Dₖ,H,scaling,dt)
        M⁻¹!(M⁻¹,qₖ)
        q̇ₖ .= M⁻¹*pₖ
        if na != 0
            update_contacts!(cₖ[contacts_bits],cₖ₋₁[contacts_bits],Dₖ*q̇ₖ,Λₘ./(scaling*dt))
        end
        if !isconverged
            @warn "Newton max iterations $maxiters, at timestep=$timestep, normRes=$(normRes)"
            if exception
                @error "Not converged!"
                break
            else
                # intor.convergence = false
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
