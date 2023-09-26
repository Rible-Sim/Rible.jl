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

function Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,Λₘ,Dₖ₋₁,Dₖ,scalingΛ,h)
    pₖ = -pₖ₋₁ .+ 
        2/h.*Mₘ*(qₖ.-qₖ₋₁) .+ 
        scalingΛ/h.*(transpose(A(qₖ))-transpose(A(qₖ₋₁)))*λₘ .+
        scalingΛ/h.*(transpose(Dₖ)-transpose(Dₖ₋₁))*Λₘ
end

function make_zhongccp_ns_stepk(
        nq,nλ,na,qₖ₋₁,vₖ₋₁,pₖ₋₁,tₖ₋₁,pₖ,vₖ,
        F!,Jac_F!,get_directions_and_positions,get_∂Dq̇∂q,get_∂DᵀΛ∂q,
        cache,h,scalingΛ,persistent_indices
    )
    (;M!,Jac_M!,M⁻¹!,Jac_M⁻¹!,Φ,A,∂Aᵀλ∂q) = cache
    ∂Fₘ∂qₘ = cache.∂F∂q
    ∂Fₘ∂q̇ₘ = cache.∂F∂q̇
    ∂Mₘqₖ∂qₘ = cache.∂Mq̇∂q
    ∂M⁻¹ₖpₖ∂qₖ = cache.∂M⁻¹p∂q
    Mₘ = cache.M
    M⁻¹ₘ = cache.M⁻¹
    M⁻¹ₖ = deepcopy(M⁻¹ₘ)
    n1 = nq
    n2 = nq+nλ
    nΛ = 3na
    nx = n2
    function ns_stepk!(𝐫𝐞𝐬,𝐉,𝐁,𝐛,𝐜ᵀ,𝐍,𝐫,x,𝚲ₘ,Dₖ₋₁,ŕₖ₋₁,H,filtered_gaps,es,timestep,iteration)
        # @show timestep, iteration, na, persistent_indices
        qₖ = @view x[   1:n1]
        λₘ = @view x[n1+1:n2]
        T = eltype(qₖ)
        qₘ = (qₖ.+qₖ₋₁)./2
        q̇ₘ = (qₖ.-qₖ₋₁)./h
        vₘ = q̇ₘ
        tₘ = tₖ₋₁+h/2
        Fₘ = zeros(T,nq)
        M!(Mₘ,qₘ)
        Jac_M!(∂Mₘqₖ∂qₘ,qₘ,qₖ)
        F!(Fₘ,qₘ,q̇ₘ,tₘ)
        Jac_F!(∂Fₘ∂qₘ,∂Fₘ∂q̇ₘ,qₘ,q̇ₘ,tₘ)
        Aₖ₋₁ = A(qₖ₋₁)
        Aₖ   = A(qₖ)

        𝐫𝐞𝐬[   1:n1] .= Mₘ*(qₖ.-qₖ₋₁) .- 
                        h.*pₖ₋₁ .-
                        (h^2)/2 .*Fₘ .-
                        scalingΛ .*transpose(Aₖ₋₁)*λₘ .-
                        scalingΛ .*transpose(Dₖ₋₁)*H*𝚲ₘ 
        𝐫𝐞𝐬[n1+1:n2] .= scalingΛ .*Φ(qₖ)
        
        𝐉 .= 0.0
        𝐉[   1:n1,   1:n1] .=  Mₘ .+ 1/2 .*∂Mₘqₖ∂qₘ .-h^2/2 .*(1/2 .*∂Fₘ∂qₘ .+ 1/h.*∂Fₘ∂q̇ₘ)
        𝐉[   1:n1,n1+1:n2] .= -scalingΛ .*transpose(Aₖ₋₁)
        𝐉[n1+1:n2,   1:n1] .=  scalingΛ .*Aₖ
        
        if na != 0
            Dₖ,ŕₖ = get_directions_and_positions(qₖ)
            pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,𝚲ₘ,Dₖ₋₁,Dₖ,scalingΛ,h)
            M⁻¹!(M⁻¹ₖ,qₖ) 
            vₖ .= M⁻¹ₖ*pₖ
            M⁻¹!(M⁻¹ₘ,qₘ)    
            Jac_M⁻¹!(∂M⁻¹ₖpₖ∂qₖ,qₖ,vₖ)
            ∂Aᵀₖλₘ∂qₖ = ∂Aᵀλ∂q(qₖ,λₘ)
            ∂DᵀₖHΛₘ∂qₖ = get_∂DᵀΛ∂q(qₖ,H*𝚲ₘ)
            ∂Mₘq̇ₘ∂qₘ = zero(∂Mₘqₖ∂qₘ)
            Jac_M!(∂Mₘq̇ₘ∂qₘ,qₘ,q̇ₘ)
            ∂pₖ∂qₖ = 2/h.*Mₘ + 
                    ∂Mₘq̇ₘ∂qₘ .+
                    1/(h).*∂Aᵀₖλₘ∂qₖ .+ 
                    1/(h).*∂DᵀₖHΛₘ∂qₖ
            ∂vₖ∂qₖ = M⁻¹ₖ*∂pₖ∂qₖ .+ ∂M⁻¹ₖpₖ∂qₖ
            ∂vₖ∂λₘ = M⁻¹ₘ*transpose(Aₖ-Aₖ₋₁)/(h)
            𝐁 .= 0
            𝐁[  1:n1,   1:nΛ] .= scalingΛ .*transpose(Dₖ₋₁)*H
            ∂Dₖvₖ∂qₖ = get_∂Dq̇∂q(qₖ,vₖ)
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
                # @show timestep,iteration, vₙⁱₖ₋₁, vₙⁱ, vₜⁱₖ₋₁, vₜⁱ, 𝚲ₘ
                v́ₜⁱ = vₜⁱ⁺ + es[i]*min(vₙⁱₖ₋₁,zero(T))
                𝐛[is+1:is+3] .= [v́ₜⁱ+filtered_gaps[i],0,0]
                
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
            # @show 𝚲ₘ, D*vₖ, 𝐛
            # @show v́ₖ
            # @show 𝚲ₘ[1:3]⋅(v́ₖ + 𝐛)[1:3]
            𝐫 .= (v́⁺ + 𝐛) - 𝐜ᵀinv𝐉*(𝐫𝐞𝐬 + 𝐁*𝚲ₘ)
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
    F!, Jac_F!, prepare_contacts!,get_directions_and_positions,get_∂Dq̇∂q,get_∂DᵀΛ∂q = dynfuncs
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
    prepare_contacts!(contacts_traj[end],q0)
    nx = nq + nλ
    Δx = zeros(T,nx)
    x = zero(Δx)
    Res = zero(Δx)
    Jac = zeros(T,nx,nx)
    mr = norm(M,Inf)
    scaling = mr

    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Time Step k Control-----------
        # control!(intor,cache)
        #---------Time Step k Control-----------
        push!(contacts_traj,deepcopy(contacts_traj[end]))
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
        active_contacts,gaps,H,es = prepare_contacts!(cₖ,qₖ₋½)
        na = length(active_contacts)
        Dₖ₋₁,ŕₖ₋₁ = get_directions_and_positions(active_contacts,qₖ₋₁)
        persistent_indices = findall((c)->c.state.persistent,active_contacts)
        Dₖ = deepcopy(Dₖ₋₁)
        ŕₖ = deepcopy(ŕₖ₋₁)
        # @show ŕₖ₋₁
        # Dₘ = copy(D)
        # Dₖ = zero(D)
        filtered_gaps = zero(gaps)
        # empty!(persistent_indices)
        # persistent_indices = [1]
        if (na !== 0) && !isempty(persistent_indices)
            # @show timestep,persistent_indices
            epi = reduce(vcat,[collect(3(i-1)+1:3i) for i in persistent_indices])
            # Dₘ[epi,:] .= D[epi,:]
            # Dₖ[epi,:] .= 0
            # filtered_gaps[persistent_indices] = gaps[persistent_indices]
        end
        isconverged = false
        normRes = typemax(T)
        iteration_break = 0
        nΛ = 3na
        𝚲ₘ = zeros(T,nΛ)
        𝚲ʳₖ = copy(𝚲ₘ)
        Δ𝚲ₖ = copy(𝚲ₘ)
        𝐁 = zeros(T,nx,nΛ)
        𝐛 = zeros(T,nΛ)
        𝐜ᵀ = zeros(T,nΛ,nx)
        𝐍 = zeros(T,nΛ,nΛ)
        𝐫 = zeros(T,nΛ)
        scalingΛ = dt
        get_directions_and_positions_active(q) = get_directions_and_positions(active_contacts,q)
        get_∂Dq̇∂q_active(q,q̇) = get_∂Dq̇∂q(active_contacts,q,q̇)
        get_∂DᵀΛ∂q_active(q,Λ) = get_∂DᵀΛ∂q(active_contacts,q,Λ)
        ns_stepk! = make_zhongccp_ns_stepk(
            nq,nλ,na,qₖ₋₁,q̇ₖ₋₁,pₖ₋₁,tₖ₋₁,pₖ,q̇ₖ,
            F!,Jac_F!,get_directions_and_positions_active,get_∂Dq̇∂q_active,get_∂DᵀΛ∂q_active,
            cache,dt,scalingΛ,persistent_indices
        )
        restart_count = 0
        𝚲_guess = 10.0
        while restart_count < 10
            x[      1:nq]          .= qₖ
            x[   nq+1:nq+nλ]       .= 0.0
            𝚲ₘ .= repeat([𝚲_guess,0,0],na)
            𝚲ʳₖ .= 0.0
            Nmax = 50
            for iteration = 1:maxiters
                # @show iteration,D,ηs,es,gaps
                ns_stepk!(Res,Jac,𝐁,𝐛,𝐜ᵀ,𝐍,𝐫,x,𝚲ₘ,Dₖ₋₁,ŕₖ₋₁,H,filtered_gaps,es,timestep,iteration)
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
                    if timestep == 308 
                        @show timestep,iteration,normRes,norm(Res),𝚲ₘ
                        # # 𝚲ₘini = repeat([1.0,0,0],na)
                        # 𝚲ₘini = deepcopy(𝚲ʳₖ)
                        Nmax = 1000
                    end
                    # else
                    # end
                    𝚲ₘini = repeat([𝚲_guess,0,0],na)
                    𝚲ₘini[begin+1:3:end] .= 0.0
                    𝚲ₘini[begin+2:3:end] .= 0.0
                    # yini = deepcopy(𝚲ₘini)
                    yini = 𝐍*𝚲ₘ + 𝐫
                    yini .= abs.(yini)
                    yini[begin+1:3:end] .= 0.0
                    yini[begin+2:3:end] .= 0.0
                    IPM!(𝚲ₘ,na,nΛ,𝚲ₘini,yini,𝐍,𝐫;ftol=1e-14,Nmax)
                    Δ𝚲ₖ .= 𝚲ₘ - 𝚲ʳₖ
                    minusRes𝚲 = -Res + 𝐁*(Δ𝚲ₖ)
                    normRes = norm(minusRes𝚲)
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
                    Δx .= Jac\minusRes𝚲
                    𝚲ʳₖ .= 𝚲ₘ
                    x .+= Δx
                    # @show timestep, iteration, normRes, norm(Δx), norm(Δ𝚲ₖ),persistent_indices
                end
            end
            if isconverged
                break
            end
            restart_count += 1
            𝚲_guess /= 10
            # @warn "restarting step: $timestep, count: $restart_count, 𝚲_guess = $𝚲_guess"
        end
        qₖ .= x[      1:nq]
        λₘ .= x[   nq+1:nq+nλ]
        qₖ₋½ .= (qₖ.+qₖ₋₁)./2
        M!(M,qₖ₋½)
        pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,M,A,𝚲ₘ,Dₖ₋₁,Dₖ,scalingΛ,dt)
        M⁻¹!(M⁻¹,qₖ)
        q̇ₖ .= M⁻¹*pₖ
        if na != 0
            update_contacts!(active_contacts,(ŕₖ.-ŕₖ₋₁)./dt.+Dₖ*q̇ₖ,𝚲ₘ./scalingΛ)
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
