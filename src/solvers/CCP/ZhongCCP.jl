struct ZhongCCP <: AbstractSolver end

struct ZhongCCPCache{CacheType}
    cache::CacheType
end

function generate_cache(::ZhongCCP,intor;dt,kargs...)
    (;prob) = intor
    (;bot,dynfuncs) = prob
    (;tg) = bot
    M = Matrix(build_M(tg))
    Φ = make_Φ(bot)
    A = make_A(bot)

    nq = size(M,2)
    T = get_numbertype(bot)
    Ψ(q,q̇) = Vector{T}()
    ∂Ψ∂q(q,q̇) = Matrix{T}(undef,0,nq)
    B(q) = Matrix{T}(undef,0,nq)

    # ∂𝐌𝐚∂𝐪(q,a) = zeros(T,nq,nq)
    ∂Aᵀλ∂q(q,λ) = ∂Aᵀλ∂q̌(tg,λ)
    # ∂𝚽𝐪𝐯∂𝒒(q,v) = TR.∂Aq̇∂q(tg,v)
    ∂Bᵀμ∂q(q,μ) = zeros(T,nq,nq)
    cache = @eponymtuple(M,Φ,A,Ψ,B,∂Ψ∂q,∂Aᵀλ∂q,∂Bᵀμ∂q)
    ZhongCCPCache(cache)
end

function Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,h)
    pₖ = -pₖ₋₁ .+ 2/h.*Mₘ*(qₖ.-qₖ₋₁) .+ 1/(h).*(transpose(A(qₖ))-transpose(A(qₖ₋₁)))*λₘ
end

function make_zhongccp_ns_stepk(nq,nλ,na,qₖ₋₁,vₖ₋₁,pₖ₋₁,tₖ₋₁,pₖ,vₖ,dynfuncs,cache,invM,h,scalingΛ,persistent_indices)
    F!,Jac_F!,_ = dynfuncs
    (;M,Φ,A,∂Aᵀλ∂q) = cache

    n1 = nq
    n2 = nq+nλ
    nΛ = 3na
    nx = n2
    function ns_stepk!(𝐫𝐞𝐬,𝐉,𝐁,𝐛,𝐜ᵀ,𝐍,𝐫,x,𝚲ₖ,D,Dₘ,Dₖ,H,filtered_gaps,es,timestep,iteration)
        # @show timestep, iteration, na, persistent_indices
        qₖ = @view x[   1:n1]
        λₘ = @view x[n1+1:n2]
        qₘ = (qₖ.+qₖ₋₁)./2
        q̇ₘ = (qₖ.-qₖ₋₁)./h
        vₘ = q̇ₘ
        tₘ = tₖ₋₁+h/2
        T = eltype(qₖ)
        Fₘ = zeros(T,nq)
        F!(Fₘ,qₘ,q̇ₘ,tₘ)
        ∂F∂q = zeros(T,nq,nq)
        ∂F∂q̇ = zeros(T,nq,nq)
        Jac_F!(∂F∂q,∂F∂q̇,qₘ,q̇ₘ,tₘ)

        Aₖ₋₁ = A(qₖ₋₁)
        Aₖ   = A(qₖ)

        ∂vₘ∂qₖ = 1/h*I

        ∂vₖ∂qₖ = 2/h*I + 1/(h).*invM*(∂Aᵀλ∂q(qₖ,λₘ))
        ∂vₖ∂λₘ = invM*transpose(Aₖ-Aₖ₋₁)/(h)

        ∂DᵀHΛₖ∂qₖ = zeros(T,nq,nq)

        𝐫𝐞𝐬[   1:n1] .= -h.*pₖ₋₁ .+ M*(qₖ.-qₖ₋₁) .-
                        scalingΛ .*transpose(D)*H*𝚲ₖ .-
                        transpose(Aₖ₋₁)*λₘ .-
                        (h^2)/2 .*Fₘ

        𝐫𝐞𝐬[n1+1:n2] .= Φ(qₖ)
        
        𝐉 .= 0.0
        𝐉[   1:n1,   1:n1] .=  M .- scalingΛ .*∂DᵀHΛₖ∂qₖ .-h^2/2 .*(1/2 .*∂F∂q .+ 1/h.*∂F∂q̇)
        𝐉[   1:n1,n1+1:n2] .= -transpose(Aₖ₋₁)

        𝐉[n1+1:n2,   1:n1] .=  Aₖ

        𝐁 .= 0
        𝐁[   1:n1,1:nΛ] .= scalingΛ .*transpose(D)*H

        
        pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,M,A,h)
        vₖ .= invM*pₖ        
        if na != 0
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
                # @show timestep,iteration, vₙⁱₖ₋₁, vₙⁱ, vₜⁱₖ₋₁, vₜⁱ, 𝚲ₖ
                v́ₜⁱ = vₜⁱ⁺ + es[i]*min(vₙⁱₖ₋₁,zero(T))
                𝐛[is+1:is+3] .= [v́ₜⁱ+filtered_gaps[i],0,0]
                
                Dⁱₘ = @view Dₘ[is+1:is+3,:]
                Dⁱₖ = @view Dₖ[is+1:is+3,:]
                𝐜ᵀ[is+1     ,   1:n1] .= 1/(norm(v́⁺[is+2:is+3])+1e-14)*(v́⁺[is+2]*∂v́⁺∂qₖ[is+2,:] .+ v́⁺[is+3]*∂v́⁺∂qₖ[is+3,:])
                𝐜ᵀ[is+1:is+3,   1:n1] .+= ∂v́⁺∂qₖ[is+1:is+3,:]
                𝐜ᵀ[is+1:is+3,n1+1:n2] .= Dⁱₖ*∂vₖ∂λₘ
            end


            𝐜ᵀinv𝐉 = 𝐜ᵀ*inv(𝐉)
            𝐍 .= 𝐜ᵀinv𝐉*𝐁
            𝐫 .= (v́⁺ + 𝐛) - 𝐜ᵀinv𝐉*(𝐫𝐞𝐬 + 𝐁*𝚲ₖ)
        end
        # debug
        # @show norm(D*vₖ + 𝐛), norm(𝐫𝐞𝐬)
        # @show 𝚲ₖ, D*vₖ, 𝐛
        # @show 𝚲ₖ[1:3]⋅(D*vₖ + 𝐛)[1:3]

    end
    ns_stepk!
end

function solve!(intor::Integrator,solvercache::ZhongCCPCache;
                dt,ftol=1e-14,xtol=ftol,verbose=false,maxiters=50,
                progress=true,exception=true)
    (;prob,controller,tspan,restart,totalstep) = intor
    (;bot,dynfuncs) = prob
    (;traj,contacts_traj) = bot
    F!, Jac_F!, prepare_contacts!,get_directions_and_positions,get_∂Dq̇∂q,get_∂DᵀΛ∂q = dynfuncs
    (;cache) = solvercache
    (;M,Φ,A,Ψ,B,∂Ψ∂q,∂Aᵀλ∂q,∂Bᵀμ∂q) = cache
    invM = inv(M)
    q0 = traj.q[begin]
    λ0 = traj.λ[begin]
    q̇0 = traj.q̇[begin]
    pₖ₋₁ = M*q̇0
    pₖ   = zero(pₖ₋₁)
    T = eltype(q0)
    nq = length(q0)
    nλ = length(λ0)
    ∂F∂q = zeros(T,nq,nq)
    ∂F∂q̇ = zeros(T,nq,nq)
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
        active_contacts,gaps,H,es = prepare_contacts!(cₖ,qˣ)
        na = length(active_contacts)
        D,_ = get_directions_and_positions(active_contacts,qˣ)        
        persistent_indices = findall((c)->c.state.persistent,active_contacts)
        Dₘ = zero(D)
        Dₖ = copy(D)
        # Dₘ = copy(D)
        # Dₖ = zero(D)
        filtered_gaps = zero(gaps)
        if (na !== 0) && !isempty(persistent_indices)
            epi = reduce(vcat,[collect(3(i-1)+1:3i) for i in persistent_indices])
            Dₘ[epi,:] .= D[epi,:]
            Dₖ[epi,:] .= 0
            # filtered_gaps[persistent_indices] = gaps[persistent_indices]
        end
        isconverged = false
        normRes = typemax(T)
        iteration_break = 0
        x[      1:nq]          .= qₖ
        x[   nq+1:nq+nλ]       .= 0.0
        isconverged = false
        nΛ = 3na
        𝚲ₖ = zeros(T,nΛ)
        𝚲ₖ .= repeat([0.1,0,0],na)
        𝚲ʳₖ = copy(𝚲ₖ)
        Δ𝚲ₖ = copy(𝚲ₖ)
        𝐁 = zeros(T,nx,nΛ)
        𝐛 = zeros(T,nΛ)
        𝐜ᵀ = zeros(T,nΛ,nx)
        𝐍 = zeros(T,nΛ,nΛ)
        𝐫 = zeros(T,nΛ)
        scalingΛ = dt
        ns_stepk! = make_zhongccp_ns_stepk(nq,nλ,na,qₖ₋₁,q̇ₖ₋₁,pₖ₋₁,tₖ₋₁,pₖ,q̇ₖ,dynfuncs,cache,invM,dt,scalingΛ,persistent_indices)

        for iteration = 1:maxiters
            # @show iteration,D,ηs,es,gaps
            ns_stepk!(Res,Jac,𝐁,𝐛,𝐜ᵀ,𝐍,𝐫,x,𝚲ₖ,D,Dₘ,Dₖ,H,filtered_gaps,es,timestep,iteration)
            normRes = norm(Res)
            if na == 0
                if normRes < ftol
                    isconverged = true
                    iteration_break = iteration-1
                    break
                end
                Δx .= -Jac\Res
                x .+= Δx
            else # na!=0
                if iteration < 4
                    Nmax = 50
                else
                    Nmax = 50
                end
                𝚲ₖini = deepcopy(𝚲ₖ)
                𝚲ₖini[begin+1:3:end] .= 0.0
                𝚲ₖini[begin+2:3:end] .= 0.0
                𝚲ₖini .*= 10
                yini = deepcopy(𝚲ₖini)
                IPM!(𝚲ₖ,na,nΛ,𝚲ₖini,yini,𝐍,𝐫;ftol=1e-14,Nmax)
                
                Δ𝚲ₖ .= 𝚲ₖ - 𝚲ʳₖ
                minusRes𝚲 = -Res + 𝐁*(Δ𝚲ₖ)
                normRes = norm(minusRes𝚲)
                if  normRes < ftol
                    isconverged = true
                    iteration_break = iteration-1
                    break
                end
                Δx .= Jac\minusRes𝚲
                𝚲ʳₖ .= 𝚲ₖ
                x .+= Δx
                # @show timestep, iteration, normRes, norm(Δx), norm(Δ𝚲ₖ),persistent_indices
            end
        end
        qₖ .= x[      1:nq]
        λₘ .= x[   nq+1:nq+nλ]
        pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,M,A,dt)
        q̇ₖ .= invM*pₖ
        if na != 0
            update_contacts!(active_contacts,Dₘ*(qₖ.-qₖ₋₁).+Dₖ*q̇ₖ,𝚲ₖ./scalingΛ)
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
        if verbose
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
