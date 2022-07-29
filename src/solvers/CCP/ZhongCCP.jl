struct ZhongCCP end

struct ZhongCCPCache{CacheType}
    cache::CacheType
end

function generate_cache(::ZhongCCP,intor;dt,kargs...)
    (;prob,state) = intor
    (;bot,dynfuncs) = prob
    (;tg) = bot
    (;q,q̇) = state.now
    # F!,_ = dynfuncs
    # mm = TR.build_MassMatrices(bot)
    M = Matrix(build_M(tg))
    # (;M) = mm
    A = make_A(bot)
    Φ = make_Φ(bot)

    nq = length(q)
    T = eltype(q)
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

function Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₖ,M,A,h)
    pₖ = -pₖ₋₁ .+ 2/h.*M*(qₖ.-qₖ₋₁) .+ 1/(h).*(transpose(A(qₖ))-transpose(A(qₖ₋₁)))*λₖ
       # .+ 1/(h).*(transpose(B(qₖ))-transpose(B(qₖ₋₁)))*μₖ
end

function make_zhongccp_ns_stepk(nq,nλ,na,qₖ₋₁,vₖ₋₁,pₖ₋₁,tₖ₋₁,pₖ,vₖ,dynfuncs,cache,invM,h,scalingΛ,persistent_indices)
    F!,Jac_F!,_ = dynfuncs
    (;M,Φ,A,Ψ,B,∂Ψ∂q,∂Aᵀλ∂q,∂Bᵀμ∂q) = cache

    n1 = nq
    n2 = nq+nλ
    nΛ = 3na
    nx = n2
    function ns_stepk!(𝐫𝐞𝐬,𝐉,𝐁,𝐛,𝐜ᵀ,𝐍,𝐫,x,𝚲ₖ,D,Dₘ,Dₖ,H,filtered_gaps,es,timestep,iteration)
        # @show timestep, iteration, na, persistent_indices
        qₖ = @view x[   1:n1]
        λₖ = @view x[n1+1:n2]

        q = (qₖ.+qₖ₋₁)./2
        q̇ = (qₖ.-qₖ₋₁)./h
        vₘ = q̇
        t = tₖ₋₁+h/2
        T = eltype(qₖ)
        F⁺ = zeros(T,nq)
        F!(F⁺,q,q̇,t)
        ∂F∂q = zeros(T,nq,nq)
        ∂F∂q̇ = zeros(T,nq,nq)
        Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)

        Aₖ₋₁ = A(qₖ₋₁)
        Aₖ   = A(qₖ)
        # Bₖ₋₁ = B(qₖ₋₁)
        # Bₖ = B(qₖ)

        ∂vₘ∂qₖ = 1/h*I
        # ∂v1∂λₖ = 0

        ∂vₖ∂qₖ = 2/h*I + 1/(h).*invM*(∂Aᵀλ∂q(qₖ,λₖ)) #+ ∂Bᵀμ∂q(qₖ,μₖ))
        ∂vₖ∂λₖ = invM*transpose(Aₖ-Aₖ₋₁)/(h)
        # ∂vₖ∂μₖ = invM*transpose(Bₖ-Bₖ₋₁)/(h)

        ∂DᵀHΛₖ∂qₖ = zeros(T,nq,nq)

        𝐫𝐞𝐬[   1:n1] .= -h.*pₖ₋₁ .+ M*(qₖ.-qₖ₋₁) .-
                        scalingΛ .*transpose(D)*H*𝚲ₖ .-
                        transpose(Aₖ₋₁)*λₖ .-
                        # transpose(Bₖ₋₁)*μₖ .-
                        (h^2)/2 .*F⁺

        𝐫𝐞𝐬[n1+1:n2] .= Φ(qₖ)
        # 𝐫𝐞𝐬[n2+1:n3] .= Ψ(qₖ,vₖ)

        # ∂Aᵀλ∂q⁺ = ∂Aᵀλ∂q(q,λₖ)
        𝐉 .= 0.0
        𝐉[   1:n1,   1:n1] .=  M .- scalingΛ .*∂DᵀHΛₖ∂qₖ .-h^2/2 .*(1/2 .*∂F∂q .+ 1/h.*∂F∂q̇)
        𝐉[   1:n1,n1+1:n2] .= -transpose(Aₖ₋₁)
        # 𝐉[   1:n1,n2+1:n3] .= -transpose(Bₖ₋₁)

        𝐉[n1+1:n2,   1:n1] .=  Aₖ

        # 𝐉[n2+1:n3,   1:n1] .=  ∂Ψ∂q(qₖ,vₖ) .+ Bₖ*∂vₖ∂qₖ
        # 𝐉[n2+1:n3,n1+1:n2] .=  Bₖ*∂vₖ∂λₖ
        # 𝐉[n2+1:n3,n2+1:n3] .=  Bₖ*∂vₖ∂μₖ

        𝐁 .= 0
        𝐁[   1:n1,1:nΛ] .= scalingΛ .*transpose(D)*H

        # ∂Dvₖ∂qₖ = zeros(T,nΛ,nq)
        # ∂bₖ∂qₖ = zeros(T,nΛ,nq)
        # v́ₖ₋₁ = D*vₖ₋₁
        # v́ₖ   = D*vₖ
        # S = zeros(Int,2na,3na)
        # for i = 1:na
        #     for j = 1:2
        #         S[2(i-1)+j,3(i-1)+1+j] = 1
        #     end
        # end
        # v́ₜ = S*v́ₖ

        𝐜ᵀ .= 0
        if na != 0
            v́ₖ₋₁ = Dₖ*vₖ₋₁
            function make_b(Dₘ,Dₖ,qₖ₋₁,h)
                function inner_b(q)
                    # @show D*vₖ
                    vₘ .= (q.-qₖ₋₁)./h
                    pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,q,λₖ,M,A,h)
                    vₖ .= invM*pₖ
                    v́ = Dₘ*vₘ .+ Dₖ*vₖ
                    v́ₜ = [
                        begin
                            is = 3(i-1)
                            vⁱₖ₋₁ = @view v́ₖ₋₁[is+1:is+3]
                            vⁱ   = @view v́[is+1:is+3]
                            vₜⁱₖ₋₁ = norm(vⁱₖ₋₁[2:3])
                            vₜⁱ   = norm(vⁱ[2:3])
                            vₙⁱₖ₋₁ = vⁱₖ₋₁[1]
                            vₙⁱ   = vⁱ[1]
                            # @show timestep,iteration, vₙⁱₖ₋₁, vₙⁱ, vₜⁱₖ₋₁, vₜⁱ, 𝚲ₖ
                            vₜⁱ + es[i]*min(vₙⁱₖ₋₁,zero(T))
                        end
                        for i = 1:na
                    ]
                    reduce(
                        vcat,[
                            [v́ₜ[i]+filtered_gaps[i],0,0]
                            for i = 1:na
                        ]
                    )
                end
                function inner_b_Jac(q)
                    FiniteDiff.finite_difference_jacobian(inner_b,q,Val{:central})
                end
                inner_b, inner_b_Jac
            end
            𝒃, 𝒋 = make_b(Dₘ,Dₖ,qₖ₋₁,h)
            𝐛 .= 𝒃(qₖ)
            𝐜ᵀ[:,1:n1].= 𝒋(qₖ)
            # for i = 1:na
            #
            #     𝐛[is+1] = vₜⁱₖ₋₁ + es[i]*min(vₙⁱₖ₋₁,zero(T)) #+ gaps[i]/h
            # end

            for i = 1:na
                is = 3(i-1)
                D̃ₘⁱ = @view Dₘ[is+1:is+3,:]
                D̃ₖⁱ = @view Dₖ[is+1:is+3,:]
                𝐜ᵀ[is+1:is+3,   1:n1] .+= D̃ₘⁱ*∂vₘ∂qₖ .+ D̃ₖⁱ*∂vₖ∂qₖ
                𝐜ᵀ[is+1:is+3,n1+1:n2] .= D̃ₖⁱ*∂vₖ∂λₖ
                # 𝐜ᵀ[is+1:is+3,n2+1:n3] .= D̃i*∂vₖ∂μₖ
            end


            𝐜ᵀinv𝐉 = 𝐜ᵀ*inv(𝐉)
            𝐍 .= 𝐜ᵀinv𝐉*𝐁
            vₘ .= (qₖ.-qₖ₋₁)./h
            pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₖ,M,A,h)
            vₖ .= invM*pₖ
            𝐫 .= (Dₖ*vₖ + Dₘ*vₘ + 𝐛) - 𝐜ᵀinv𝐉*(𝐫𝐞𝐬 + 𝐁*𝚲ₖ)
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
    (;prob,state,control!,tspan,restart,totalstep) = intor
    (;bot,dynfuncs) = prob
    (;traj,contacts_traj) = bot
    # @unpack t,q,q̇,tprev,qprev,q̇prev = state
    F!, Jac_F!, prepare_contacts! = dynfuncs
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
    step = 0
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
        #---------Step k Control-----------
        # control!(intor,cache)
        #---------Step k Control-----------
        push!(contacts_traj,deepcopy(contacts_traj[end]))
        cₖ = contacts_traj[timestep+1]
        qₖ₋₁ = traj.q[timestep]
        q̇ₖ₋₁ = traj.q̇[timestep]
        # pₖ₋₁ = traj.p[timestep]
        λₖ₋₁ = traj.λ[timestep]
        tₖ₋₁ = traj.t[timestep]
        qₖ   = traj.q[timestep+1]
        q̇ₖ   = traj.q̇[timestep+1]
        # pₖ   = traj.p[timestep+1]
        λₖ   = traj.λ[timestep+1]
        pₖ₋₁ = M*q̇ₖ₋₁
        qˣ = qₖ₋₁ .+ dt./2 .*q̇ₖ₋₁
        qₖ .= qₖ₋₁ .+ dt .*q̇ₖ₋₁
        q̇ₖ .= q̇ₖ₋₁
        active_contacts,na,gaps,D,H,es = prepare_contacts!(cₖ,qˣ)
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
        # x[nq+nλ+1:nq+nλ+nμ]    .= 0.0
        isconverged = false
        nΛ = 3na
        𝚲ₖ = zeros(T,nΛ)
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
                # r4 = make_residual4(ηs,𝐍,𝐫;gd=1e-3)
                # Jacobi_B = make_B(na,D,invM)
                # 𝚲ₖ,_ = Jacobi(Jacobi_B,r4,ηs,𝐍,𝐫;τ=1e-10,Nmax=1000)
                # 𝚲uₖ₊₁,GS_k,GS_res = GaussSeidel(u,B,r,ηs,𝐍,𝐫)
                IPM!(𝚲ₖ,na,nΛ,repeat([0.1,0,0],na),repeat([0.1,0,0],na),𝐍,𝐫;ftol=1e-14,Nmax=50)
                # @show iteration, 𝚲ₖ, y_split[1]
                # y = 𝐍*𝚲ₖ+𝐫
                # @show timestep, iteration, 𝐫, 𝚲ₖ, y, 𝚲ₖ⋅y
                # APGD_res = APGD!(𝚲ₖ,r4,ηs,𝐍,𝐫;τ=1e-10,Nmax=1000)
                # @show APGD_res
                Δ𝚲ₖ .= 𝚲ₖ - 𝚲ʳₖ
                # @show "GD", 𝚲ₖ#, 𝚲ʳₖ
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
                # normΔx = norm(Δx)
                # res = normΔx
                # @show timestep, iteration, normRes, normΔx, norm(Δ𝚲ₖ)
                # iteration_break = iteration
                # @show timestep, iteration, (D*q̇ₖ)[1:3]
            end
        end
        qₖ .= x[      1:nq]
        λₖ .= x[   nq+1:nq+nλ]
        pₖ .= Momentum_k(qₖ₋₁,pₖ₋₁,qₖ,λₖ,M,A,dt)
        q̇ₖ .= invM*pₖ
        if na != 0
            # @show [ac.id for ac in active_contacts]
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
        else
            if na == 0
                stepstring = "Smooth"
                # @info "$stepstring timestep=$timestep, iterations=$iteration_break"
            else
                stepstring = "Nonsmooth with $na contact(s)"
                # @info "$stepstring timestep=$timestep, iterations=$iteration_break, and APGD_res=$APGD_res"
            end
        end

        #---------Step k finisher-----------
        step += 1
        #---------Step k finisher-----------
        if verbose
            @printf("Progress: %5.1f%%, step: %s, time: %s, iterations: %s \n", (
            ts[timestep]/totaltime*100), timestep, ts[timestep], R_stepk_result.iterations)
        end
        next!(prog)
    end
    bot
end
