struct ZhongQCCPNMono <: AbstractSolver end

struct ZhongQCCPNMonoCache{CacheType}
    cache::CacheType
end

function generate_cache(::ZhongQCCPNMono,intor;dt,kargs...)
    (;structure) = intor.prob.bot
    Mₘ = assemble_M(structure) 
    M⁻¹ₘ = assemble_M⁻¹(structure)
    M⁻¹ₖ = deepcopy(M⁻¹ₘ)
    ∂Mₘq̇ₘ∂qₘ = assemble_∂Mq̇∂q(structure)
    ∂Mₘqₖ∂qₘ = zero(∂Mₘq̇ₘ∂qₘ)
    ∂M⁻¹ₖpₖ∂qₖ = assemble_∂M⁻¹p∂q(structure)
    M! = make_M!(structure)
    M⁻¹! = make_M⁻¹!(structure)
    Jac_M! = make_Jac_M!(structure)
    Jac_M⁻¹! = make_Jac_M⁻¹!(structure)
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
    ∂Aᵀλ∂q(q,λ) = cstr_forces_jacobian(structure,λ)
    # ∂𝚽𝐪𝐯∂𝒒(q,v) = RB.∂Aq̇∂q(st,v)
    ∂Bᵀμ∂q(q,μ) = zeros(T,nq,nq)
    cache = @eponymtuple(
        Mₘ,M⁻¹ₘ,M⁻¹ₖ,
        ∂Mₘq̇ₘ∂qₘ,∂Mₘqₖ∂qₘ,∂M⁻¹ₖpₖ∂qₖ,
        Fₘ,∂Fₘ∂qₘ,∂Fₘ∂q̇ₘ,
        M!,Jac_M!,
        M⁻¹!,Jac_M⁻¹!,
        Φ,A,Ψ,B,
        ∂Ψ∂q,
        ∂Aᵀλ∂q,∂Bᵀμ∂q,
    )
    ZhongQCCPNMonoCache(cache)
end

function Momentum_ZhongQCCPNMono_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,Λₘ,Dₖ₋₁,Dₖ,H,scaling,h)
    pₖ = -pₖ₋₁ .+ 
        2/h.*Mₘ*(qₖ.-qₖ₋₁) .+ 
        scaling/h.*(transpose(A(qₖ))-transpose(A(qₖ₋₁)))*λₘ .+
        scaling.*(transpose(Dₖ)-transpose(Dₖ₋₁))*H*Λₘ
end

function make_zhongccpn_mono_ns_stepk(
        nq,nλ,na,
        qₖ₋₁,vₖ₋₁,pₖ₋₁,tₖ₋₁,
        pₖ,vₖ,
        F!,Jac_F!,
        get_directions_and_positions!,
        cache,
        h,scaling,
        persistent_idx,bodyid2act_idx
    )
    (;
        Mₘ,M⁻¹ₘ,M⁻¹ₖ,
        ∂Mₘq̇ₘ∂qₘ,∂Mₘqₖ∂qₘ,∂M⁻¹ₖpₖ∂qₖ,
        Fₘ,∂Fₘ∂qₘ,∂Fₘ∂q̇ₘ,
        M!,Jac_M!,
        M⁻¹!,Jac_M⁻¹!,
        Φ,A,
        ∂Aᵀλ∂q,
    ) = cache
    # T = eltype(qₖ₋₁)
    n1 = nq
    n2 = nq+nλ
    nx = n2+2na
    
    function ns_stepk!(𝐫𝐞𝐬,𝐉,𝐰,
            x,Λₘ,y,∂y∂x,
            Λ_split,y_split,
            Dₖ₋₁,ŕₖ₋₁,
            Dₖ,Dper, Dimp, ∂Dₖvₖ∂qₖ, ∂DᵀₖHΛₘ∂qₖ, ŕₖ,H,
            restitution_coefficients,timestep,iteration)
        # @show timestep, iteration, na, persistent_idx
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
            pₖ .= Momentum_ZhongQCCPNMono_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,Λₘ,Dₖ₋₁,Dₖ,H,scaling,h)
            M⁻¹!(M⁻¹ₖ,qₖ) 
            vₖ .= M⁻¹ₖ*pₖ
            M⁻¹!(M⁻¹ₘ,qₘ)
            Jac_M!(∂Mₘq̇ₘ∂qₘ,qₘ,q̇ₘ)
            Jac_M⁻¹!(∂M⁻¹ₖpₖ∂qₖ,qₖ,pₖ)
            ∂Aᵀₖλₘ∂qₖ = ∂Aᵀλ∂q(qₖ,λₘ)
            get_directions_and_positions!(Dₖ,Dper, Dimp, ∂Dₖvₖ∂qₖ, ∂DᵀₖHΛₘ∂qₖ,ŕₖ,qₖ, vₖ, H*Λₘ,bodyid2act_idx)
            ∂pₖ∂qₖ = 2/h.*Mₘ + 
                    ∂Mₘq̇ₘ∂qₘ .+
                    scaling/(h).*∂Aᵀₖλₘ∂qₖ .+ 
                    scaling.*∂DᵀₖHΛₘ∂qₖ
            ∂vₖ∂qₖ = M⁻¹ₖ*∂pₖ∂qₖ .+ ∂M⁻¹ₖpₖ∂qₖ
            ∂vₖ∂λₘ = scaling/h.*M⁻¹ₘ*transpose(Aₖ-Aₖ₋₁)
            
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
                    ∂y∂x[i,n1+1:n2] .= Dⁱₖ*∂vₖ∂λₘ
                end
            end
            𝐫𝐞𝐬[(n2   +1):(n2+ na)] .= v́⁺ .+ 𝐰 .- y
            𝐫𝐞𝐬[n2+na+1:n2+2na]     .= reduce(vcat,Λ_split⊙y_split)

            𝐉[      1:n1    , n2+   1:n2+ na] .=  -scaling.*h .*transpose(Dₖ₋₁)*H
            𝐉[n2+1:n2+ na,      1:n2    ]     .=  ∂y∂x
            𝐉[n2+1:n2+ na, n2+na+1:n2+2na]    .= -I(na)
            𝐉[n2+na+1:n2+2na, n2+   1:n2+ na] .=  BlockDiagonal(mat.(y_split))
            𝐉[n2+na+1:n2+2na, n2+na+1:n2+2na] .=  BlockDiagonal(mat.(Λ_split))
        end
    end
    ns_stepk!
end

function solve!(intor::Integrator,solvercache::ZhongQCCPNMonoCache;
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
    (;Mₘ,M⁻¹ₘ,M!,M⁻¹!,A) = cache
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
    prepare_contacts!(q0)
    mr = norm(Mₘ,Inf)
    scaling = mr
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
        na,bodyid2act_idx,persistent_idx,contacts_bits,
        H,restitution_coefficients,Dₖ₋₁, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕₖ₋₁, 
        L = prepare_contacts!(qₖ₋½)
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

        get_directions_and_positions!(Dₖ₋₁, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕₖ₋₁, qₖ₋₁, q̇ₖ₋₁, Λₘ, bodyid2act_idx,)
        Dₖ = deepcopy(Dₖ₋₁)
        ŕₖ = deepcopy(ŕₖ₋₁)
        ns_stepk! = make_zhongccpn_mono_ns_stepk(
            nq,nλ,na,qₖ₋₁,q̇ₖ₋₁,pₖ₋₁,tₖ₋₁,pₖ,q̇ₖ,
            F!,Jac_F!,
            get_directions_and_positions!,
            cache,dt,scaling,persistent_idx,bodyid2act_idx
        )
        restart_count = 0
        Λ_guess = 0.1
        while restart_count < 10
            Λₘ .= repeat([Λ_guess],na)
            y .= Λₘ
            x[      1:nq]          .= qₖ
            x[   nq+1:nq+nλ]       .= 0.0
            Nmax = 50
            for iteration = 1:maxiters
                # @show iteration,D,ηs,restitution_coefficients,gaps
                ns_stepk!(Res,Jac,
                    𝐰,x,Λₘ,y,∂y∂x,
                    Λ_split,y_split,
                    Dₖ₋₁,ŕₖ₋₁,Dₖ,Dper,Dimp,∂Dq̇∂q,∂DᵀΛ∂q,ŕₖ,H,
                    restitution_coefficients,timestep,iteration
                )
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
                    Δxc .= lu𝐉\(-Res)
                    # η = exp(-0.1μ) + 0.9
                    α_Λ = find_nonnegative_step_length(Λ_split,ΔΛc_split)
                    # @show Λ_split,ΔΛc_split
                    α_y = find_nonnegative_step_length(y_split,Δyc_split)
                    αmax = min(α_Λ,α_y)
                    α = min(1,0.99αmax)
                    # Λₘ .+= α.*ΔΛc
                    # y .+= α.*Δyc
                    x .+= α.*Δxc
                    μ = transpose(y)*Λₘ/na
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
        get_directions_and_positions!(Dₖ, Dper, Dimp, ∂Dq̇∂q, ∂DᵀΛ∂q, ŕₖ, qₖ, q̇ₖ, Λₘ, bodyid2act_idx)
        pₖ .= Momentum_ZhongQCCPNMono_k(qₖ₋₁,pₖ₋₁,qₖ,λₘ,Mₘ,A,Λₘ,Dₖ₋₁,Dₖ,H,scaling,dt)
        M⁻¹!(M⁻¹ₘ,qₖ)
        q̇ₖ .= M⁻¹ₘ*pₖ
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
