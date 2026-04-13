

struct AlphaCCPCache{CacheType}
    cache::CacheType
end

function generate_cache(
        simulator::Simulator{DynamicsProblem{
            RobotType,
            RestitutionFrictionCombined{NewtonRestitution,CoulombFriction}
        }},
        solver::DynamicsSolver{
            GeneralizedAlpha,
            <:AbstractBodySolver,
            <:AbstractApparatusSolver,
            InnerLayerContactSolver
        };
        dt,kargs...
    )   where RobotType
    (;prob,state) = sim
    (;bot,dynfuncs) = prob
    (;st) = bot
    (;q,q̇) = state.now
    (;ρ∞) = solver
    coeffs = generalized_alpha(ρ∞,dt)
    # F!,_ = dynfuncs
    # mm = RB.build_mass_matrices(bot)
    M = Matrix(assemble_M(st))
    # (;M) = mm
    A = make_cstr_jacobian(bot)
    Φ = make_cstr_function(bot)

    nq = length(q)
    T = eltype(q)
    Ψ(q,q̇) = Vector{T}()
    ∂Ψ∂q(q,q̇) = Matrix{T}(undef,0,nq)
    B(q) = Matrix{T}(undef,0,nq)

    # ∂𝐌𝐚∂𝐪(q,a) = zeros(T,nq,nq)
    cstr_forces_jacobian(q,λ) = cstr_forces_jacobian(st,λ)
    # ∂𝚽𝐪𝐯∂𝒒(q,v) = RB.∂Aq̇∂q(st,v)
    ∂Bᵀμ∂q(q,μ) = zeros(T,nq,nq)
    cache = @eponymtuple(M,Φ,A,Ψ,B,∂Ψ∂q,cstr_forces_jacobian,∂Bᵀμ∂q,coeffs)
    AlphaCCPCache(cache)
end

function update_nonsmooth!(Res,Jac,nq,nλ,na,xe,vₛ,gₙ,Dₛ₊₁,H,restitution_coefficients,t,p,h,scaling,dynfuncs,cache)
    qₛ₊₁, vₛ₊₁, ṽₛ₊₁, ṽ̇ₛ₊₁, 𝛌bₛ₊₁, 𝚲uₛ₊₁ = xe
    F!,Jac_F!,_ = dynfuncs
    (;M,Φ,A,Ψ,B,∂Ψ∂q,cstr_forces_jacobian,∂Bᵀμ∂q) = cache
    # 𝐠,get_idx,get_D = contact_funcs

    (;αm,αf,γ,β,βₜ,γₜ) = p
    T = eltype(qₛ₊₁)
    v̂ₛ₊₁ = vₛ₊₁ - ṽₛ₊₁
    Mₛ₊₁ = M
    Fₛ₊₁ = copy(qₛ₊₁)
    F!(Fₛ₊₁,qₛ₊₁,vₛ₊₁,t)
    Fqₛ₊₁ = copy(Mₛ₊₁)
    Fvₛ₊₁ = copy(Mₛ₊₁)
    Jac_F!(Fqₛ₊₁,Fvₛ₊₁,qₛ₊₁,vₛ₊₁,t)
    Φₛ₊₁ = Φ(qₛ₊₁)
    Φqₛ₊₁ = A(qₛ₊₁)
    cstr_forces_jacobianₛ₊₁ = cstr_forces_jacobian(qₛ₊₁,𝛌bₛ₊₁)
    # D̃ₛ₊₁ = D̃(qₛ₊₁)
    # D̃qₛ₊₁ = D̃q(qₛ₊₁)
    # D̂qₛ₊₁ = D̂q(qₛ₊₁)
    # Dₛ₊₁ = vcat(D̃ₛ₊₁,D̂ₛ₊₁)
    # Dqₛ₊₁ = vcat(D̃qₛ₊₁,D̂qₛ₊₁)
    nΛ = 3na
    ∂Dv∂qₛ₊₁ = zeros(T,nΛ,nq)
    h² = h^2
    Res .= vcat(
        Mₛ₊₁*ṽ̇ₛ₊₁ .+ scaling.*transpose(Φqₛ₊₁)*𝛌bₛ₊₁ .- Fₛ₊₁,
        Mₛ₊₁*v̂ₛ₊₁- transpose(Dₛ₊₁)*H*𝚲uₛ₊₁,
        scaling.*Φₛ₊₁
    )

    Kₜ = scaling.*cstr_forces_jacobianₛ₊₁ .- Fqₛ₊₁
    Jac .= [
        γₜ.*Mₛ₊₁.+βₜ.*Kₜ (h/2).*(Kₜ.-Fvₛ₊₁) scaling.*transpose(Φqₛ₊₁);
         -Mₛ₊₁          Mₛ₊₁      zeros(T,nq,nλ);
         scaling.*βₜ.*Φqₛ₊₁        scaling.*(h/2).*Φqₛ₊₁ 0I
    ]
    𝐎nu = zeros(T,nq,nΛ)
    𝐎bu = zeros(T,nλ,nΛ)
    𝐁 = [
        𝐎nu;
        transpose(Dₛ₊₁)*H;
        𝐎bu
    ]
    𝐜ᵀ = [
        βₜ.*∂Dv∂qₛ₊₁ (h/2).*∂Dv∂qₛ₊₁ .+ Dₛ₊₁ zeros(T,nΛ,nλ);
    ]
    𝐜ᵀinv𝐉 = 𝐜ᵀ*inv(Jac)
    𝐍 = 𝐜ᵀinv𝐉*𝐁
    𝐛 = zeros(T,nΛ)
    v́ₛ = Dₛ₊₁*vₛ
    v́ₛ₊₁ = Dₛ₊₁*vₛ₊₁
    for i = 1:na
        is = 3(i-1)
        vⁱₛ = v́ₛ[is+1:is+3]
        vⁱₛ₊₁ = v́ₛ₊₁[is+1:is+3]
        vₜⁱₛ₊₁ = norm(vⁱₛ₊₁[2:3])
        vₙⁱₛ = vⁱₛ[1]
        # @show i, vₙⁱₛ, vₜⁱₛ₊₁
        # @show vₛ, vⁱₛ
        # @show Dₛ₊₁
        # @show i,friction_coefficients[i],restitution_coefficients[i],vₙⁱₛ
        # @show i, 𝚲uₛ₊₁[is+1]
        𝐛[is+1] = vₜⁱₛ₊₁ + restitution_coefficients[i]*vₙⁱₛ
    end
    𝐫 = (Dₛ₊₁*vₛ₊₁ + 𝐛) - 𝐜ᵀinv𝐉*(Res + 𝐁*𝚲uₛ₊₁)
    # @show Dₛ₊₁*vₛ₊₁

    𝐁,𝐜ᵀ,𝐍,𝐫
end

function solve!(sim::Simulator,solvercache::AlphaCCPCache;
                dt,ftol=1e-14,xtol=ftol,verbose=false,maxiters=50,
                progress=true,exception=true)
    (;prob,state,controller,tspan,restart,totalstep) = sim
    (;bot,dynfuncs) = prob
    (;traj) = bot
    # (;t,q,q̇,tprev,qprev,q̇prev) = state
    F!, Jac_F!, (contacts, prepare_contacts!, update_contacts!) = dynfuncs
    (;cache) = solvercache
    (;M,Φ,A,Ψ,B,∂Ψ∂q,cstr_forces_jacobian,∂Bᵀμ∂q,coeffs) = cache
    (;αm,αf,γ,β,γₜ,βₜ)= coeffs
    contacts_traj = [deepcopy(contacts) for i in 1:totalstep]
    q0 = traj.q[begin]
    v0 = traj.q̇[begin]
    λ0 = traj.λ[begin]
    v0 = traj.q̇[begin]
    t0 = traj.t[begin]
    h = dt
    T = eltype(q0)
    invM = inv(M)
    # initial conditions
    M0 = M
    F0 = copy(q0)
    F!(F0,q0,v0,t0)
    ṽ̇0 = M0\F0
    a0 = ṽ̇0
    # allocations
    as = [copy(a0) for i = 1:totalstep+1]
    ṽ̇s = [copy(ṽ̇0) for i = 1:totalstep+1]
    ṽs = [copy(v0) for i = 1:totalstep+1]
    mr = norm(M0,Inf)
    scaling = mr/h
    qˣ = copy(q0)
    nq = length(q0)
    nλ = length(λ0)
    nx = nq+nq+nλ
    Δx = ones(T,nx)
    Res = zeros(T,nx)
    Jac = zeros(T,nx,nx)
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        cₛ₊₁ = contacts_traj[timestep]
        tₛ   = traj.t[timestep]
        qₛ   = traj.q[timestep]
        qₛ₊₁ = traj.q[timestep+1]
        vₛ   = traj.q̇[timestep]
        vₛ₊₁ = traj.q̇[timestep+1]
        𝛌bₛ₊₁ = traj.λ[timestep+1]
        tₛ₊₁ = traj.t[timestep+1]
        aₛ   = as[timestep]
        aₛ₊₁ = as[timestep+1]
        ṽ̇ₛ   = ṽ̇s[timestep]
        ṽ̇ₛ₊₁ = ṽ̇s[timestep+1]
        ṽₛ   = ṽs[timestep]
        ṽₛ₊₁ = ṽs[timestep+1]

        # initial_guesses
        aₛ₊₁ .= (αf*ṽ̇ₛ - αm*aₛ)/(1-αm)
        ṽₛ₊₁ .= vₛ + h*(1-γ)*aₛ + h*γ*aₛ₊₁
        qₛ₊₁ .= qₛ + h*vₛ + h^2*(0.5-β)*aₛ + h^2*β*aₛ₊₁
        vₛ₊₁ .= ṽₛ₊₁
        ṽ̇ₛ₊₁ .= 0
        𝛌bₛ₊₁ .= 0

        # contact detection
        qˣ .= qₛ₊₁
        active_contacts,na,gₙ,Dₛ₊₁,H,restitution_coefficients = prepare_contacts!(cₛ₊₁,qˣ)
        nΛ = 3na
        𝚲uₛ₊₁ = zeros(T,nΛ)
        𝚲uʳₛ₊₁ = copy(𝚲uₛ₊₁)
        # Newton iteration
        isconverged = false
        for iteration = 1:maxiters
            xe = (qₛ₊₁, vₛ₊₁, ṽₛ₊₁, ṽ̇ₛ₊₁, 𝛌bₛ₊₁, 𝚲uₛ₊₁)

            𝐁,𝐜ᵀ,𝐍,𝐫 = update_nonsmooth!(Res,Jac,nq,nλ,na,xe,vₛ,gₙ,Dₛ₊₁,H,restitution_coefficients,tₛ₊₁,coeffs,dt,scaling,dynfuncs,cache)

            normRes = norm(Res)
            # @show normRes
            if  normRes < ftol
                isconverged = true
                break
            end
            if na != 0
                # B = make_B(u,Dₛ₊₁,invM)
                # r4 = make_residual4(friction_coefficients,𝐍,𝐫)
                # 𝚲uₛ₊₁,_ = Jacobi(B,r,friction_coefficients,𝐍,𝐫;τ=1e-13,Nmax=1000)
                # 𝚲uₛ₊₁,GS_k,GS_res = GaussSeidel(u,B,r,friction_coefficients,𝐍,𝐫)
                # @show GS_k,GS_res
                # 𝚲uₛ₊₁,_ = APGD(r,friction_coefficients,𝐍,𝐫)
                # @show 𝚲uₛ₊₁, vₛ₊₁
                # 𝚲uₛ₊₁,_  = APGD(r,friction_coefficients,𝐍,𝐫;τ=1e-10,Nmax=1000)
                # APGD!(𝚲uₛ₊₁,r4,friction_coefficients,𝐍,𝐫;τ=1e-10,Nmax=1000)
                IPM!(𝚲uₛ₊₁,na,nΛ,repeat([0.01,0,0],na),repeat([0.01,0,0],na),𝐍,𝐫;ftol=1e-14,Nmax=50)

                # @show sum(𝚲uₛ₊₁/h)
                Δ𝚲uₛ₊₁ = 𝚲uₛ₊₁ - 𝚲uʳₛ₊₁
                # @show iteration,abs(𝚲uₛ₊₁[3]/𝚲uₛ₊₁[1])
                # @show iteration,𝚲uʳₛ₊₁
                # @show iteration,Δ𝚲uₛ₊₁
                𝚲uʳₛ₊₁ .= 𝚲uₛ₊₁
                Δx .= Jac\(-Res + 𝐁*(Δ𝚲uₛ₊₁))
            else
                Δx .= Jac\(-Res)
            end
            Δṽ,Δv,Δ𝛌b = split_by_lengths(Δx,[nq,nq,nλ])
            ṽₛ₊₁  .+= Δṽ
            vₛ₊₁  .+= Δv
            ṽ̇ₛ₊₁  .+= γₜ*Δṽ
            qₛ₊₁  .+= βₜ*Δṽ + h/2*Δv
            𝛌bₛ₊₁ .+= Δ𝛌b
            # @show ṽ̇ₛ₊₁,vₛ₊₁,ṽₛ₊₁
            # @show 𝚲uₛ₊₁
        end
        aₛ₊₁ .+= (1-αf)/(1-αm)*ṽ̇ₛ₊₁
        if !isconverged
            @warn "Newton max iterations $maxiters, at timestep=$timestep, Res=$(res)"
            if exception
                @error "Not converged!"
                break
            else
                # sim.convergence = false
                # break
            end
        end
        next!(prog)

    end
    contacts_traj
end
