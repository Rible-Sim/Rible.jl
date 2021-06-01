mutable struct ActiveSets{T}
    𝒞::T
    𝒰::T
    𝒰c::T
    𝒜::T
    𝒜c::T
    ℬ::T
    ℬc::T
end

function generalized_α(ρ∞)
    αm = (2ρ∞-1)/(ρ∞+1)
    αf = ρ∞/(ρ∞+1)
    γ = 1/2 + αf - αm
    β = 1/4*(γ+1/2)^2
    (αm=αm,αf=αf,γ=γ,β=β)
end

function Newmark(ρ∞)
    αf = αm = 0.0
    γ = 1/2
    β = 1/4
    (αm=αm,αf=αf,γ=γ,β=β)
end

function find_corrections(ΔW,Δṽ,ΔU,p,h)
    @unpack αm,αf,γ,β = p
    Δv = Δṽ + ΔW
    Δṽ̇ = (1-αm)/((1-αf)*γ*h)*Δṽ
    Δq = h*β/γ*Δṽ + ΔU
    Δv,Δṽ,Δq
end

function initial_guesses(c,ū,qₙ,vₙ,aₙ,ṽ̇ₙ,p,h)
    @unpack αm,αf,γ,β = p
    aₙ₊₁0 = (αf*ṽ̇ₙ - αm*aₙ)/(1-αm)
    qₙ₊₁0 = qₙ + h*vₙ + h^2*(0.5-β)*aₙ + h^2*β*aₙ₊₁0
    ṽₙ₊₁0 = vₙ + h*(1-γ)*aₙ + h*γ*aₙ₊₁0
    vₙ₊₁0 = ṽₙ₊₁0
    Uₙ₊₁0 = zero(qₙ₊₁0)
    Wₙ₊₁0 = zero(qₙ₊₁0)
    𝛌ₙ₊₁0 = zeros(eltype(qₙ₊₁0),ū)
    𝛎ₙ₊₁0 = zeros(eltype(qₙ₊₁0),c)
    𝚲ₙ₊₁0 = zeros(eltype(qₙ₊₁0),c)
    ṽ̇ₙ₊₁0 = zero(qₙ₊₁0)
    qₙ₊₁0,vₙ₊₁0,aₙ₊₁0,ṽₙ₊₁0,ṽ̇ₙ₊₁0,Uₙ₊₁0,Wₙ₊₁0,𝛌ₙ₊₁0,𝛎ₙ₊₁0,𝚲ₙ₊₁0
end

function get_residuals(qₙ,vₙ,xe,t,h,active_sets,𝐞,𝐌,𝒈,𝒈𝒒,𝐟)
    @unpack 𝒰,𝒰c,𝒜,𝒜c,ℬ,ℬc = active_sets
    qₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁ = xe

    Mₙ₊₁ = 𝐌(qₙ₊₁)
    𝐠ₙ₊₁ = 𝒈(qₙ₊₁)
    𝐠𝐪ₙ = 𝒈𝒒(qₙ)
    𝐠𝐪ₙ₊₁ = 𝒈𝒒(qₙ₊₁)

    𝒰c_indice = sort(collect(𝒰c))
    𝐠𝐪ₙ₊₁𝒰c = @view 𝐠𝐪ₙ₊₁[𝒰c_indice,:]
    𝐫ˢ = vcat(
        Mₙ₊₁*h*ṽ̇ₙ₊₁ - transpose(𝐠𝐪ₙ₊₁𝒰c)*h*𝛌ₙ₊₁ - h*𝐟(qₙ₊₁,vₙ₊₁,t),
        𝐠𝐪ₙ₊₁𝒰c*ṽₙ₊₁
    )

    𝒜_indice = sort(collect(𝒜))
    𝒜c_indice = sort(collect(𝒜c))
    𝐠𝐪ₙ₊₁𝒜 = @view 𝐠𝐪ₙ₊₁[𝒜_indice,:]
    𝐠ₙ₊₁𝒜 = @view 𝐠ₙ₊₁[𝒜_indice]
    𝐫ᵖ = 1/h*vcat(
        Mₙ₊₁*Uₙ₊₁ - transpose(𝐠𝐪ₙ₊₁𝒜)*𝛎ₙ₊₁[𝒜_indice],
        𝐠ₙ₊₁𝒜,
        𝛎ₙ₊₁[𝒜c_indice]
    )
    # @show 𝛎ₙ₊₁[𝒜c_indice]

    ℬ_indice = sort(collect(ℬ))
    ℬc_indice = sort(collect(ℬc))
    𝐠ₙ₊₁ℬ = @view 𝐠ₙ₊₁[ℬ_indice]
    𝐠𝐪ₙℬ = @view 𝐠𝐪ₙ[ℬ_indice,:]
    𝐠𝐪ₙ₊₁ℬ = @view 𝐠𝐪ₙ₊₁[ℬ_indice,:]
    𝐫ᵛ = vcat(
        Mₙ₊₁*Wₙ₊₁ - transpose(𝐠𝐪ₙ₊₁ℬ)*𝚲ₙ₊₁[ℬ_indice],
        𝐠𝐪ₙ₊₁ℬ*vₙ₊₁ + Diagonal(𝐞[ℬ_indice])*𝐠𝐪ₙℬ*vₙ,
        𝚲ₙ₊₁[ℬc_indice]
    )

    # @show 𝚲ₙ₊₁[ℬc_indice]
    𝐫ˢ,𝐫ᵖ,𝐫ᵛ
end

function initialize_active_sets(𝒞,𝒰c)
    𝒰 = setdiff(𝒞,𝒰c)
    𝒜 = copy(𝒰c)
    𝒜c = copy(𝒰)
    ℬ = copy(𝒰c)
    ℬc = copy(𝒰)
    ActiveSets(𝒞,𝒰,𝒰c,𝒜,𝒜c,ℬ,ℬc)
end

function update_active_sets!(active_sets,qₙ,vₙ,q̃ₙ₊₁,xe,𝐞,𝒈,𝒈𝒒,r)
    @unpack 𝒞,𝒰,𝒰c,𝒜,𝒜c,ℬ,ℬc = active_sets
    qₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁ = xe
    𝐠ₙ = 𝒈(qₙ)
    𝐠ₙ₊₁ = 𝒈(qₙ₊₁)
    𝐠̃ₙ₊₁ = 𝒈(q̃ₙ₊₁)
    𝐠𝐪ₙ = 𝒈𝒒(qₙ)
    𝐠𝐪ₙ₊₁ = 𝒈𝒒(qₙ₊₁)
    # @show 𝛎ₙ₊₁
    # @show qₙ
    # @show qₙ₊₁
    # @show 𝐠ₙ
    # @show 𝐠ₙ₊₁
    𝒜 = 𝒰c ∪ Set([j for j ∈ 𝒰 if 𝛎ₙ₊₁[j] - r*𝐠ₙ₊₁[j] ≥ 0])
    𝒜c = setdiff(𝒞,𝒜)
    ℬ = 𝒰c ∪ Set([j for j ∈ 𝒰 if 𝐠̃ₙ₊₁[j] ≤ 0 && 𝚲ₙ₊₁[j] - r*(𝐠𝐪ₙ₊₁[j,:]⋅vₙ₊₁+𝐞[j]*𝐠𝐪ₙ[j,:]⋅vₙ) ≥ 0])
    ℬc = setdiff(𝒞,ℬ)
    @pack! active_sets = 𝒞,𝒰,𝒰c,𝒜,𝒜c,ℬ,ℬc
    active_sets
end

function pad_λ(c,S,λ)
    ret = zeros(eltype(λ),c)
    ret[S] = λ
    ret
end

function compute_Sₜˢ(n,c,ū,xe,t,p,h,active_sets,dyfuncs)
    qₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁ = xe
    @unpack αm,αf,γ,β = p
    @unpack 𝒞,𝒰,𝒰c,𝒜,𝒜c,ℬ,ℬc = active_sets
    𝐌,𝒈,𝒈𝒒,𝐟,Jacobians = dyfuncs
    ∂𝒈𝒒T𝛌∂𝒒,∂𝒈𝒒𝐯∂𝒒,∂𝐌𝐚∂𝐪,∂𝐟∂𝐪,∂𝐟∂𝐯 = Jacobians

    𝒰c_indice = sort(collect(𝒰c))

    Mₙ₊₁ = 𝐌(qₙ₊₁)

    𝛌 = pad_λ(c,𝒰c_indice,𝛌ₙ₊₁)

    Kₜ = ∂𝐌𝐚∂𝐪(qₙ₊₁,ṽ̇ₙ₊₁) - ∂𝐟∂𝐪(qₙ₊₁,vₙ₊₁,t) - ∂𝒈𝒒T𝛌∂𝒒(qₙ₊₁,𝛌)
    Cₜ = -∂𝐟∂𝐯(qₙ₊₁,vₙ₊₁,t)
    Sₜˣ = (1-αm)/(1-αf)/γ*Mₙ₊₁ + h*Cₜ + β*h^2/γ*Kₜ

    𝐠𝐪ₙ₊₁ = 𝒈𝒒(qₙ₊₁)
    𝐠𝐪ₙ₊₁𝒰c = 𝐠𝐪ₙ₊₁[𝒰c_indice,:]


    Gˣˢ = ∂𝒈𝒒𝐯∂𝒒(qₙ₊₁,ṽₙ₊₁)[𝒰c_indice,:]

    # 𝐎 = zeros(eltype(qₙ₊₁),ū,ū)
    Sₜˢ = [
        Sₜˣ transpose(𝐠𝐪ₙ₊₁𝒰c);
        𝐠𝐪ₙ₊₁𝒰c+h*β/γ*Gˣˢ 0I
    ]
end

function compute_Sₜᵖ(n,c,ū,xe,t,p,h,active_sets,dyfuncs)
    qₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁ = xe
    @unpack αm,αf,γ,β = p
    @unpack 𝒞,𝒰,𝒰c,𝒜,𝒜c,ℬ,ℬc = active_sets
    𝐌,𝒈,𝒈𝒒,𝐟,Jacobians = dyfuncs
    ∂𝒈𝒒T𝛌∂𝒒,∂𝒈𝒒𝐯∂𝒒,∂𝐌𝐚∂𝐪,∂𝐟∂𝐪,∂𝐟∂𝐯 = Jacobians
    𝒜_indice = sort(collect(𝒜))
    𝒜c_indice = sort(collect(𝒜c))
    ā = length(𝒜c_indice)
    Mₙ₊₁ = 𝐌(qₙ₊₁)
    𝐠𝐪ₙ₊₁ = 𝒈𝒒(qₙ₊₁)
    𝐠𝐪ₙ₊₁𝒜 = 𝐠𝐪ₙ₊₁[𝒜_indice,:]
    𝛎ₙ₊₁𝒜 = pad_λ(c,𝒜_indice,𝛎ₙ₊₁[𝒜_indice])
    Gᵖ = ∂𝐌𝐚∂𝐪(qₙ₊₁,Uₙ₊₁) - ∂𝒈𝒒T𝛌∂𝒒(qₙ₊₁,𝛎ₙ₊₁𝒜)
    Sₜᵖ = Array(BlockDiagonal(
        [
            [
                Mₙ₊₁+Gᵖ transpose(𝐠𝐪ₙ₊₁𝒜);
                𝐠𝐪ₙ₊₁𝒜 0I
            ],
            Matrix(1.0I,ā,ā)
        ]
    ))
end

function compute_Sₜᵛ(n,c,ū,xe,t,p,h,active_sets,dyfuncs)
    qₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁ = xe
    @unpack αm,αf,γ,β = p
    @unpack 𝒞,𝒰,𝒰c,𝒜,𝒜c,ℬ,ℬc = active_sets
    𝐌,𝒈,𝒈𝒒,𝐟,Jacobians = dyfuncs
    ∂𝒈𝒒T𝛌∂𝒒,∂𝒈𝒒𝐯∂𝒒,∂𝐌𝐚∂𝐪,∂𝐟∂𝐪,∂𝐟∂𝐯 = Jacobians
    ℬ_indice = sort(collect(ℬ))
    ℬc_indice = sort(collect(ℬc))
    b̄ = length(ℬc_indice)
    Mₙ₊₁ = 𝐌(qₙ₊₁)
    𝐠𝐪ₙ₊₁ = 𝒈𝒒(qₙ₊₁)
    𝐠𝐪ₙ₊₁ℬ = 𝐠𝐪ₙ₊₁[ℬ_indice,:]

    Sₜᵛ = Array(BlockDiagonal(
        [
            [
                Mₙ₊₁ transpose(𝐠𝐪ₙ₊₁ℬ);
                𝐠𝐪ₙ₊₁ℬ 0I
            ],
            Matrix(1.0I,b̄,b̄)
        ]
    ))

end

function initialize_St(n,c,ū,q0,v0,t,p,h,𝒞,𝒰c,𝐞,dyfuncs;tol=1e-14,imax=100)
    𝐌,𝒈,𝒈𝒒,𝐟,Jacobians = dyfuncs
    # ∂𝒈𝒒T𝛌∂𝒒,∂𝒈𝒒𝐯∂𝒒,∂𝐌𝐚∂𝐪,∂𝐟∂𝐪,∂𝐟∂𝐯 = Jacobians
    @unpack αm,αf,γ,β = p
    r = 1
    qₙ = q0
    vₙ = v0
    Mₙ = 𝐌(qₙ)
    fₙ = 𝐟(qₙ,vₙ,t)
    ṽ̇ₙ = Mₙ\fₙ
    aₙ = ṽ̇ₙ
    # @show ṽ̇0
    qₙ₊₁ ,vₙ₊₁ ,aₙ₊₁ ,ṽₙ₊₁ ,ṽ̇ₙ₊₁ ,Uₙ₊₁ ,Wₙ₊₁ ,𝛌ₙ₊₁ ,𝛎ₙ₊₁ ,𝚲ₙ₊₁  = initial_guesses(c,ū,qₙ,vₙ,aₙ,ṽ̇ₙ,p,h)
    xe = (qₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁)
    q̃ₙ₊₁ = copy(qₙ₊₁)
    active_sets = initialize_active_sets(𝒞,𝒰c)
    update_active_sets!(active_sets,qₙ,vₙ,q̃ₙ₊₁,xe,𝐞,𝒈,𝒈𝒒,r)
    for i = 1:imax
        xe = (qₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁)
        𝐫ˢ, 𝐫ᵖ, 𝐫ᵛ = get_residuals(qₙ,vₙ,xe,t,h,active_sets,𝐞,𝐌,𝒈,𝒈𝒒,𝐟)
        if maximum(norm.([𝐫ˢ, 𝐫ᵖ, 𝐫ᵛ])) < tol
            @show i
            break
        end
        # Step 1
        Sₜˢ = compute_Sₜˢ(n,c,ū,xe,t,p,h,active_sets,dyfuncs)
        Δxˢ = -inv(Sₜˢ)*𝐫ˢ
        Δṽ,Δ𝛌 = split_by_lengths(Δxˢ,[n,ū])
        Δ𝛌 /= -h
        ṽₙ₊₁ += Δṽ
        ṽ̇ₙ₊₁ += (1-αm)/(1-αf)/(γ*h)*Δṽ
        vₙ₊₁ = ṽₙ₊₁ + Wₙ₊₁
        qₙ₊₁ += h*β/γ*Δṽ
        𝛌ₙ₊₁ += Δ𝛌
        # step 2
        xe = (qₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁)
        _, 𝐫ᵖ, _ = get_residuals(qₙ,vₙ,xe,t,h,active_sets,𝐞,𝐌,𝒈,𝒈𝒒,𝐟)
        Sₜᵖ = compute_Sₜᵖ(n,c,ū,xe,t,p,h,active_sets,dyfuncs)
        Δxᵖ = -inv(Sₜᵖ)*𝐫ᵖ
        ΔU,Δ𝛎 = split_by_lengths(Δxᵖ,[n,c])
        ΔU *=  h
        Δ𝛎 *= -h
        Uₙ₊₁ += ΔU
        qₙ₊₁ += ΔU
        𝛎ₙ₊₁ += Δ𝛎
        # step 3
        xe = (qₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁)
        _, _, 𝐫ᵛ = get_residuals(qₙ,vₙ,xe,t,h,active_sets,𝐞,𝐌,𝒈,𝒈𝒒,𝐟)
        Sₜᵛ = compute_Sₜᵛ(n,c,ū,xe,t,p,h,active_sets,dyfuncs)
        Δxᵛ = -inv(Sₜᵛ)*𝐫ᵛ
        ΔW,Δ𝚲 = split_by_lengths(Δxᵛ,[n,c])
        Δ𝚲 *= -1
        Wₙ₊₁ += ΔW
        vₙ₊₁ = ṽₙ₊₁ + Wₙ₊₁
        𝚲ₙ₊₁ += Δ𝚲
        # if i == imax
        #     @error "Max iteration"
        # end
    end
    aₙ₊₁ += (1-αf)/(1-αm)*ṽ̇ₙ₊₁
end
