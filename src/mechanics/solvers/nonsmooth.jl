
mutable struct ActiveSets{T}
    𝒞::T
    𝒰::T
    𝒰c::T
    𝒜::T
    𝒜c::T
    ℬ::T
    ℬc::T
end

function initialize_active_sets(𝒞,𝒰c)
    𝒰 = setdiff(𝒞,𝒰c)
    𝒜 = copy(𝒰c)
    𝒜c = copy(𝒰)
    ℬ = copy(𝒰c)
    ℬc = copy(𝒰)
    ActiveSets(𝒞,𝒰,𝒰c,𝒜,𝒜c,ℬ,ℬc)
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

module NSGA
import ..Rible as TR
using LinearAlgebra
using BlockDiagonals
using Parameters

function find_corrections(ΔW,Δṽ,ΔU,p,h)
    @unpack αm,αf,γ,β = p
    Δv = Δṽ + ΔW
    Δṽ̇ = (1-αm)/((1-αf)*γ*h)*Δṽ
    Δq = h*β/γ*Δṽ + ΔU
    Δv,Δṽ,Δq
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
    𝐫ᵖ = vcat(
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

    Gᵖ = ∂𝐌𝐚∂𝐪(qₙ₊₁,Uₙ₊₁) - ∂𝒈𝒒T𝛌∂𝒒(qₙ₊₁,𝛎ₙ₊₁)
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

function initialize_St(n,c,ū,q0,v0,t,p,h,𝒞,𝒰c,𝐞,dyfuncs,tspan;tol=1e-14,imax=100)
    𝐌,𝒈,𝒈𝒒,𝐟,Jacobians = dyfuncs
    # ∂𝒈𝒒T𝛌∂𝒒,∂𝒈𝒒𝐯∂𝒒,∂𝐌𝐚∂𝐪,∂𝐟∂𝐪,∂𝐟∂𝐯 = Jacobians
    @unpack αm,αf,γ,β = p
    r = 1
    M0 = 𝐌(q0)
    f0 = 𝐟(q0,v0,t)
    ṽ̇0 = M0\f0
    a0 = ṽ̇0
    # @show ṽ̇0
    tstart,tend = tspan
    totaltime = tend - tstart
    totalstep = Int(ceil(totaltime/h))
    qs = [copy(q0) for i = 1:totalstep+1]
    vs = [copy(v0) for i = 1:totalstep+1]
    as = [copy(a0) for i = 1:totalstep+1]
    ṽ̇s = [copy(ṽ̇0) for i = 1:totalstep+1]
    for timestep = 1:totalstep
        qₙ = qs[timestep]
        vₙ = vs[timestep]
        aₙ = as[timestep]
        ṽ̇ₙ = ṽ̇s[timestep]
        qₙ₊₁ ,vₙ₊₁ ,aₙ₊₁ ,ṽₙ₊₁ ,ṽ̇ₙ₊₁ ,Uₙ₊₁ ,Wₙ₊₁ ,𝛌ₙ₊₁ ,𝛎ₙ₊₁ ,𝚲ₙ₊₁  = RB.initial_guesses(c,ū,qₙ,vₙ,aₙ,ṽ̇ₙ,p,h)
        xe = (qₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁)
        q̃ₙ₊₁ = copy(qₙ₊₁)
        active_sets = RB.initialize_active_sets(𝒞,𝒰c)
        update_active_sets!(active_sets,qₙ,vₙ,q̃ₙ₊₁,xe,𝐞,𝒈,𝒈𝒒,r)
        for i = 1:imax
            xe = (qₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁)
            residuals = get_residuals(qₙ,vₙ,xe,t,h,active_sets,𝐞,𝐌,𝒈,𝒈𝒒,𝐟)
            max_err,isub = findmax(norm.(residuals))
            if  max_err < tol
                # @show i
                break
            elseif i == imax
                @error "Reach max iteration $i, err=$max_err for the $isub subproblem"
                # @show residuals[isub]
                @show abs.(residuals[isub]) .> tol
            end
            # Step 1
            𝐫ˢ, _, _ = residuals
            Sₜˢ = compute_Sₜˢ(n,c,ū,xe,t,p,h,active_sets,dyfuncs)
            Δxˢ = -Sₜˢ\𝐫ˢ
            Δṽ,Δ𝛌 = RB.split_by_lengths(Δxˢ,[n,ū])
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
            Δxᵖ = -Sₜᵖ\𝐫ᵖ
            ΔU,Δ𝛎 = RB.split_by_lengths(Δxᵖ,[n,c])
            # ΔU *=  h
            Δ𝛎 *= -1
            Uₙ₊₁ += ΔU
            qₙ₊₁ += ΔU
            𝛎ₙ₊₁ += Δ𝛎
            # step 3
            xe = (qₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁)
            _, _, 𝐫ᵛ = get_residuals(qₙ,vₙ,xe,t,h,active_sets,𝐞,𝐌,𝒈,𝒈𝒒,𝐟)
            Sₜᵛ = compute_Sₜᵛ(n,c,ū,xe,t,p,h,active_sets,dyfuncs)
            Δxᵛ = -Sₜᵛ\𝐫ᵛ
            ΔW,Δ𝚲 = RB.split_by_lengths(Δxᵛ,[n,c])
            Δ𝚲 *= -1
            Wₙ₊₁ += ΔW
            vₙ₊₁ = ṽₙ₊₁ + Wₙ₊₁
            𝚲ₙ₊₁ += Δ𝚲
        end
        aₙ₊₁ += (1-αf)/(1-αm)*ṽ̇ₙ₊₁
        qs[timestep+1] .= qₙ₊₁
        vs[timestep+1] .= vₙ₊₁
        as[timestep+1] .= aₙ₊₁
        ṽ̇s[timestep+1] .= ṽ̇ₙ₊₁
    end
    qs
end

end

module RobustNSGA
import ..Rible as TR
using BlockDiagonals
using Parameters

@inline function 𝐟ᵖ(qₙ₊₁,q̃ₙ₊₁,ṽₙ₊₁,𝛌ₙ₊₁,t,𝒈𝒒,𝐟)
    𝐟(qₙ₊₁,ṽₙ₊₁,t) - 𝐟(q̃ₙ₊₁,ṽₙ₊₁,t) +
    transpose(𝒈𝒒(qₙ₊₁)-𝒈𝒒(q̃ₙ₊₁))*𝛌ₙ₊₁
end

@inline function 𝐟ˣ(qₙ₊₁,vₙ₊₁,q̃ₙ₊₁,ṽₙ₊₁,ṽ̇ₙ₊₁,𝛌ₙ₊₁,t,𝐌,𝒈𝒒,𝐟)
    𝐟(qₙ₊₁,vₙ₊₁,t) - 𝐟(q̃ₙ₊₁,ṽₙ₊₁,t) +
    transpose(𝒈𝒒(qₙ₊₁)-𝒈𝒒(q̃ₙ₊₁))*𝛌ₙ₊₁ +
    (𝐌(q̃ₙ₊₁)-𝐌(qₙ₊₁))*ṽ̇ₙ₊₁
end

function get_𝐫ˢ(qₙ,vₙ,xe,t,h,kₛ,pₛ,active_sets,𝐞,𝐌,𝒈,𝒈𝒒,𝐟)
    @unpack 𝒰,𝒰c,𝒜,𝒜c,ℬ,ℬc = active_sets
    qₙ₊₁, q̃ₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁ = xe
    M̃ₙ₊₁ = 𝐌(q̃ₙ₊₁)
    f̃ₙ₊₁ = 𝐟(q̃ₙ₊₁,ṽₙ₊₁,t)
    𝒰c_indice = sort(collect(𝒰c))
    𝐠𝐪ₙ₊₁𝒰c = @view 𝐠𝐪ₙ₊₁[𝒰c_indice,:]
    𝐫ˢ = vcat(
        M̃ₙ₊₁*ṽ̇ₙ₊₁ - f̃ₙ₊₁ - transpose(𝐠𝐪ₙ₊₁𝒰c)*(kₛ*𝛌ₙ₊₁-pₛ*𝐠𝐪ₙ₊₁𝒰c*ṽₙ₊₁),
        -kₛ*𝐠𝐪ₙ₊₁𝒰c*ṽₙ₊₁
    )
end

function get_𝐫ᵖ(qₙ,vₙ,xe,t,h,kₚ,pₚ,active_sets,𝐞,𝐌,𝒈,𝒈𝒒,𝐟)
    @unpack 𝒰,𝒰c,𝒜,𝒜c,ℬ,ℬc = active_sets
    qₙ₊₁, q̃ₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁ = xe
    M̃ₙ₊₁ = 𝐌(q̃ₙ₊₁)
    fᵖₙ₊₁ = 𝐟ᵖ(qₙ₊₁,q̃ₙ₊₁,ṽₙ₊₁,𝛌ₙ₊₁,t,𝒈𝒒,𝐟)
    𝛏ₙ₊₁ = kₚ*𝛎ₙ₊₁ - pₚ*𝐠ₙ₊₁
    𝒜 = 𝒰c ∪ Set([j for j ∈ 𝒰 if 𝛏ₙ₊₁[j] ≥ 0])
    𝒜c = setdiff(𝒞,𝒜)
    @pack! active_sets = 𝒜,𝒜c
    𝒜_indice = sort(collect(𝒜))
    𝒜c_indice = sort(collect(𝒜c))
    𝐠𝐪ₙ₊₁𝒜 = @view 𝐠𝐪ₙ₊₁[𝒜_indice,:]
    𝐠ₙ₊₁𝒜 = @view 𝐠ₙ₊₁[𝒜_indice]
    𝐫ᵖ = vcat(
        M̃ₙ₊₁*Uₙ₊₁ - h^2*fᵖₙ₊₁ - transpose(𝐠𝐪ₙ₊₁𝒜)*𝛏ₙ₊₁[𝒜_indice],
        -kₚ*𝐠ₙ₊₁𝒜,
        -kₚ^2/pₚ*𝛎ₙ₊₁[𝒜c_indice]
    )
end

function get_𝐫ᵛ(qₙ,vₙ,xe,t,h,kᵥ,pᵥ,active_sets,𝐞,𝐌,𝒈,𝒈𝒒,𝐟)
    @unpack 𝒰,𝒰c,𝒜,𝒜c,ℬ,ℬc = active_sets
    qₙ₊₁, q̃ₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁ = xe
    Mₙ₊₁ = 𝐌(qₙ₊₁)
    fˣₙ₊₁ = 𝐟ˣ(qₙ₊₁,vₙ₊₁,q̃ₙ₊₁,ṽₙ₊₁,ṽ̇ₙ₊₁,𝛌ₙ₊₁,t,𝐌,𝒈𝒒,𝐟)
    g̊ₙ₊₁ = 𝐠𝐪ₙ₊₁*vₙ₊₁ + Diagonal(𝐞)*𝐠𝐪ₙ*vₙ
    𝛔ₙ₊₁ = kᵥ*𝚲ₙ₊₁ - pᵥ*g̊ₙ₊₁
    ℬ = 𝒰c ∪ Set([j for j ∈ 𝒜 if 𝛔ₙ₊₁[j] ≥ 0])
    ℬc = setdiff(𝒞,ℬ)
    @pack! active_sets = ℬ,ℬc
    ℬ_indice = sort(collect(ℬ))
    ℬc_indice = sort(collect(ℬc))
    𝐠𝐪ₙ₊₁ℬ = @view 𝐠𝐪ₙ₊₁[ℬ_indice,:]
    𝐠ₙ₊₁ℬ = @view 𝐠ₙ₊₁[ℬ_indice]
    𝐫ᵛ = vcat(
        Mₙ₊₁*Wₙ₊₁ - h*fˣₙ₊₁ - transpose(𝐠𝐪ₙ₊₁ℬ)*𝛔ₙ₊₁[ℬ_indice],
        -kᵥ*g̊ₙ₊₁[ℬ_indice],
        -kᵥ^2/pᵥ*𝚲ₙ₊₁[ℬc_indice]
    )
end

function compute_Sₜˢ(n,c,ū,xe,t,p,h,kₛ,pₛ,active_sets,dyfuncs)
    qₙ₊₁, q̃ₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁ = xe
    @unpack αm,αf,γ,β = p
    @unpack 𝒞,𝒰,𝒰c,𝒜,𝒜c,ℬ,ℬc = active_sets
    𝐌,𝒈,𝒈𝒒,𝐟,Jacobians = dyfuncs
    ∂𝒈𝒒T𝛌∂𝒒,∂𝒈𝒒𝐯∂𝒒,∂𝐌𝐚∂𝐪,∂𝐟∂𝐪,∂𝐟∂𝐯 = Jacobians

    𝒰c_indice = sort(collect(𝒰c))

    M̃ₙ₊₁ = 𝐌(q̃ₙ₊₁)

    𝛌 = pad_λ(c,𝒰c_indice,𝛌ₙ₊₁)

    𝐠𝐪ₙ₊₁ = 𝒈𝒒(qₙ₊₁)
    𝐠𝐪ₙ₊₁𝒰c = 𝐠𝐪ₙ₊₁[𝒰c_indice,:]
    ∂𝐟∂𝐪̃ₙ₊₁ = ∂𝐟∂𝐪(q̃ₙ₊₁,ṽₙ₊₁,t)
    Kₜ = ∂𝐌𝐚∂𝐪(q̃ₙ₊₁,ṽ̇ₙ₊₁) - ∂𝒈𝒒T𝛌∂𝒒(q̃ₙ₊₁,kₛ*𝛌-pₛ*𝐠𝐪ₙ₊₁𝒰c*ṽₙ₊₁) - ∂𝐟∂𝐪̃ₙ₊₁
    Cₜ = pₛ*transpose(𝐠𝐪ₙ₊₁𝒰c)*𝐠𝐪ₙ₊₁𝒰c  - ∂𝐟∂𝐪̃ₙ₊₁
    Sₜˣ = (1-αm)/h/(1-αf)/γ*Mₙ₊₁ + Cₜ + h*β/γ*Kₜ

    Gˣˢ = ∂𝒈𝒒𝐯∂𝒒(q̃ₙ₊₁,ṽₙ₊₁)[𝒰c_indice,:]

    Sₜˢ = [
        Sₜˣ -kₛ*transpose(𝐠𝐪ₙ₊₁𝒰c);
        -kₛ*(𝐠𝐪ₙ₊₁𝒰c+h*β/γ*Gˣˢ) 0I
    ]
end

function compute_Sₜᵖ(n,c,ū,xe,t,p,h,kₚ,pₚ,active_sets,dyfuncs)
    qₙ₊₁, q̃ₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁ = xe
    @unpack αm,αf,γ,β = p
    @unpack 𝒞,𝒰,𝒰c,𝒜,𝒜c,ℬ,ℬc = active_sets
    𝐌,𝒈,𝒈𝒒,𝐟,Jacobians = dyfuncs
    ∂𝒈𝒒T𝛌∂𝒒,∂𝒈𝒒𝐯∂𝒒,∂𝐌𝐚∂𝐪,∂𝐟∂𝐪,∂𝐟∂𝐯 = Jacobians
    𝒜_indice = sort(collect(𝒜))
    𝒜c_indice = sort(collect(𝒜c))
    ā = length(𝒜c_indice)
    M̃ₙ₊₁ = 𝐌(q̃ₙ₊₁)
    𝐠𝐪ₙ₊₁ = 𝒈𝒒(qₙ₊₁)
    𝐠𝐪ₙ₊₁𝒜 = 𝐠𝐪ₙ₊₁[𝒜_indice,:]

    Sₜᵖˣ = M̃ₙ₊₁ - ∂𝒈𝒒T𝛌∂𝒒(qₙ₊₁,𝛏ₙ₊₁) - h^2*∂𝐟ᵖ∂𝐪(qₙ₊₁,q̃ₙ₊₁,ṽₙ₊₁,𝛌ₙ₊₁,t)

    Sₜᵖ = Array(BlockDiagonal(
        [
            [
                Sₜᵖˣ -kₚ*transpose(𝐠𝐪ₙ₊₁𝒜);
                -kₚ*𝐠𝐪ₙ₊₁𝒜 0I
            ],
            Matrix(-kₚ^2/pₚ*I,ā,ā)
        ]
    ))
end

function compute_Sₜᵛ(n,c,ū,xe,t,p,h,kᵥ,pᵥ,active_sets,dyfuncs)
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
    𝐠𝐪ₙ₊₁ℬc = 𝐠𝐪ₙ₊₁[ℬc_indice,:]

    Sₜᵛˣ = 𝐌(qₙ₊₁) + pᵥ*transpose(𝐠𝐪ₙ₊₁ℬc)*𝐠𝐪ₙ₊₁ℬc - h*∂𝐟∂𝐯(qₙ₊₁,vₙ₊₁,t)

    Sₜᵛ = Array(BlockDiagonal(
        [
            [
                Sₜᵛˣ -kᵥ*transpose(𝐠𝐪ₙ₊₁ℬ);
                -kᵥ*𝐠𝐪ₙ₊₁ℬ 0I
            ],
            -kᵥ^2/pᵥ*Matrix(1.0I,b̄,b̄)
        ]
    ))
end

function robustnsga(n,c,ū,q0,v0,t,p,h,𝒞,𝒰c,𝐞,dyfuncs,tspan;tol=1e-14,imax=100)
    𝐌,𝒈,𝒈𝒒,𝐟,Jacobians = dyfuncs
    # ∂𝒈𝒒T𝛌∂𝒒,∂𝒈𝒒𝐯∂𝒒,∂𝐌𝐚∂𝐪,∂𝐟∂𝐪,∂𝐟∂𝐯 = Jacobians
    @unpack αm,αf,γ,β = p
    M0 = 𝐌(q0)
    m̄ = maximum(M0)
    kₛ = pₛ = m̄/h; kₚ = pₚ = m̄; kᵥ = pᵥ = m̄
    f0 = 𝐟(q0,v0,t)
    ṽ̇0 = M0\f0
    a0 = ṽ̇0
    # @show ṽ̇0
    tstart,tend = tspan
    totaltime = tend - tstart
    totalstep = Int(ceil(totaltime/h))
    qs = [copy(q0) for i = 1:totalstep+1]
    vs = [copy(v0) for i = 1:totalstep+1]
    as = [copy(a0) for i = 1:totalstep+1]
    ṽ̇s = [copy(ṽ̇0) for i = 1:totalstep+1]
    for timestep = 1:totalstep
        qₙ = qs[timestep]
        vₙ = vs[timestep]
        aₙ = as[timestep]
        ṽ̇ₙ = ṽ̇s[timestep]
        qₙ₊₁ ,vₙ₊₁ ,aₙ₊₁ ,ṽₙ₊₁ ,ṽ̇ₙ₊₁ ,Uₙ₊₁ ,Wₙ₊₁ ,𝛌ₙ₊₁ ,𝛎ₙ₊₁ ,𝚲ₙ₊₁  = initial_guesses(c,ū,qₙ,vₙ,aₙ,ṽ̇ₙ,p,h)
        q̃ₙ₊₁ = copy(qₙ₊₁)
        xe = (qₙ₊₁, q̃ₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁)

        active_sets = RB.initialize_active_sets(𝒞,𝒰c)
        update_active_sets!(active_sets,qₙ,vₙ,q̃ₙ₊₁,xe,𝐞,𝒈,𝒈𝒒,r)
        # Step 1
        for i = 1:imax
            xe = (qₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁)
            𝐫ˢ = get_𝐫ˢ(qₙ,vₙ,xe,t,h,active_sets,𝐞,𝐌,𝒈,𝒈𝒒,𝐟)
            norm_𝐫ˢ = norm(𝐫ˢ)
            if  norm(𝐫ˢ) < tol
                break
            elseif i == imax
                @error "Reach max iteration $i, err=$norm_𝐫ˢ for the first subproblem"
                @show abs.(𝐫ˢ) .> tol
            end
            Sₜˢ = compute_Sₜˢ(n,c,ū,xe,t,p,h,active_sets,dyfuncs)
            Δxˢ = -Sₜˢ\𝐫ˢ
            Δṽ,Δ𝛌 = RB.split_by_lengths(Δxˢ,[n,ū])
            ṽₙ₊₁ += Δṽ
            ṽ̇ₙ₊₁ += (1-αm)/(1-αf)/(γ*h)*Δṽ
            qₙ₊₁ += h*β/γ*Δṽ
            𝛌ₙ₊₁ += Δ𝛌
        end
        # step 2
        for i = 1:imax
            xe = (qₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁)
            𝐫ᵖ = get_𝐫ᵖ(qₙ,vₙ,xe,t,h,active_sets,𝐞,𝐌,𝒈,𝒈𝒒,𝐟)
            norm_𝐫ᵖ = norm(𝐫ᵖ)
            if  norm(𝐫ᵖ) < tol
                break
            elseif i == imax
                @error "Reach max iteration $i, err=$norm_𝐫ᵖ for the first subproblem"
                @show abs.(𝐫ᵖ) .> tol
            end
            Sₜᵖ = compute_Sₜᵖ(n,c,ū,xe,t,p,h,active_sets,dyfuncs)
            Δxᵖ = -Sₜᵖ\𝐫ᵖ
            ΔU,Δ𝛎 = RB.split_by_lengths(Δxᵖ,[n,c])
            Uₙ₊₁ += ΔU
            qₙ₊₁ += ΔU
            𝛎ₙ₊₁ += Δ𝛎
        end
        # step 3
        for i = 1:imax
            xe = (qₙ₊₁, vₙ₊₁, ṽ̇ₙ₊₁, ṽₙ₊₁, 𝛌ₙ₊₁, Uₙ₊₁, 𝛎ₙ₊₁, Wₙ₊₁, 𝚲ₙ₊₁)
            𝐫ᵛ = get_𝐫ᵛ(qₙ,vₙ,xe,t,h,active_sets,𝐞,𝐌,𝒈,𝒈𝒒,𝐟)
            norm_𝐫ᵛ = norm(𝐫ᵛ)
            if  norm(𝐫ᵛ) < tol
                break
            elseif i == imax
                @error "Reach max iteration $i, err=$norm_𝐫ᵛ for the first subproblem"
                @show abs.(𝐫ᵛ) .> tol
            end
            Sₜᵛ = compute_Sₜᵛ(n,c,ū,xe,t,p,h,active_sets,dyfuncs)
            Δxᵛ = -Sₜᵛ\𝐫ᵛ
            ΔW,Δ𝚲 = RB.split_by_lengths(Δxᵛ,[n,c])
            Wₙ₊₁ += ΔW
            vₙ₊₁ = ṽₙ₊₁ + Wₙ₊₁
            𝚲ₙ₊₁ += Δ𝚲
        end
        aₙ₊₁ += (1-αf)/(1-αm)*ṽ̇ₙ₊₁
        qs[timestep+1] .= qₙ₊₁
        vs[timestep+1] .= vₙ₊₁
        as[timestep+1] .= aₙ₊₁
        ṽ̇s[timestep+1] .= ṽ̇ₙ₊₁
    end
    qs
end

end #module

include("nssfc.jl")
