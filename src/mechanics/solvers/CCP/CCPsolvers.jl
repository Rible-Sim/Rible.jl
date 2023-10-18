
@inline next_θ(θ) = θ*(-θ+sqrt(θ^2+4))/2
@inline next_β(θₖ,θₖ₊₁) = (1-θₖ)*θₖ/(θₖ^2+θₖ₊₁)

function ∏c(μ::Real,γ::AbstractVector)
    γₙ = γ[1]
    γᵤ = γ[2]
    γᵥ = γ[3]
    o = zero(μ)
    γₜ = sqrt(γᵤ^2+γᵥ^2)
    if γₜ < μ*γₙ # in the cone
        return SVector(γₙ,γᵤ,γᵥ)
    elseif γₜ < -γₙ/μ # in the polar cone
        return SVector(o,o,o)
    else
        ∏ₙ = (γₜ*μ + γₙ)/(μ^2 + 1)
        ∏ᵤ = γᵤ*(μ*∏ₙ/γₜ)
        ∏ᵥ = γᵥ*(μ*∏ₙ/γₜ)
        return SVector(∏ₙ,∏ᵤ,∏ᵥ)
    end
end

function ∏c(μs::AbstractVector,γs::AbstractVector)
    ret = similar(γs)
    is = 0
    for i in 1:length(μs)
        is = 3(i-1)
        ret[is+1:is+3] = ∏c(μs[i],γs[is+1:is+3])
    end
    ret
end

function make_residual1(μs,𝐍,𝐫)
    nc = length(μs)
    function r(𝛄)
        𝐱 = 𝐍*𝛄 + 𝐫
        𝒇 = [min(𝐱[3(i-1)+1],0) for i in 1:nc]
        norm(𝒇,Inf)
    end
    r
end

function make_residual4(μs,𝐍,𝐫;gd=1e-6)
    nr = length(𝐫)
    𝛙(𝛄) = 1/(nr*gd)*(𝛄 - ∏c(μs,𝛄-gd*(𝐍*𝛄 + 𝐫)))
    r(𝛄) = norm(𝛙(𝛄))
    r
end

function make_B(nc,D,invM)
    g = zeros(eltype(invM),nc)
    for i = 1:nc
        Di = @view D[3(i-1)+1:3(i-1)+3,:]
        g[i] = 1#tr(Di*invM*transpose(Di))/3
    end
    BlockDiagonal([Matrix(gi*I,3,3) for gi in g])
end

function Jacobi(B,r,μs,𝐍,𝐫,𝛄0=zero(𝐫);τ=1e-5,Nmax=20,ω=0.3,λ=1.1)
    𝛄ₖ₊₁ = 𝛄0
    𝛄̂ₖ₊₁ = copy(𝛄ₖ₊₁)
    j = 0
    res = zero(eltype(𝛄ₖ₊₁))
    for k = 1:Nmax
        𝛄̂ₖ₊₁ .= ∏c(μs,𝛄ₖ₊₁ - ω*B*(𝐍*𝛄ₖ₊₁ + 𝐫))
        𝛄ₖ₊₁ .= λ.*𝛄̂ₖ₊₁ .+ (1-λ).*𝛄ₖ₊₁
        res = r(𝛄ₖ₊₁)
        # @info "k=$k, res=$res"
        # @show 𝛄ₖ₊₁
        if res < τ
            j = k
            break
        end
        if k == Nmax
            @error "Jacobi: Max iteration $k, res=$res"
        end
    end
    𝛄ₖ₊₁,j,res
end

function GaussSeidel(nc,B,r,μs,𝐍,𝐫,𝛄0=zero(𝐫);τ=1e-5,Nmax=20,ω=1.0,λ=1.0)
    𝛄ₖ = 𝛄0
    𝛄ₖ₊₁ = copy(𝛄ₖ)
    j = 0
    res = zero(eltype(𝛄ₖ₊₁))
    gd = 1e-6
    nr = length(𝐫)
    𝛙(𝛄) = 1/(nr*gd)*(𝛄 - ∏c(μs,𝛄-gd*(𝐍*𝛄 + 𝐫)))
    for k = 1:Nmax
        for i = 1:nc
            is = 3(i-1)
            # Bᵢ = @view B[is+1:is+3,is+1:is+3]
            𝛄̂ᵢₖ₊₁ = ∏c(μs[i],𝛄ₖ[is+1:is+3] - ω*((𝐍[is+1:is+3,:]*𝛄ₖ + 𝐫[is+1:is+3])))
            𝛄ₖ₊₁[is+1:is+3] .= λ.*𝛄̂ᵢₖ₊₁ .+ (1-λ).*𝛄ₖ[is+1:is+3]
            𝛄ₖ[is+1:is+3] = 𝛄ₖ₊₁[is+1:is+3]
        end
        𝛄ₖ .= 𝛄ₖ₊₁
        @show 𝛙(𝛄ₖ₊₁)
        res = r(𝛄ₖ₊₁)
        j = k
        if res < τ
            break
        end
        if k == Nmax
            @error "GaussSeidel: Max iteration, res=$res"
        end
    end
    𝛄ₖ₊₁,j,res
end

function APGD!(output,r,μs,𝐍,𝐫;τ=1e-5,Nmax=20)
    f(𝛄) = 1/2*transpose(𝛄)*𝐍*𝛄 + transpose(𝛄)*𝐫
    ∇f(𝛄) = 𝐍*𝛄 + 𝐫
    nr = length(𝐫)
    θₖ = θ0 = 1
    𝛄0 = zero(𝐫)
    𝛄ₖ = copy(𝛄0)
    𝐲ₖ = copy(𝛄0)
    𝛄ₕ = ones(eltype(𝛄0),length(𝛄0))
    𝐲ₖ₊₁ = copy(𝛄0)
    𝛄ₖ₊₁ = copy(𝛄0)
    Lₖ = norm(∇f(𝛄0)-∇f(𝛄ₕ))/norm(𝛄0-𝛄ₕ)
    tₖ = 1/Lₖ
    rmin = 1e14
    for k = 1:Nmax
        # solve
        gₖ = ∇f(𝐲ₖ)
        𝛄ₖ₊₁ .= ∏c(μs,𝐲ₖ-tₖ*gₖ)
        while f(𝛄ₖ₊₁) > f(𝐲ₖ) + transpose(gₖ)*(𝛄ₖ₊₁-𝐲ₖ) + 1/2*Lₖ*sum((𝛄ₖ₊₁.-𝐲ₖ).^2)
            Lₖ *= 2
            tₖ = 1/Lₖ
            𝛄ₖ₊₁ .= ∏c(μs,𝐲ₖ-tₖ*gₖ)
        end
        # update
        θₖ₊₁ = next_θ(θₖ)
        βₖ₊₁ = next_β(θₖ,θₖ₊₁)
        𝐲ₖ₊₁ .= 𝛄ₖ₊₁ .+ βₖ₊₁.*(𝛄ₖ₊₁ .- 𝛄ₖ)
        # exit
        rₖ₊₁ = r(𝛄ₖ₊₁)
        # @show k,rₖ₊₁
        if rₖ₊₁ < rmin
            rmin = rₖ₊₁
            𝛄ₕ .= 𝛄ₖ₊₁
        end
        if rₖ₊₁ < τ
            break
        end
        if transpose(gₖ)*(𝛄ₖ₊₁ - 𝛄ₖ) > 0
            𝐲ₖ₊₁ .= 𝛄ₖ₊₁
            θₖ₊₁ = 1
        end
        if k == Nmax
            # @error "APGD: Max iteration $k, res=$(rₖ₊₁)"
        end
        # update
        𝐲ₖ .= 𝐲ₖ₊₁
        𝛄ₖ .= 𝛄ₖ₊₁
        Lₖ = 0.9Lₖ
        tₖ = 1/Lₖ

        # μ = transpose(∇f(𝛄ₖ))*𝛄ₖ/nr
        # @show μ
    end
    output .= 𝛄ₕ
    return rmin
end

function IPM!(output,nu,nΛ,Λ,y,N,r;ftol=1e-14,Nmax=50)
    T = eltype(Λ)
    e = [one(T),zero(T),zero(T)]
    J = Diagonal([one(T),-one(T),-one(T)])
    𝐞_split = [e for i = 1:nu]
    Λ_split = split_by_lengths(Λ,3)
    y_split = split_by_lengths(y,3)
    # Λ_cone = [transpose(Λi)*J*Λi for Λi in Λ_split]
    # y_cone = [transpose(yi)*J*yi for yi in y_split]
    # @show Λ_cone
    # @show y_cone
    # W_blocks = NTScale.(Λ_split,y_split)
    # z_split = W_blocks.*Λ_split
    # z = reduce(vcat,z_split)
    # W = BlockDiagonal(W_blocks)
    # WᵀW = transpose(W)*W
    # @show z_split
    # @show z_split⊙z_split
    n1 = nΛ
    n2 = 2nΛ
    nx = n2
    𝐫𝐞𝐬 = zeros(T,nx)
    𝐉 = zeros(T,nx,nx)
    μ = transpose(y)*Λ/nΛ
    for k = 1:Nmax

        𝐫𝐞𝐬[   1:n1] .= N*Λ .+ r .- y
        𝐫𝐞𝐬[n1+1:n2] .= reduce(vcat,Λ_split⊙y_split)

        res = norm(𝐫𝐞𝐬)
        if res < ftol
            # @show k, Λ_split[1],y_split[1]
            # @show Λ_split[1]⊙y_split[1]
            # @show k
            break
        elseif k == Nmax
            # @warn "IPM: Max iteration $k reached"
        end

        𝐉[   1:n1,   1:n1] .=  N
        𝐉[   1:n1,n1+1:n2] .= -Matrix(1I,nΛ,nΛ)

        𝐉[n1+1:n2,   1:n1] .=  BlockDiagonal(mat.(y_split))
        𝐉[n1+1:n2,n1+1:n2] .=  BlockDiagonal(mat.(Λ_split))

        # @show Λₛ,yₛ
        # η = 1.0
        lu𝐉 = lu(𝐉)
        Δxp = lu𝐉\(-𝐫𝐞𝐬)
        ΔΛp = @view Δxp[   1:n1]
        Δyp = @view Δxp[n1+1:n2]
        ΔΛp_split = split_by_lengths(ΔΛp,3)
        Δyp_split = split_by_lengths(Δyp,3)
        # @show ΔΛp, Δyp
        # @show z_split,W_blocks,Δyp_split,ΔΛp_split,J
        αp_Λ = find_cone_step_length(Λ_split,ΔΛp_split,J)
        αp_y = find_cone_step_length(y_split,Δyp_split,J)
        αpmax = min(αp_Λ,αp_y)
        # αpmax = find_cone_step_length(z_split,W_blocks,Δyp_split,ΔΛp_split,J)
        αp = min(one(αpmax),0.99αpmax)
        # Λp_split = Λ_split .+ αp.*ΔΛp_split
        # yp_split = y_split .+ αp.*Δyp_split
        # Λp_cone = [transpose(Λi)*J*Λi for Λi in Λp_split]
        # yp_cone = [transpose(yi)*J*yi for yi in yp_split]
        # @show Λp_cone
        # @show yp_cone
        Λp = Λ .+ αp.*ΔΛp
        yp = y .+ αp.*Δyp
        μp = transpose(yp)*Λp/nΛ
        σ = (μp/μ)^3
        if σ == NaN || μ == 0
            break
        end
        τ = σ*μp
        # @show "Prediction",αpmax,αp,τ,σ,μ,μp
        # @show Δxp
        # @show αp.*ΔΛp,αp.*Δyp
        # @show Λp,yp
        𝐫𝐞𝐬_c_split = -τ.*𝐞_split.+((Δyp_split)⊙(ΔΛp_split))
        𝐫𝐞𝐬[n1+1:n2] .+= reduce(vcat,𝐫𝐞𝐬_c_split)
        # res = norm(𝐫𝐞𝐬)
        # @show 𝐫𝐞𝐬
        # @show res
        Δxc = lu𝐉\(-𝐫𝐞𝐬)
        # η = exp(-0.1μ) + 0.9
        ΔΛc = @view Δxc[   1:n1]
        Δyc = @view Δxc[n1+1:n2]
        ΔΛc_split = split_by_lengths(ΔΛc,3)
        Δyc_split = split_by_lengths(Δyc,3)
        # αmax = find_cone_step_length(z_split,W_blocks,Δyc_split,ΔΛc_split,J)
        α_Λ = find_cone_step_length(Λ_split,ΔΛc_split,J)
        # @show Λ_split,ΔΛc_split
        α_y = find_cone_step_length(y_split,Δyc_split,J)
        αmax = min(α_Λ,α_y)
        α = min(1,0.9αmax)
        Λ_split .+= α.*ΔΛc_split
        y_split .+= α.*Δyc_split
        # @show ΔΛc_split, Δyc_split
        # Λ_cone = [transpose(Λi)*J*Λi for Λi in Λ_split]
        # y_cone = [transpose(yi)*J*yi for yi in y_split]
        Λ .+= α.*ΔΛc
        y .+= α.*Δyc
        μ = transpose(y)*Λ/nΛ
        # @show μ
        # @show "after",y,Λ
    end
    # @show Λ_split
    # @show y_split
    output .= Λ
    y_split
end

function ⊙(x::AbstractVector{T},y::AbstractVector{T}) where {T<:Real}
    x0 = x[begin]
    x1 = @view x[begin+1:end]
    y0 = y[begin]
    y1 = @view y[begin+1:end]
    vcat(transpose(x)*y, x0.*y1 .+ y0.*x1)
end

function mat(x::AbstractVector{T}) where {T<:Real}
    x0 = x[begin]
    x1 = @view x[begin+1:end]
    nx = length(x)
    X = similar(x,nx,nx)
    X[begin,:] .= x
    X[:,begin] .= x
    X[begin+1:end,begin+1:end] .= Matrix(x0*I,nx-1,nx-1)
    X
end

function ⊘(x::AbstractVector{T},y::AbstractVector{T}) where {T<:Real}
    X = mat(x)
    inv(X)*y
end

function NTScale(x::AbstractVector{T},y::AbstractVector{T}) where {T<:Real}
    J = Diagonal(vcat(one(x[begin]),-one.(x[begin+1:end])))
    x̌ = transpose(x)*J*x
    y̌ = transpose(y)*J*y
    x̄ = x./sqrt(x̌)
    ȳ = y./sqrt(y̌)
    γ = sqrt((1+transpose(x̄)*ȳ)/2)
    w̄ = (ȳ + J*x̄)./(2γ)
    w̄0 = w̄[begin]
    w̄1 = @view w̄[begin+1:end]
    W̄ = zeros(eltype(w̄),length(w̄),length(w̄))
    W̄[1,:] .= w̄
    W̄[:,1] .= w̄
    W̄[2:end,2:end] = I+w̄1*transpose(w̄1)/(w̄0+1)
    θ = (y̌/x̌)^(1/4)
    W = θ*W̄
    # @show sqrt(transpose(x)*J*x),sqrt(transpose(y)*J*y)
    # @show W
end

function NTScale_Anderson(x::AbstractVector{T},y::AbstractVector{T}) where {T<:Real}
    Q = Diagonal(vcat(one(x[begin]),-one.(x[begin+1:end])))
    x̌ = transpose(x)*Q*x
    y̌ = transpose(y)*Q*y
    θ² = sqrt(y̌/x̌)
    θ = sqrt(θ²)
    w = (y./θ .+ θ.*(Q*x))/(
    (√2)*sqrt(transpose(x)*y+sqrt(x̌*y̌))
    )
    w0 = w[begin]
    w1 = @view w[begin+1:end]
    W = zeros(eltype(w),length(w),length(w))
    W[1,:] .= w
    W[:,1] .= w
    W[2:end,2:end] = I+w1*transpose(w1)/(w0+1)
    Θ = θ*I
    ΘW = Θ*W
end

function ⊙(x::AbstractVector{T},y::AbstractVector{T}) where {T<:AbstractVector}
    x.⊙y
end

function ⊘(x::AbstractVector{T},y::AbstractVector{T}) where {T<:AbstractVector}
    x.⊘y
end

function NTScale(x::AbstractVector{T},y::AbstractVector{T}) where {T<:AbstractVector}
     BlockDiagonal(NTScale.(x,y))
end

function find_cone_step_length(z::AbstractVector{T},Δz::AbstractVector{T},z̄::T,J) where {T<:Real}
    ẑ = z./z̄
    head = transpose(ẑ)*J*Δz
    tail = Δz[begin+1:end] .- (head+Δz[begin])./(ẑ[begin] + 1).*ẑ[begin+1:end]
    ρ = vcat(head,tail)/z̄
    α = max(0,norm(ρ[begin+1:end])-ρ[begin])^(-1)
end

function find_cone_step_length(z::AbstractVector{T},Δz::AbstractVector{T},J) where {T<:Real}
    z̄² = transpose(z)*J*z
    if z̄² ≤ 0
        α = zero(z̄²)
        # @warn "Zero Step Length: On Boundary"
    else
        z̄ = sqrt(z̄²)
        α = find_cone_step_length(z,Δz,z̄,J)
    end
    α
end

function find_cone_step_length(z_split::AbstractVector{T},Δz_split::AbstractVector{T},J) where {T<:AbstractVector}
    α = [
        find_cone_step_length(zi,Δzi,J)
        for (zi,Δzi) in zip(z_split,Δz_split)
    ]
    αmin = minimum(α)
end

# function find_cone_step_length(z_split,Δz_split,J)
#     z_split_norm = [sqrt(transpose(zi)*J*zi) for zi in z_split]
#     α = [
#         find_cone_step_length(zi,Δzi,zi_norm,J)
#         for (zi,Δzi,zi_norm) in zip(z_split,Δz_split,z_split_norm)
#     ]
#     αmin = minimum(α)
# end

function find_cone_step_length(z_split,W_blocks,Δy_split,Δx_split,J)
    Δx̃_split = W_blocks.*Δx_split
    Δỹ_split = inv.(W_blocks).*Δy_split
    z_split_norm = [sqrt(transpose(zi)*J*zi) for zi in z_split]
    # z̄_split = z_split./z_split_norm
    # ρ = [
    #     begin
    #         head = transpose(z̄i)*J*Δỹi
    #         tail = Δỹi[begin+1:end] .- (head+Δỹi[begin])./(z̄i[begin] + 1).*z̄i[begin+1:end]
    #         vcat(head,tail)/zi_norm
    #     end
    #     for (z̄i,Δỹi,zi_norm) in zip(z̄_split,Δỹ_split,z_split_norm)
    # ]
    #
    # σ = [
    #     begin
    #         head = transpose(z̄i)*J*Δx̃i
    #         tail = Δx̃i[begin+1:end] .- (head+Δx̃i[begin])./(z̄i[begin] + 1).*z̄i[begin+1:end]
    #         vcat(head,tail)/zi_norm
    #     end
    #     for (z̄i,Δx̃i,zi_norm) in zip(z̄_split,Δx̃_split,z_split_norm)
    # ]
    α_y = [
        find_cone_step_length(zi,Δỹi,zi_norm,J)
        for (zi,Δỹi,zi_norm) in zip(z_split,Δỹ_split,z_split_norm)
    ]
    α_x = [
        find_cone_step_length(zi,Δx̃i,zi_norm,J)
        for (zi,Δx̃i,zi_norm) in zip(z_split,Δx̃_split,z_split_norm)
    ]
    # @show [max(0,norm(ρi[begin+1:end])-ρi[begin])^(-1) for ρi in ρ]
    # @show [max(0,norm(σi[begin+1:end])-σi[begin])^(-1) for σi in σ]
    αmin = min(minimum(α_y),minimum(α_x))
    # αmin = minimum(max(0,norm(ρi[begin+1:end])-ρi[begin],norm(σi[begin+1:end])-σi[begin])^(-1) for (ρi,σi) in zip(ρ,σ))
    # @show αmin
end
