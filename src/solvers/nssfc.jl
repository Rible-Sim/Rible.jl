
module NSSFC
import ..TensegrityRobots as TR
using Parameters
using LinearAlgebra
using StaticArrays
using BlockDiagonals
using Printf
using ProgressMeter
using OffsetArrays

function initial_guesses(b,qₛ,vₛ,ṽ̇ₛ,aₛ,𝛌bₛ,p,h)
    @unpack αm,αf,γ,β = p
    aₛ₊₁0 = (αf*ṽ̇ₛ - αm*aₛ)/(1-αm)
    ṽₛ₊₁0 = vₛ + h*(1-γ)*aₛ + h*γ*aₛ₊₁0
    qₛ₊₁0 = qₛ + h*vₛ + h^2*(0.5-β)*aₛ + h^2*β*aₛ₊₁0
    vₛ₊₁0 = ṽₛ₊₁0
    ṽ̇ₛ₊₁0 = zero(qₛ₊₁0)
    𝛌bₛ₊₁0 = copy(𝛌bₛ)
    qₛ₊₁0,vₛ₊₁0,ṽₛ₊₁0,ṽ̇ₛ₊₁0,aₛ₊₁0,𝛌bₛ₊₁0
end

function update_smooth(n,b,xe,t,p,h,scaling,dyfuncs)
    qₛ₊₁, vₛ₊₁, ṽₛ₊₁, ṽ̇ₛ₊₁, 𝛌bₛ₊₁, 𝛌uₛ₊₁, 𝚲uₛ₊₁ = xe
    𝐌,𝚽,𝚽𝐪,𝐅,Jacobians,contact_funcs = dyfuncs
    ∂𝚽𝐪T𝛌∂𝒒,∂𝚽𝐪𝐯∂𝒒,∂𝐌𝐚∂𝐪,∂𝐅∂𝐪,∂𝐅∂𝐯 = Jacobians
    𝐠,get_indices,get_FCs,get_D = contact_funcs

    @unpack αm,αf,γ,β,βₜ,γₜ = p
    T = eltype(qₛ₊₁)
    v̂ₛ₊₁ = vₛ₊₁ - ṽₛ₊₁
    Mₛ₊₁ = 𝐌(qₛ₊₁)
    Fₛ₊₁ = 𝐅(qₛ₊₁,vₛ₊₁,t)
    Fqₛ₊₁ = ∂𝐅∂𝐪(qₛ₊₁,vₛ₊₁,t)
    Fvₛ₊₁ = ∂𝐅∂𝐯(qₛ₊₁,vₛ₊₁,t)
    Φₛ₊₁ = 𝚽(qₛ₊₁)
    Φqₛ₊₁ = 𝚽𝐪(qₛ₊₁)
    ∂ΦqTλ∂qₛ₊₁ = ∂𝚽𝐪T𝛌∂𝒒(qₛ₊₁,𝛌bₛ₊₁)
    h² = h^2
    𝐫𝐞𝐬 = vcat(
        Mₛ₊₁*ṽ̇ₛ₊₁ .+ scaling.*transpose(Φqₛ₊₁)*𝛌bₛ₊₁ .- Fₛ₊₁,
                 # - transpose(D̃ₛ₊₁)*𝛌uₛ₊₁,
        Mₛ₊₁*v̂ₛ₊₁, #- transpose(D̂ₛ₊₁)*𝚲uₛ₊₁,
        scaling.*Φₛ₊₁
    )
    Kₜ = scaling.*∂ΦqTλ∂qₛ₊₁ .- Fqₛ₊₁
    𝐉 = [
        γₜ.*Mₛ₊₁.+βₜ.*Kₜ (h/2).*(Kₜ.-Fvₛ₊₁) scaling.*transpose(Φqₛ₊₁);
         -Mₛ₊₁          Mₛ₊₁      zeros(T,n,b);
         scaling.*βₜ.*Φqₛ₊₁        scaling.*(h/2).*Φqₛ₊₁ 0I
    ]
    𝐉,𝐫𝐞𝐬
end

function update_nonsmooth(n,b,u,xe,vₛ,gₙ,Dₛ₊₁,μs,es,t,p,h,scaling,dyfuncs)
    qₛ₊₁, vₛ₊₁, ṽₛ₊₁, ṽ̇ₛ₊₁, 𝛌bₛ₊₁, 𝛌uₛ₊₁, 𝚲uₛ₊₁ = xe
    𝐌,𝚽,𝚽𝐪,𝐅,Jacobians,contact_funcs = dyfuncs
    ∂𝚽𝐪T𝛌∂𝒒,∂𝚽𝐪𝐯∂𝒒,∂𝐌𝐚∂𝐪,∂𝐅∂𝐪,∂𝐅∂𝐯 = Jacobians
    𝐠,get_indices,get_FCs,get_D = contact_funcs

    @unpack αm,αf,γ,β,βₜ,γₜ = p
    T = eltype(qₛ₊₁)
    v̂ₛ₊₁ = vₛ₊₁ - ṽₛ₊₁
    Mₛ₊₁ = 𝐌(qₛ₊₁)
    Fₛ₊₁ = 𝐅(qₛ₊₁,vₛ₊₁,t)
    Fqₛ₊₁ = ∂𝐅∂𝐪(qₛ₊₁,vₛ₊₁,t)
    Fvₛ₊₁ = ∂𝐅∂𝐯(qₛ₊₁,vₛ₊₁,t)
    Φₛ₊₁ = 𝚽(qₛ₊₁)
    Φqₛ₊₁ = 𝚽𝐪(qₛ₊₁)
    ∂ΦqTλ∂qₛ₊₁ = ∂𝚽𝐪T𝛌∂𝒒(qₛ₊₁,𝛌bₛ₊₁)
    # D̃ₛ₊₁ = D̃(qₛ₊₁)
    # D̃qₛ₊₁ = D̃q(qₛ₊₁)
    # D̂qₛ₊₁ = D̂q(qₛ₊₁)
    # Dₛ₊₁ = vcat(D̃ₛ₊₁,D̂ₛ₊₁)
    # Dqₛ₊₁ = vcat(D̃qₛ₊₁,D̂qₛ₊₁)
    D̂ₛ₊₁ = Dₛ₊₁
    ∂Dv∂qₛ₊₁ = zeros(T,size(Dₛ₊₁,1),length(qₛ₊₁))
    h² = h^2
    𝐫𝐞𝐬 = vcat(
        Mₛ₊₁*ṽ̇ₛ₊₁ .+ scaling.*transpose(Φqₛ₊₁)*𝛌bₛ₊₁ .- Fₛ₊₁,
        Mₛ₊₁*v̂ₛ₊₁- transpose(D̂ₛ₊₁)*𝚲uₛ₊₁,
        scaling.*Φₛ₊₁
    )
    # @show transpose(D̂ₛ₊₁)
    Kₜ = scaling.*∂ΦqTλ∂qₛ₊₁ .- Fqₛ₊₁
    𝐉 = [
        γₜ.*Mₛ₊₁.+βₜ.*Kₜ (h/2).*(Kₜ.-Fvₛ₊₁) scaling.*transpose(Φqₛ₊₁);
         -Mₛ₊₁          Mₛ₊₁      zeros(T,n,b);
         scaling.*βₜ.*Φqₛ₊₁        scaling.*(h/2).*Φqₛ₊₁ 0I
    ]
    𝐎nu = zeros(T,n,3u)
    𝐎bu = zeros(T,b,3u)
    𝐁 = [
        𝐎nu;
        transpose(D̂ₛ₊₁);
        𝐎bu
    ]
    𝐜ᵀ = [
        βₜ.*∂Dv∂qₛ₊₁ (h/2).*∂Dv∂qₛ₊₁ .+ Dₛ₊₁ zeros(T,3u,b);
    ]
    𝐜ᵀinv𝐉 = 𝐜ᵀ*inv(𝐉)
    𝐍 = 𝐜ᵀinv𝐉*𝐁
    𝐛 = zeros(T,3u)
    v́ₛ = Dₛ₊₁*vₛ
    v́ₛ₊₁ = Dₛ₊₁*vₛ₊₁
    for i = 1:u
        is = 3(i-1)
        vⁱₛ = v́ₛ[is+1:is+3]
        vⁱₛ₊₁ = v́ₛ₊₁[is+1:is+3]
        vₜⁱₛ₊₁ = norm(vⁱₛ₊₁[2:3])
        vₙⁱₛ = vⁱₛ[1]
        # @show i, vₙⁱₛ, vₜⁱₛ₊₁
        # @show vₛ, vⁱₛ
        # @show Dₛ₊₁
        # @show μs[i],es[i],vₙⁱₛ
        # @show i, 𝚲uₛ₊₁[is+1]
        𝐛[is+1] = μs[i]*vₜⁱₛ₊₁ + es[i]*vₙⁱₛ
    end
    # 𝐫 = (D*v + 𝐛) - 𝐜ᵀinv𝐉*(𝐫𝐞𝐬 - B*vcat(𝛌uₛ₊₁,𝚲uₛ₊₁))
    𝐫 = (Dₛ₊₁*vₛ₊₁ + 𝐛) - 𝐜ᵀinv𝐉*(𝐫𝐞𝐬 + 𝐁*𝚲uₛ₊₁)
    # @show Dₛ₊₁*vₛ₊₁

    𝐉,𝐫𝐞𝐬,𝐁,𝐜ᵀ,𝐍,𝐫
    # 𝐉,𝐫𝐞𝐬
end

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
    end
    output .= 𝛄ₕ
    return rmin
end

function nssfc(n,b,q0,v0,p,h,dyfuncs,tspan;tol=1e-14,imax=100)
    𝐌,𝚽,𝚽𝐪,𝐅,Jacobians,contact_funcs = dyfuncs
    ∂𝚽𝐪T𝛌∂𝒒,∂𝚽𝐪𝐯∂𝒒,∂𝐌𝐚∂𝐪,∂𝐅∂𝐪,∂𝐅∂𝐯 = Jacobians
    𝐠,get_indices,get_FCs,get_D = contact_funcs
    @unpack αm,αf,γ,β,γₜ,βₜ = p
    tstart,tend = tspan
    totaltime = tend - tstart
    totalstep = round(Int,totaltime/h,RoundUp)
    M0 = 𝐌(q0)
    F0 = 𝐅(q0,v0,tstart)
    ṽ̇0 = M0\F0
    a0 = ṽ̇0
    mr = norm(M0,Inf)
    scaling = mr/h
    # @show ṽ̇0
    ts = collect(tstart:h:tend)
    qs = [copy(q0) for i = 1:totalstep+1]
    vs = [copy(v0) for i = 1:totalstep+1]
    as = [copy(a0) for i = 1:totalstep+1]
    ṽ̇s = [copy(ṽ̇0) for i = 1:totalstep+1]
    λbs = [zeros(eltype(q0),b) for i = 1:totalstep+1]
    qˣ = copy(q0)
    nx = n+n+b
    Δx = ones(eltype(q0),nx)
    J = zeros(eltype(q0),nx,nx)
    res = zeros(eltype(q0),nx)
    converged = true
    for timestep = 1:totalstep
        tₛ = ts[timestep]
        qₛ = qs[timestep]
        vₛ = vs[timestep]
        aₛ = as[timestep]
        ṽ̇ₛ = ṽ̇s[timestep]
        λbₛ = λbs[timestep]
        tₛ₊₁ = ts[timestep+1]
        qₛ₊₁,vₛ₊₁,ṽₛ₊₁,ṽ̇ₛ₊₁,aₛ₊₁,𝛌bₛ₊₁ = initial_guesses(b,qₛ,vₛ,ṽ̇ₛ,aₛ,λbₛ,p,h)
        qˣ .= qₛ₊₁
        u,active_indices,g = get_indices(qˣ)
        gₙ = g[active_indices]
        𝛌uₛ₊₁ = zeros(eltype(qₛ₊₁),3u)
        𝚲uₛ₊₁ = zeros(eltype(qₛ₊₁),3u)
        𝚲uʳₛ₊₁ = copy(𝚲uₛ₊₁)
        # Newton iteration
        invM = inv(𝐌(qₛ₊₁))
        for iteration = 1:imax
            xe = (qₛ₊₁, vₛ₊₁, ṽₛ₊₁, ṽ̇ₛ₊₁, 𝛌bₛ₊₁, 𝛌uₛ₊₁, 𝚲uₛ₊₁)

            if u == 0
                # @warn "Smooth, timestep=$timestep"
                𝐉,𝐫𝐞𝐬 = update_smooth(n,b,xe,tₛ₊₁,p,h,scaling,dyfuncs)
                Δx = 𝐉\(-𝐫𝐞𝐬)
            else
                @show timestep,iteration,u,active_indices
                Dₛ₊₁,μs,es = get_D(active_indices,qₛ₊₁)
                # @show Dₛ₊₁
                # _,_,g = get_indices(qₛ₊₁)
                gₙ = g[active_indices]
                𝐉,𝐫𝐞𝐬,𝐁,𝐜ᵀ,𝐍,𝐫 = update_nonsmooth(n,b,u,xe,vₛ,gₙ,Dₛ₊₁,μs,es,tₛ₊₁,p,h,scaling,dyfuncs)
                B = make_B(u,Dₛ₊₁,invM)
                r = make_residual4(μs,𝐍,𝐫)
                # 𝚲uₛ₊₁,_ = Jacobi(B,r,μs,𝐍,𝐫;τ=1e-13,Nmax=1000)
                # 𝚲uₛ₊₁,GS_k,GS_res = GaussSeidel(u,B,r,μs,𝐍,𝐫)
                # @show GS_k,GS_res
                # 𝚲uₛ₊₁,_ = APGD(r,μs,𝐍,𝐫)
                # @show 𝚲uₛ₊₁, vₛ₊₁
                𝚲uₛ₊₁,_  = APGD(r,μs,𝐍,𝐫;τ=1e-10,Nmax=1000)
                @show sum(𝚲uₛ₊₁/h)
                Δ𝚲uₛ₊₁ = 𝚲uₛ₊₁ - 𝚲uʳₛ₊₁
                # @show iteration,abs(𝚲uₛ₊₁[3]/𝚲uₛ₊₁[1])
                # @show iteration,𝚲uʳₛ₊₁
                # @show iteration,Δ𝚲uₛ₊₁
                Δx = 𝐉\(-𝐫𝐞𝐬 + 𝐁*(Δ𝚲uₛ₊₁))
                𝚲uʳₛ₊₁ .= 𝚲uₛ₊₁
            end
            normΔx = norm(Δx)
            # @show normΔx
            if  normΔx < tol
                break
            end
            if iteration == imax
                # converged = false
                @error "Newton: Max iteration $iteration, res=$(normΔx)"
                @info "timestep=$timestep, u=$u,active_indices=$active_indices"
            end
            Δṽ,Δv,Δ𝛌b = TR.split_by_lengths(Δx,[n,n,b])
            ṽₛ₊₁ += Δṽ
            vₛ₊₁ = vₛ₊₁ + Δv
            ṽ̇ₛ₊₁ += γₜ*Δṽ
            qₛ₊₁ += βₜ*Δṽ + 1/2*h*Δv
            𝛌bₛ₊₁ += Δ𝛌b
            # @show ṽ̇ₛ₊₁,vₛ₊₁,ṽₛ₊₁
            # @show 𝚲uₛ₊₁
        end
        if !converged
            break
        end
        aₛ₊₁ += (1-αf)/(1-αm)*ṽ̇ₛ₊₁
        qs[timestep+1] .= qₛ₊₁
        vs[timestep+1] .= vₛ₊₁
        as[timestep+1] .= aₛ₊₁
        ṽ̇s[timestep+1] .= ṽ̇ₛ₊₁
    end
    ts,qs,vs,as,ṽ̇s
end


function Momentum_k(qₛ₋₁,pₛ₋₁,qₛ,λₛ,μₛ,M,A,B,h)
    pₛ = -pₛ₋₁ .+ 2/h.*M*(qₛ.-qₛ₋₁) .+
        1/(2h).*(transpose(A(qₛ))-transpose(A(qₛ₋₁)))*λₛ .+
        1/(2h).*(transpose(B(qₛ))-transpose(B(qₛ₋₁)))*μₛ
end

function stepk_maker(nq,nλ,nμ,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,h)
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dyfuncs
    Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q = jacobians
    # E,𝐠,𝐠𝐪,∂𝐠𝐪ᵀΛ∂q,∂𝐠𝐪q̇∂q = contact_funcs
    n1 =  nq
    n2 = nq+nλ
    n3 = nq+nλ+nμ
    function stepk!(R,J,x,q̇ₛ)
        qₛ = @view x[   1:n1]
        λₛ = @view x[n1+1:n2]
        μₛ = @view x[n2+1:n3]

        q = (qₛ.+qₛ₋₁)./2
        q̇ = (qₛ.-qₛ₋₁)./h
        t = tₛ₋₁+h/2

        F⁺ = zeros(eltype(qₛ),nq)
        F!(F⁺,q,q̇,t)
        ∂F∂q = zeros(eltype(qₛ),nq,nq)
        ∂F∂q̇ = zeros(eltype(qₛ),nq,nq)
        Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)

        # q̇ₛ = invM*Momentum_k(qₛ₋₁,pₛ₋₁,qₛ,λₛ,μₛ,M,A,B,h)

        Aₛ₋₁ = A(qₛ₋₁)
        Bₛ₋₁ = B(qₛ₋₁)
        Aₛ = A(qₛ)
        Bₛ = B(qₛ)
        ∂q̇ₛ∂qₛ = 2/h*I + 1/(2h).*invM*(∂Aᵀλ∂q(qₛ,λₛ) + ∂Bᵀμ∂q(qₛ,μₛ))
        ∂q̇ₛ∂λₛ = invM*transpose(Aₛ-Aₛ₋₁)/(2h)
        ∂q̇ₛ∂μₛ = invM*transpose(Bₛ-Bₛ₋₁)/(2h)

        R[   1:n1] .= -h.*pₛ₋₁ .+ M*(qₛ.-qₛ₋₁) .-
                        1/2 .*transpose(Aₛ₋₁)*λₛ .-
                        1/2 .*transpose(Bₛ₋₁)*μₛ .-
                        (h^2)/2 .*F⁺
        R[n1+1:n2] .= Φ(qₛ)
        R[n2+1:n3] .= Ψ(qₛ,q̇ₛ)

        J .= 0.0
        J[   1:n1,   1:n1] .=  M .-h^2/2 .*(1/2 .*∂F∂q .+ 1/h.*∂F∂q̇)
        J[   1:n1,n1+1:n2] .= -1/2 .*transpose(Aₛ₋₁)
        J[   1:n1,n2+1:n3] .= -1/2 .*transpose(Bₛ₋₁)

        J[n1+1:n2,   1:n1] .=  Aₛ

        J[n2+1:n3,   1:n1] .=  Ψq(qₛ,q̇ₛ)+Bₛ*∂q̇ₛ∂qₛ
        J[n2+1:n3,n1+1:n2] .=  Bₛ*∂q̇ₛ∂λₛ
        J[n2+1:n3,n2+1:n3] .=  Bₛ*∂q̇ₛ∂μₛ
    end

    stepk!
end

function ns_stepk_maker(nq,nλ,nμ,nu,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,h)
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dyfuncs
    Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q = jacobians
    𝐠,get_indices,get_FCs,get_D = contact_funcs

    stepk! = stepk_maker(nq,nλ,nμ,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,h)

    n1 =  nq
    n2 = 2nq
    n3 = 2nq+nλ
    n4 = 2nq+nλ+nμ
    nΛ = 3nu
    nx = n4
    T = eltype(qₛ₋₁)
    𝐁 = zeros(T,nx,nΛ)
    𝐜ᵀ = zeros(T,nΛ,nx)
    # 𝐜ᵀinv𝐉 = zeros(T,nΛ,nx)
    𝐍 = zeros(T,nΛ,nΛ)
    𝐫 = zeros(T,nΛ)
    function ns_stepk!(𝐫𝐞𝐬,𝐉,𝐁,𝐜ᵀ,𝐍,𝐫,x,𝚲ₛ,vₛ,Dₛ,ηs,es,H)
        q̃ₛ = @view x[   1:n1]
        qₛ = @view x[n1+1:n2]
        λₛ = @view x[n2+1:n3]
        μₛ = @view x[n3+1:n4]
        vₛ₋₁ = q̇ₛ₋₁

        q = (qₛ.+qₛ₋₁)./2
        q̇ = (qₛ.-qₛ₋₁)./h
        t = tₛ₋₁+h/2

        F⁺ = zeros(eltype(qₛ),nq)
        F!(F⁺,q,q̇,t)
        ∂F∂q = zeros(eltype(qₛ),nq,nq)
        ∂F∂q̇ = zeros(eltype(qₛ),nq,nq)
        Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)

        Aₛ₋₁ = A(qₛ₋₁)
        Bₛ₋₁ = B(qₛ₋₁)
        Aₛ = A(qₛ)
        Bₛ = B(qₛ)
        ∂vₛ∂qₛ = 2/h*I + 1/(2h).*invM*(∂Aᵀλ∂q(qₛ,λₛ) + ∂Bᵀμ∂q(qₛ,μₛ))
        ∂vₛ∂λₛ = invM*transpose(Aₛ-Aₛ₋₁)/(2h)
        ∂vₛ∂μₛ = invM*transpose(Bₛ-Bₛ₋₁)/(2h)

        ∂DₛᵀHΛₛ∂qₛ = zeros(eltype(x),nq,nq)

        𝐫𝐞𝐬[   1:n1] .= -h.*pₛ₋₁ .+ M*(q̃ₛ.-qₛ₋₁) .-
                        1/2 .*transpose(Aₛ₋₁)*λₛ .-
                        1/2 .*transpose(Bₛ₋₁)*μₛ .-
                        (h^2)/2 .*F⁺
        𝐫𝐞𝐬[n1+1:n2] .= (2/h).*M*(qₛ - q̃ₛ) - transpose(Dₛ)*H*𝚲ₛ
        𝐫𝐞𝐬[n2+1:n3] .= Φ(qₛ)
        𝐫𝐞𝐬[n3+1:n4] .= Ψ(qₛ,vₛ)

        𝐉 .= 0.0
        𝐉[   1:n1,   1:n1] .=  M
        𝐉[   1:n1,n1+1:n2] .= -h^2/2 .*(1/2 .*∂F∂q .+ 1/h.*∂F∂q̇)
        𝐉[   1:n1,n2+1:n3] .= -1/2 .*transpose(Aₛ₋₁)
        𝐉[   1:n1,n3+1:n4] .= -1/2 .*transpose(Bₛ₋₁)

        𝐉[n1+1:n2,   1:n1] .= -(2/h).*M
        𝐉[n1+1:n2,n1+1:n2] .=  (2/h).*M .- ∂DₛᵀHΛₛ∂qₛ

        𝐉[n2+1:n3,n1+1:n2] .=  Aₛ

        𝐉[n3+1:n4,n1+1:n2] .=  Ψq(qₛ,vₛ) .+ Bₛ*∂vₛ∂qₛ
        𝐉[n3+1:n4,n2+1:n3] .=  Bₛ*∂vₛ∂λₛ
        𝐉[n3+1:n4,n3+1:n4] .=  Bₛ*∂vₛ∂μₛ

        𝐁 .= 0
        𝐁[n1+1:n2,1:nΛ] .= transpose(Dₛ)*H

        𝐛 = zeros(eltype(x),nΛ)
        v́ₛ₋₁ = Dₛ*vₛ₋₁
        v́ₛ = Dₛ*vₛ
        for i = 1:nu
            is = 3(i-1)
            vⁱₛ₋₁ = v́ₛ₋₁[is+1:is+3]
            vⁱₛ   = v́ₛ[is+1:is+3]
            vₜⁱₛ = norm(vⁱₛ[2:3])
            vₙⁱₛ₋₁ = vⁱₛ₋₁[1]
            vₙⁱₛ = vⁱₛ[1]
            # @show vₜⁱₛ, vₙⁱₛ₋₁, vₙⁱₛ
            𝐛[is+1] = ηs[i]*vₜⁱₛ + es[i]*vₙⁱₛ₋₁
        end

        ∂Dₛvₛ∂qₛ = zeros(eltype(x),nΛ,nq)
        ∂bₛ∂qₛ = zeros(eltype(x),nΛ,nq)

        𝐜ᵀ[1:nΛ,   1:n1] .= 0
        𝐜ᵀ[1:nΛ,n1+1:n2] .= ∂Dₛvₛ∂qₛ + ∂bₛ∂qₛ + Dₛ*∂vₛ∂qₛ
        𝐜ᵀ[1:nΛ,n2+1:n3] .= Dₛ*∂vₛ∂λₛ
        𝐜ᵀ[1:nΛ,n3+1:n4] .= Dₛ*∂vₛ∂μₛ

        𝐜ᵀinv𝐉 = 𝐜ᵀ*inv(𝐉)
        𝐍 .= 𝐜ᵀinv𝐉*𝐁
        @show Dₛ*vₛ + 𝐛
        𝐫 .= (Dₛ*vₛ + 𝐛) - 𝐜ᵀinv𝐉*(𝐫𝐞𝐬 + 𝐁*𝚲ₛ)
    end
    ns_stepk!,𝐁,𝐜ᵀ,𝐍,𝐫
end

function ⊙(x::AbstractVector{T},y::AbstractVector{T}) where {T<:Real}
    x0 = x[begin]
    x1 = @view x[begin+1:end]
    y0 = y[begin]
    y1 = @view y[begin+1:end]
    vcat(transpose(x)*y, x0.*y1 .+ y0.*x1)
end

function ⊘(x::AbstractVector{T},y::AbstractVector{T}) where {T<:Real}
    x0 = x[begin]
    x1 = @view x[begin+1:end]
    X = hvcat(
            (2,2),
            x0,transpose(x1),
            x1,Diagonal(x0.*one.(x1))
        )
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

function find_cone_step_length(z_split,W_blocks,Δy_split,Δx_split,J)
    Δx̃_split = W_blocks.*Δx_split
    Δỹ_split = inv.(W_blocks).*Δy_split
    z_split_norm = [sqrt(transpose(zi)*J*zi) for zi in z_split]
    z̄_split = z_split./z_split_norm
    ρ = [
        begin
            head = transpose(z̄i)*J*Δỹi
            tail = Δỹi[begin+1:end] .- (head+Δỹi[begin])./(z̄i[begin] + 1).*z̄i[begin+1:end]
            vcat(head,tail)/zi_norm
        end
        for (z̄i,Δỹi,zi_norm) in zip(z̄_split,Δỹ_split,z_split_norm)
    ]

    σ = [
        begin
            head = transpose(z̄i)*J*Δx̃i
            tail = Δx̃i[begin+1:end] .- (head+Δx̃i[begin])./(z̄i[begin] + 1).*z̄i[begin+1:end]
            vcat(head,tail)/zi_norm
        end
        for (z̄i,Δx̃i,zi_norm) in zip(z̄_split,Δx̃_split,z_split_norm)
    ]

    αmin = minimum(max(0,norm(ρi[begin+1:end])-ρi[begin],norm(σi[begin+1:end])-σi[begin])^(-1) for (ρi,σi) in zip(ρ,σ))
end


function ip_ns_stepk_maker(nq,nλ,nμ,nu,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,h)
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dyfuncs
    Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q = jacobians
    𝐠,get_indices,get_FCs,get_D = contact_funcs

    stepk! = stepk_maker(nq,nλ,nμ,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,h)

    n1 = nq
    n2 = n1 + nq
    n3 = n2 + nλ
    n4 = n3 + nμ
    nΛ = 3nu
    n5 = n4 + nΛ
    nx = n5
    T = eltype(qₛ₋₁)
    e = [one(T),zero(T),zero(T)]
    J = Diagonal([one(T),-one(T),-one(T)])
    𝐞_split = [e for i = 1:nu]
    function ip_ns_stepk!(𝐫𝐞𝐬,𝐉,x,yₛ,vₛ,Dₛ,ηs,es,H,μ)
        q̃ₛ = @view x[   1:n1]
        qₛ = @view x[n1+1:n2]
        λₛ = @view x[n2+1:n3]
        μₛ = @view x[n3+1:n4]
        Λₛ = @view x[n4+1:n5]

        𝐉 .= 0.0
        vₛ₋₁ = q̇ₛ₋₁

        q = (qₛ.+qₛ₋₁)./2
        q̇ = (qₛ.-qₛ₋₁)./h
        t = tₛ₋₁+h/2

        F⁺ = zeros(eltype(qₛ),nq)
        F!(F⁺,q,q̇,t)
        ∂F∂q = zeros(eltype(qₛ),nq,nq)
        ∂F∂q̇ = zeros(eltype(qₛ),nq,nq)
        Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)

        Aₛ₋₁ = A(qₛ₋₁)
        Bₛ₋₁ = B(qₛ₋₁)
        Aₛ = A(qₛ)
        Bₛ = B(qₛ)
        ∂vₛ∂qₛ = 2/h*I + 1/(2h).*invM*(∂Aᵀλ∂q(qₛ,λₛ) + ∂Bᵀμ∂q(qₛ,μₛ))
        ∂vₛ∂λₛ = invM*transpose(Aₛ-Aₛ₋₁)/(2h)
        ∂vₛ∂μₛ = invM*transpose(Bₛ-Bₛ₋₁)/(2h)

        ∂Dₛvₛ∂qₛ = zeros(eltype(x),nΛ,nq)
        ∂bₛ∂qₛ = zeros(eltype(x),nΛ,nq)
        ∂DₛᵀHΛₛ∂qₛ = zeros(eltype(x),nq,nq)

        𝐛 = zeros(eltype(x),nΛ)
        ∂𝐛∂𝐱 = @view 𝐉[n4+1:n5,1:nx]
        v́ₛ₋₁ = Dₛ*vₛ₋₁
        v́ₛ = Dₛ*vₛ
        for i = 1:nu
            is = 3(i-1)
            vⁱₛ₋₁ = v́ₛ₋₁[is+1:is+3]
            vⁱₛ   = v́ₛ[is+1:is+3]
            vₜⁱₛ = norm(vⁱₛ[2:3])
            vₙⁱₛ₋₁ = vⁱₛ₋₁[1]
            vₙⁱₛ = vⁱₛ[1]
            # @show vₜⁱₛ, vₙⁱₛ₋₁, vₙⁱₛ
            𝐛[is+1] = ηs[i]*vₜⁱₛ + es[i]*vₙⁱₛ₋₁
            D̃i = @view Dₛ[is+1:is+3,:]
            D̃i[1,:] += (vⁱₛ[2]*Dₛ[is+2,:]+vⁱₛ[3]*Dₛ[is+3,:])/vₜⁱₛ
            ∂𝐛∂𝐱[is+1:is+3,n1+1:n2] .= D̃i*∂vₛ∂qₛ
            ∂𝐛∂𝐱[is+1:is+3,n2+1:n3] .= D̃i*∂vₛ∂λₛ
            ∂𝐛∂𝐱[is+1:is+3,n3+1:n4] .= D̃i*∂vₛ∂μₛ
        end
        # @show "before",yₛ,Λₛ
        y_split = TR.split_by_lengths(yₛ,3)
        Λ_split = TR.split_by_lengths(Λₛ,3)
        W_blocks = NTScale.(y_split,Λ_split)
        ΘW_blocks = NTScale_Anderson.(y_split,Λ_split)
        # @show W_blocks
        # @show ΘW_blocks
        z_split = W_blocks.*y_split
        z = reduce(vcat,z_split)
        W = BlockDiagonal(W_blocks)
        WᵀW = transpose(W)*W
        @show z_split
        @show z_split⊙z_split
        𝐫𝐞𝐬[   1:n1] .= -h.*pₛ₋₁ .+ M*(q̃ₛ.-qₛ₋₁) .-
                        1/2 .*transpose(Aₛ₋₁)*λₛ .-
                        1/2 .*transpose(Bₛ₋₁)*μₛ .-
                        (h^2)/2 .*F⁺
        𝐫𝐞𝐬[n1+1:n2] .= (2/h).*M*(qₛ - q̃ₛ) - transpose(Dₛ)*H*Λₛ
        𝐫𝐞𝐬[n2+1:n3] .= Φ(qₛ)
        𝐫𝐞𝐬[n3+1:n4] .= Ψ(qₛ,vₛ)
        𝐫𝐞𝐬[n4+1:n5] .= Dₛ*vₛ + 𝐛
        # res = norm(𝐫𝐞𝐬)
        @show (2/h).*M*(qₛ - q̃ₛ)
        @show transpose(Dₛ)*H*Λₛ
        @show Λₛ
        @show 𝐫𝐞𝐬
        𝐉[   1:n1,   1:n1] .=  M
        𝐉[   1:n1,n1+1:n2] .= -h^2/2 .*(1/2 .*∂F∂q .+ 1/h.*∂F∂q̇)
        𝐉[   1:n1,n2+1:n3] .= -1/2 .*transpose(Aₛ₋₁)
        𝐉[   1:n1,n3+1:n4] .= -1/2 .*transpose(Bₛ₋₁)

        𝐉[n1+1:n2,   1:n1] .= -(2/h).*M
        𝐉[n1+1:n2,n1+1:n2] .=  (2/h).*M .- ∂DₛᵀHΛₛ∂qₛ
        𝐉[n1+1:n2,n4+1:n5] .= -transpose(Dₛ)*H

        𝐉[n2+1:n3,n1+1:n2] .=  Aₛ

        𝐉[n3+1:n4,n1+1:n2] .=  Ψq(qₛ,vₛ) .+ Bₛ*∂vₛ∂qₛ
        𝐉[n3+1:n4,n2+1:n3] .=  Bₛ*∂vₛ∂λₛ
        𝐉[n3+1:n4,n3+1:n4] .=  Bₛ*∂vₛ∂μₛ

        𝐉[n4+1:n5,n4+1:n5] .= WᵀW

        η = 1.0
        Δxp = 𝐉\(-𝐫𝐞𝐬)
        # @show 𝐉
        # @show 𝐫𝐞𝐬
        # @show Δxp
        ΔΛp = @view Δxp[n4+1:n5]
        Δyp = -yₛ-WᵀW*ΔΛp
        ΔΛp_split = TR.split_by_lengths(ΔΛp,3)
        Δyp_split = TR.split_by_lengths(Δyp,3)
        # @show ΔΛp_split, Δyp_split
        αmax = find_cone_step_length(z_split,W_blocks,Δyp_split,ΔΛp_split,J)
        α = 1#min(one(αmax),0.99αmax)
        yp_split = y_split .+ α.*Δyp_split
        Λp_split = Λ_split .+ α.*ΔΛp_split
        # yp_cone = [transpose(yi)*J*yi for yi in yp_split]
        # Λp_cone = [transpose(Λi)*J*Λi for Λi in Λp_split]
        # @show yp_cone
        # @show Λp_cone
        yp = reduce(vcat,yp_split)
        Λp = reduce(vcat,Λp_split)
        μp = transpose(yp)*Λp/nΛ
        σ = (μp/μ)^3
        τ = σ*μp
        # @show αmax,α,τ,σ,μ,μp
        @show α
        @show Δxp
        @show ΔΛp,Δyp
        𝐫𝐞𝐬_c_split = W_blocks.*(z_split⊘(-σ.*μp.*𝐞_split.+((inv.(W_blocks).*Δyp_split)⊙(W_blocks.*ΔΛp_split))))
        𝐫𝐞𝐬_c = reduce(vcat,𝐫𝐞𝐬_c_split)
        𝐫𝐞𝐬[n4+1:n5] .+= 𝐫𝐞𝐬_c
        Δxc = 𝐉\(-𝐫𝐞𝐬)
        # η = exp(-0.1μ) + 0.9
        ΔΛc = @view Δxc[n4+1:n5]
        Δyc = -yₛ-𝐫𝐞𝐬_c-WᵀW*ΔΛp
        ΔΛc_split = TR.split_by_lengths(ΔΛc,3)
        Δyc_split = TR.split_by_lengths(Δyc,3)
        αmax = find_cone_step_length(z_split,W_blocks,Δyc_split,ΔΛc_split,J)
        α = min(1,0.99αmax)
        # @show α
        q̃ₛ .+= α.*Δxc[   1:n1]
        qₛ .+= α.*Δxc[n1+1:n2]
        λₛ .+= α.*Δxc[n2+1:n3]
        μₛ .+= α.*Δxc[n3+1:n4]
        Λ_split .+= α.*ΔΛc_split
        y_split .+= α.*Δyc_split
        # @show ΔΛc_split, Δyc_split
        y_cone = [transpose(yi)*J*yi for yi in y_split]
        Λ_cone = [transpose(Λi)*J*Λi for Λi in Λ_split]
        # @show y_cone
        # @show Λ_cone
        Λₛ .= reduce(vcat,Λ_split)
        yₛ .= reduce(vcat,y_split)
        # @show "after",yₛ,Λₛ
        μ = transpose(yₛ)*Λₛ/nΛ
        @show σ,μ,μp
        μ,Δxc
    end
    ip_ns_stepk!
end

function nhsolve(nq,nλ,nμ,q0,q̇0,dyfuncs,tspan;dt=0.01,ftol=1e-14,xtol=ftol,verbose=false,imax=50,
                progress=true,exception=true)
    # @unpack bot,tspan,dyfuncs,control!,restart = prob
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dyfuncs
    Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q = jacobians
    𝐠,get_indices,get_FCs,get_D = contact_funcs
    totaltime = tspan[end] - tspan[begin]
    totalstep = ceil(Int,totaltime/dt)
    cs = OffsetArray([-1 for i in 1:totalstep+1],0:totalstep)
    ts = OffsetArray([tspan[begin]+(i-1)*dt for i in 1:totalstep+1],0:totalstep)
    q̃s = OffsetArray([copy(q0) for i in 1:totalstep+1],0:totalstep)
    qs = OffsetArray([copy(q0) for i in 1:totalstep+1],0:totalstep)
    q̇s = OffsetArray([copy(q̇0) for i in 1:totalstep+1],0:totalstep)
    ps = OffsetArray([M*copy(q̇0) for i in 1:totalstep+1],0:totalstep)
    λs = OffsetArray([zeros(eltype(q0),nλ) for i in 1:totalstep+1],0:totalstep)
    μs = OffsetArray([zeros(eltype(q0),nμ) for i in 1:totalstep+1],0:totalstep)
    # Λs = [zeros(eltype(q0),nu) for i in 1:totalstep+1]
    invM = inv(M)
    # F⁺ = zero(q0)
    step = 0
    smooth_nx = nq + nλ + nμ
    nonsmooth_nx = nq + nq + nλ + nμ
    smooth_Δx = zeros(eltype(q0),smooth_nx)
    nonsmooth_Δx = zeros(eltype(q0),nonsmooth_nx)
    smooth_x = zero(smooth_Δx)
    smooth_R = zero(smooth_Δx)
    smooth_J = zeros(eltype(smooth_x),smooth_nx,smooth_nx)
    nonsmooth_x = zero(nonsmooth_Δx)
    nonsmooth_R = zero(nonsmooth_Δx)
    nonsmooth_J = zeros(eltype(nonsmooth_x),nonsmooth_nx,nonsmooth_nx)
    mr = norm(M,Inf)
    scaling = 1

    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        # control!(intor,cache)
        #---------Step k Control-----------
        q̃ₛ₋₁ = q̃s[timestep-1]
        qₛ₋₁ = qs[timestep-1]
        q̇ₛ₋₁ = q̇s[timestep-1]
        pₛ₋₁ = ps[timestep-1]
        λₛ₋₁ = λs[timestep-1]
        μₛ₋₁ = μs[timestep-1]
        tₛ₋₁ = ts[timestep-1]
        q̃ₛ   = q̃s[timestep]
        qₛ   = qs[timestep]
        q̇ₛ   = q̇s[timestep]
        pₛ   = ps[timestep]
        λₛ   = λs[timestep]
        μₛ   = μs[timestep]
        qˣ = qₛ₋₁ .+ dt./2 .*q̇ₛ₋₁
        qₛ .= qˣ
        q̇ₛ .= q̇ₛ₋₁
        nu,active_indices,g = get_indices(qˣ)
        gₙ = g[active_indices]
        cs[timestep] = nu
        # ns_stepk! = ns_stepk_maker(nq,nλ,nμ,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,Aset,dyfuncs,invM,dt)
        # stepk! = stepk_maker(nq,nλ,nμ,nu,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,dt)
        isconverged = false
        res = typemax(eltype(qₛ))
        iteration_break = 0
        APGD_res = typemax(eltype(qₛ))
        if nu == 0
            smooth_x[      1:nq]          .= qₛ
            smooth_x[   nq+1:nq+nλ]       .= 0.0
            smooth_x[nq+nλ+1:nq+nλ+nμ]    .= 0.0
            stepk! = stepk_maker(nq,nλ,nμ,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,dt)

            for iteration = 1:imax
                    stepk!(smooth_R,smooth_J,smooth_x,q̇ₛ)
                    res = norm(smooth_R)
                    if res < ftol
                        isconverged = true
                        iteration_break = iteration-1
                        break
                    end
                    smooth_Δx .= -smooth_J\smooth_R
                    smooth_x .+= smooth_Δx
                    qₛ .= smooth_x[      1:nq]
                    λₛ .= smooth_x[   nq+1:nq+nλ]
                    μₛ .= smooth_x[nq+nλ+1:nq+nλ+nμ]
                    pₛ .= Momentum_k(qₛ₋₁,pₛ₋₁,qₛ,λₛ,μₛ,M,A,B,dt)
                    q̇ₛ .= invM*pₛ
            end
            q̃ₛ .= smooth_x[      1:nq]
        else # u!=0
            nonsmooth_x[         1:nq]          .= qₛ
            nonsmooth_x[      nq+1:nq+nq]       .= qₛ
            nonsmooth_x[   nq+nq+1:nq+nq+nλ]    .= 0.0
            nonsmooth_x[nq+nq+nλ+1:nq+nq+nλ+nμ] .= 0.0
            isconverged = false
            𝚲ₛ = zeros(eltype(qₛ₋₁),3nu)
            𝚲ʳₛ = copy(𝚲ₛ)
            Δ𝚲ₛ = copy(𝚲ₛ)
            # @show timestep, nu
            ns_stepk!,𝐁,𝐜ᵀ,𝐍,𝐫 = ns_stepk_maker(nq,nλ,nμ,nu,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,dt)
            for iteration = 1:imax
                Dₛ,ηs,es,H = get_D(active_indices,qₛ)
                # ηs .= 1
                _,_,g = get_indices(qₛ)
                gₙ = g[active_indices]
                # @show iteration,Dₛ,ηs,es,gₙ
                ns_stepk!(nonsmooth_R,nonsmooth_J,
                            𝐁,𝐜ᵀ,𝐍,𝐫,nonsmooth_x,𝚲ₛ,q̇ₛ,Dₛ,ηs,es,H)
                res = norm(nonsmooth_R)
                # if res < ftol
                #     isconverged = true
                #     iteration_break = iteration-1
                #     break
                # end

                r4 = make_residual4(ηs,𝐍,𝐫;gd=1e-3)
                # Jacobi_B = make_B(nu,Dₛ,invM)
                # 𝚲ₛ,_ = Jacobi(Jacobi_B,r,ηs,𝐍,𝐫;τ=1e-13,Nmax=1000)
                # 𝚲uₛ₊₁,GS_k,GS_res = GaussSeidel(u,B,r,ηs,𝐍,𝐫)
                APGD_res = APGD!(𝚲ₛ,r4,ηs,𝐍,𝐫;τ=1e-10,Nmax=100)
                # @show APGD_res
                Δ𝚲ₛ .= 𝚲ₛ - 𝚲ʳₛ
                # @show 𝚲ₛ, 𝚲ʳₛ
                nonsmooth_Δx .= nonsmooth_J\(-nonsmooth_R + 𝐁*(Δ𝚲ₛ))
                𝚲ʳₛ .= 𝚲ₛ
                nonsmooth_x .+= nonsmooth_Δx
                qₛ .= nonsmooth_x[         nq+1:nq+nq]
                λₛ .= nonsmooth_x[      nq+nq+1:nq+nq+nλ]
                μₛ .= nonsmooth_x[   nq+nq+nλ+1:nq+nq+nλ+nμ]
                pₛ .= Momentum_k(qₛ₋₁,pₛ₋₁,qₛ,λₛ,μₛ,M,A,B,dt)
                q̇ₛ .= invM*pₛ
                normΔx = norm(nonsmooth_Δx)
                res = normΔx
                # @show normΔx, norm(Δ𝚲ₛ)
                iteration_break = iteration
                if  normΔx < xtol
                    isconverged = true
                    break
                end
            end
            # @show gₙ*9.81
            # @show 𝚲ₛ./dt
            q̃ₛ .= nonsmooth_x[      1:nq]
            @show q̃ₛ
            @show qₛ
            @show q̇ₛ
            @show 𝚲ₛ
        end

        if !isconverged
            @warn "NLsolve max iterations $iteration_break, at timestep=$timestep, Res=$(res)"
            if exception
                error("Not converged!")
            else
                # intor.convergence = false
                # break
            end
        else
            if nu == 0
                stepstring = "Smooth"
                # @info "$stepstring timestep=$timestep, iterations=$iteration_break"
            else
                stepstring = "Nonsmooth with $nu contact(s)"
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
    ts,cs,qs,q̇s,ps,λs,μs
end

function ipsolve(nq,nλ,nμ,q0,q̇0,dyfuncs,tspan;dt=0.01,ftol=1e-14,xtol=ftol,verbose=false,imax=50,
                progress=true,exception=true)
    # @unpack bot,tspan,dyfuncs,control!,restart = prob
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dyfuncs
    Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q = jacobians
    𝐠,get_indices,get_FCs,get_D = contact_funcs
    totaltime = tspan[end] - tspan[begin]
    totalstep = ceil(Int,totaltime/dt)
    cs = OffsetArray([-1 for i in 1:totalstep+1],0:totalstep)
    ts = OffsetArray([tspan[begin]+(i-1)*dt for i in 1:totalstep+1],0:totalstep)
    q̃s = OffsetArray([copy(q0) for i in 1:totalstep+1],0:totalstep)
    qs = OffsetArray([copy(q0) for i in 1:totalstep+1],0:totalstep)
    q̇s = OffsetArray([copy(q̇0) for i in 1:totalstep+1],0:totalstep)
    ps = OffsetArray([M*copy(q̇0) for i in 1:totalstep+1],0:totalstep)
    λs = OffsetArray([zeros(eltype(q0),nλ) for i in 1:totalstep+1],0:totalstep)
    μs = OffsetArray([zeros(eltype(q0),nμ) for i in 1:totalstep+1],0:totalstep)
    # Λs = [zeros(eltype(q0),nu) for i in 1:totalstep+1]
    invM = inv(M)
    # F⁺ = zero(q0)
    step = 0
    smooth_nx = nq + nλ + nμ
    smooth_Δx = zeros(eltype(q0),smooth_nx)
    smooth_x = zero(smooth_Δx)
    smooth_R = zero(smooth_Δx)
    smooth_J = zeros(eltype(smooth_x),smooth_nx,smooth_nx)
    mr = norm(M,Inf)
    scaling = 1

    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        # control!(intor,cache)
        #---------Step k Control-----------
        q̃ₛ₋₁ = q̃s[timestep-1]
        qₛ₋₁ = qs[timestep-1]
        q̇ₛ₋₁ = q̇s[timestep-1]
        pₛ₋₁ = ps[timestep-1]
        λₛ₋₁ = λs[timestep-1]
        μₛ₋₁ = μs[timestep-1]
        tₛ₋₁ = ts[timestep-1]
        q̃ₛ   = q̃s[timestep]
        qₛ   = qs[timestep]
        q̇ₛ   = q̇s[timestep]
        pₛ   = ps[timestep]
        λₛ   = λs[timestep]
        μₛ   = μs[timestep]
        qˣ = qₛ₋₁ .+ dt./2 .*q̇ₛ₋₁
        qₛ .= qˣ
        q̇ₛ .= q̇ₛ₋₁
        nu,active_indices,g = get_indices(qˣ)
        gₙ = g[active_indices]
        cs[timestep] = nu
        # ns_stepk! = ns_stepk_maker(nq,nλ,nμ,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,Aset,dyfuncs,invM,dt)
        # stepk! = stepk_maker(nq,nλ,nμ,nu,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,dt)
        isconverged = false
        res = typemax(eltype(qₛ))
        iteration_break = 0
        if nu == 0
            smooth_x[      1:nq]          .= qₛ
            smooth_x[   nq+1:nq+nλ]       .= 0.0
            smooth_x[nq+nλ+1:nq+nλ+nμ]    .= 0.0
            stepk! = stepk_maker(nq,nλ,nμ,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,dt)

            for iteration = 1:imax
                    stepk!(smooth_R,smooth_J,smooth_x,q̇ₛ)
                    res = norm(smooth_R)
                    if res < ftol
                        isconverged = true
                        iteration_break = iteration-1
                        break
                    end
                    smooth_Δx .= -smooth_J\smooth_R
                    smooth_x .+= smooth_Δx
                    qₛ .= smooth_x[      1:nq]
                    λₛ .= smooth_x[   nq+1:nq+nλ]
                    μₛ .= smooth_x[nq+nλ+1:nq+nλ+nμ]
                    pₛ .= Momentum_k(qₛ₋₁,pₛ₋₁,qₛ,λₛ,μₛ,M,A,B,dt)
                    q̇ₛ .= invM*pₛ
            end
            q̃ₛ .= smooth_x[      1:nq]
        else # u!=0
            nΛ = 3nu
            nonsmooth_nx = nq + nq + nλ + nμ + nΛ
            nonsmooth_x = zeros(eltype(q0),nonsmooth_nx)
            nonsmooth_R = zeros(eltype(q0),nonsmooth_nx)
            nonsmooth_J = zeros(eltype(q0),nonsmooth_nx,nonsmooth_nx)
            nonsmooth_x[            1:nq]             .= [0.2781866379712523, 0.0, 0.07282980498953609]
            nonsmooth_x[         nq+1:nq+nq]          .= [0.2774778420780419, 0.0, 0.07397747720342729]
            nonsmooth_x[nq+nq+nλ+nμ+1:nq+nq+nλ+nμ+nΛ] .= [(1.000000001)*3.8760483232954943, 0.0, 3.8760483232954943]
            y = copy(nonsmooth_x[nq+nq+nλ+nμ+1:nq+nq+nλ+nμ+nΛ])
            y .= [(1.000000001)*0.8591722059215059, 0.0, -0.8591721975917848]
            q̇ₛ .= [0.5142598661572809, 0.0, 1.400342517987585]
            isconverged = false
            # @show timestep, nu
            ip_ns_stepk! = ip_ns_stepk_maker(nq,nλ,nμ,nu,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,dt)
            μ = 1.0
            for iteration = 1:imax
                Dₛ,ηs,es,H = get_D(active_indices,qₛ)
                # ηs .= 1
                _,_,g = get_indices(qₛ)
                gₙ = g[active_indices]
                # @show iteration,Dₛ,ηs,es,gₙ
                μ,Δxc = ip_ns_stepk!(nonsmooth_R,nonsmooth_J,
                            nonsmooth_x,y,q̇ₛ,Dₛ,ηs,es,H,μ)
                res = norm(Δxc)
                @show iteration, res
                iteration_break = iteration
                if  res < ftol
                    isconverged = true
                    break
                end
                qₛ .= nonsmooth_x[         nq+1:nq+nq]
                λₛ .= nonsmooth_x[      nq+nq+1:nq+nq+nλ]
                μₛ .= nonsmooth_x[   nq+nq+nλ+1:nq+nq+nλ+nμ]
                pₛ .= Momentum_k(qₛ₋₁,pₛ₋₁,qₛ,λₛ,μₛ,M,A,B,dt)
                q̇ₛ .= invM*pₛ
            end
            # @show gₙ*9.81
            # @show 𝚲ₛ./dt
            q̃ₛ .= nonsmooth_x[      1:nq]
        end

        if !isconverged
            @warn "NLsolve max iterations $iteration_break, at timestep=$timestep, Res=$(res)"
            if exception
                error("Not converged!")
            else
                # intor.convergence = false
                # break
            end
        else
            if nu == 0
                stepstring = "Smooth"
                # @info "$stepstring timestep=$timestep, iterations=$iteration_break"
            else
                stepstring = "Nonsmooth with $nu contact(s)"
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
    ts,cs,qs,q̇s,ps,λs,μs
end

end
