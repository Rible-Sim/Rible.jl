
module NSSFC
import ..TensegrityRobots as TR
using Parameters
using LinearAlgebra
using StaticArrays
using BlockDiagonals

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

function APGD(r,μs,𝐍,𝐫;τ=1e-5,Nmax=20)
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
            @error "APGD: Max iteration $k, res=$(rₖ₊₁)"
        end
        # update
        𝐲ₖ .= 𝐲ₖ₊₁
        𝛄ₖ .= 𝛄ₖ₊₁
        Lₖ = 0.9Lₖ
        tₖ = 1/Lₖ
    end
    return 𝛄ₕ,rmin
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

end
