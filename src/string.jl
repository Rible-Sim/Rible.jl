
mutable struct LinearLaw{T}
    k::T
    F0::T
end

function (ll::LinearLaw)(Δl)
    @unpack F0, k = ll
    F = F0 + k*Δl
end

mutable struct SStringState{N,T}
    restlen::T
    length::T
    lengthdot::T
    tension::T
    direction::MArray{Tuple{N},T,1,N}
end

function SStringState(restlen,direction)
    SStringState(restlen,restlen,zero(restlen),zero(restlen),direction)
end

mutable struct SMAStringState{N,T}
    temp::T
    restlen::T
    length::T
    lengthdot::T
    tension::T
    direction::MArray{Tuple{N},T,1,N}
end

function SMAStringState(temp,restlen,direction)
    SMAStringState(temp,restlen,restlen,zero(restlen),zero(restlen),direction)
end

struct CableSegment{N,T}
    id::Int
    k::T
    c::T
    prestress::T
    original_restlen::T
    state::SStringState{N,T}
end

struct PrestressString{N,T}
    id::Int
    k::T
    c::T
    prestress::T
    state::SStringState{N,T}
end

struct SString{N,T}
    id::Int
    k::T
    c::T
    state::SStringState{N,T}
end

function SString2DState(original_restlen::T) where T
    direction = MVector{2}(one(T),zero(T))
    SStringState(original_restlen,direction)
end

function SString3DState(original_restlen::T) where T
    direction = MVector{3}(one(T),zero(T),zero(T))
    SStringState(original_restlen,direction)
end

function SString2D(id,original_restlen::T,k::T,c=zero(k)) where T
    state = SString2DState(original_restlen)
    SString(id,k,c,state)
end

function SString3D(id,original_restlen::T,k::T,c=zero(k)) where T
    state = SString3DState(original_restlen)
    SString(id,k,c,state)
end

function CableSegment2D(id,original_restlen::T,k::T;c=zero(k),prestress=zero(k)) where T
    state = SString2DState(original_restlen)
    CableSegment(id,k,c,prestress,original_restlen,state)
end

function CableSegment3D(id,original_restlen::T,k::T;c=zero(k),prestress=zero(k)) where T
    state = SString3DState(original_restlen)
    CableSegment(id,k,c,prestress,original_restlen,state)
end

struct SMAString{N,T,F}
    id::Int
    law::F
    c::T
    state::SMAStringState{N,T}
end

function SMAString3D(id,restlen::T,law::LawT,c::T,original_temp::T=0.0) where {LawT,T}
    direction = MVector{3}(one(T),zero(T),zero(T))
    state = SMAStringState(original_temp,restlen,direction)
    SMAString(id,law,c,state)
end

mutable struct SlidingPoint{T}
    μ::T
    θ::T
    α::T
    s::T
    s⁺::T
    s⁻::T
end

function calculate_α(μ,θ)
    exp(-μ*θ)
end

function SlidingPoint(μ)
    θ = one(μ) * 2pi
    #θ = zero(μ)
    α = calculate_α(μ,θ)
    s = zero(μ)
    s⁺,s⁻ = s2s̄(s)
    SlidingPoint(μ,θ,α,s,s⁺,s⁻)
end

struct ClusterCables{spsType,segsType}
    ID::Int
    sps::spsType
    segs::segsType
end

function ClusterCables3D(id, nsp, segs;μ=0.0)
    sps = StructArray([SlidingPoint(μ) for i = 1:nsp])
    ClusterCables(id, sps, segs)
end

function s2s̄(s::Number)
    abss = abs(s)
    s⁺ = (abss + s)/2
    s⁻ = (abss - s)/2
    s⁺,s⁻
end
