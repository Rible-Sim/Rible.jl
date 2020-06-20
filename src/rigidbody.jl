abstract type AbstractRigidBody end
abstract type AbstractTensegrityStructure end

struct RigidBodyProperty{N,T,iT}
    movable::Bool
    name::Symbol
    type::Symbol
    mass::T
    inertia::iT
    CoM::SArray{Tuple{N},T,1,N}
    naps::Int
    aps::Vector{SArray{Tuple{N},T,1,N}}
end

function RigidBodyProperty(i,movable,mass,inertia,CoM,aps)
    name = Symbol("rb"*string(i))
    type = :generic
    naps = length(aps)
    RigidBodyProperty(movable,name,type,mass,inertia,CoM,naps,aps)
end


struct RigidBodyState{N,T,L,ωT,CoordinatesType,cacheType}
    r::MArray{Tuple{N},T,1,N}
    R::MArray{Tuple{N,N},T,2,L}
    ṙ::MArray{Tuple{N},T,1,N}
    ω::ωT
    p::Vector{MArray{Tuple{N},T,1,N}} # Anchor Points in global frame
    F::MArray{Tuple{N},T,1,N}
    τ::ωT
    Fanc::Vector{MArray{Tuple{N},T,1,N}}
    τanc::Vector{MArray{Tuple{N},T,1,N}}
    coords::CoordinatesType
    cache::cacheType
end

struct RigidBodyCoordinates{N,T}
    q::MArray{Tuple{N},T,1,N}
    q̇::MArray{Tuple{N},T,1,N}
    q̈::MArray{Tuple{N},T,1,N}
    Q::MArray{Tuple{N},T,1,N}
end

function RigidBodyState(prop,ri,rj)
    q = MVector{4}(vcat(ri,rj))
    q̇ = zero(q)
    q̈ = zero(q)
    Q = zero(q)
    coords = RigidBodyCoordinates(q,q̇,q̈,Q)
    r = MVector{2}(ri)
    ṙ = zero(r)
    u = rj-ri
    θ = atan(u[2],u[1])
    R = MMatrix{2,2}(cos(θ), -sin(θ), sin(θ), cos(θ))
    ω = 0.0
    F = MVector(0.0,0.0)
    τ = 0.0
    τanc = MVector(0.0,0.0)
    cache = NaturalCoordinatesCache(prop,norm(u))
    nap = prop.naps
    p = [MVector{2}(cache.Cp[i]*q) for i in 1:nap]
    Fanc = [zeros(MVector{2}) for i in 1:nap]
    τanc = [zeros(MVector{2}) for i in 1:nap]
    RigidBodyState(r,R,ṙ,ω,p,F,τ,Fanc,τanc,coords,cache)
end

struct RigidBody{N,T,iT,L,ωT,CoordinatesType,cacheType} <: AbstractRigidBody
    prop::RigidBodyProperty{N,T,iT}
    state::RigidBodyState{N,T,L,ωT,CoordinatesType,cacheType}
end

function reset_forces!(rb::RigidBody)
    for f in rb.state.Fanc
        f .= 0.0
    end
    rb.state.F .= 0.0
end
