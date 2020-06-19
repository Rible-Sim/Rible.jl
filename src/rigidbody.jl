abstract type AbstractRigidBody end
abstract type TensegrityStructure end

struct RigidBody2DProperty{T}
    movable::Bool
    name::Symbol
    type::Symbol
    mass::T
    inertia::T
    CoM::SArray{Tuple{2},T,1,2}
    number_aps::Int
    anchorpoints::Vector{SArray{Tuple{2},T,1,2}}
end

function RigidBody2DProperty(i,movable,mass,inertia,CoM,aps)
    name = Symbol("rb"*string(i))
    type = :generic
    naps = length(aps)
    RigidBody2DProperty(movable,name,type,mass,inertia,CoM,naps,aps)
end


struct RigidBody2DState{T,CoordinatesType,AuxiliariesType}
    r::MArray{Tuple{2},T,1,2}
    θ::T
    ṙ::MArray{Tuple{2},T,1,2}
    ω::T
    p::Vector{MArray{Tuple{2},T,1,2}} # Anchor Points in global frame
    F::MArray{Tuple{2},T,1,2}
    τ::T
    Fanc::Vector{MArray{Tuple{2},T,1,2}}
    coords::CoordinatesType
    auxs::AuxiliariesType
end

struct RigidBody2DNaturalCoordinates{T}
    q::MArray{Tuple{4},T,1,4}
    q̇::MArray{Tuple{4},T,1,4}
    q̈::MArray{Tuple{4},T,1,4}
end

function RigidBody2DState(prop,ri,rj)
    q = MVector{4}(vcat(ri,rj))
    q̇ = zero(q)
    q̈ = zero(q)
    coords = RigidBody2DNaturalCoordinates(q,q̇,q̈)
    r = MVector{2}(ri)
    ṙ = zero(r)
    rirj = rj-ri
    θ = atan(rirj[2],rirj[1])
    ω = 0.0
    F = MVector(0.0,0.0)
    τ = 0.0
    τanc = MVector(0.0,0.0)
    aux = NCaux(prop,ri,rj)
    nap = prop.number_aps
    p = [MVector{2}(aux.Cp[i]*q) for i in 1:nap]
    Fanc = [zeros(MVector{2}) for i in 1:nap]
    RigidBody2DState(r,θ,ṙ,ω,p,F,τ,Fanc,coords,aux)
end

struct RigidBody2D{T,CoordinatesType,AuxiliariesType} <: AbstractRigidBody
    prop::RigidBody2DProperty{T}
    state::RigidBody2DState{T,CoordinatesType,AuxiliariesType}
end

function reset_forces!(rb::RigidBody2D)
    for f in rb.state.Fanc
        f .= 0.0
    end
    rb.state.F .= 0.0
end
