struct RigidBody3DProperty{T}
    movable::Bool
    name::Symbol
    type::Symbol
    mass::T
    inertia::SArray{Tuple{3,3},T,2,9}
    CoM::SArray{Tuple{3},T,1,3}
    naps::Int
    aps::Vector{SArray{Tuple{3},T,1,3}} # Anchor Points in local frame
end

function RigidBody3DProperty(i,movable,mass,inertia,CoM,aps)
    name = Symbol("rb"*string(i))
    type = :generic
    naps = length(aps)
    RigidBody3DProperty(movable,name,type,mass,inertia,CoM,naps,aps)
end

struct RigidBody3DState{T,CoordinatesType,cacheType}
    r::MArray{Tuple{3},T,1,3}
    R::MArray{Tuple{3,3},T,2,9}
    ṙ::MArray{Tuple{3},T,1,3}
    ω::MArray{Tuple{3},T,1,3}
    p::Vector{MArray{Tuple{3},T,1,3}} # Anchor Points in global frame
    τ::MArray{Tuple{3},T,1,3}
    Faps::Vector{MArray{Tuple{3},T,1,3}}
    τaps::Vector{MArray{Tuple{3},T,1,3}}
    coords::CoordinatesType
    cache::cacheType
end

struct RigidBodyCoordinates{T,N}
    q::MArray{Tuple{N},T,1,N}
    q̇::MArray{Tuple{N},T,1,N}
    q̈::MArray{Tuple{N},T,1,N}
    Q::MArray{Tuple{N},T,1,N} # Generalized force
end


struct RigidBody3D{T,CoordinatesType,cacheType} <: AbstractRigidBody
    prop::RigidBody3DProperty{T}
    state::RigidBody3DState{T,CoordinatesType,cacheType}
end

struct NaturalCoordinatesCache{ArrayT,MT,fT}
    CG::ArrayT
    Cp::Vector{ArrayT}
    M::MT
    funcs::fT
end

function NaturalCoordinatesAuxiliaries(M,CG,Cp)
    StaticM =SMatrix{12,12}(M)
    StaticCG = SMatrix{3,12}(CG)
    numberofpoints = length(Cp)
    StaticCp = SArray{Tuple{numberofpoints},SArray{Tuple{3,12},eltype(M),2,36},1,numberofpoints}(Cp)
    NaturalCoordinatesAuxiliaries(StaticM,StaticCG,StaticCp)
end

function quat_multiply(q0,q1)
    w0,x0,y0,z0 = q0
    w1,x1,y1,z1 = q1
    [w0*w1 - x0*x1 - y0*y1 - z0*z1,
     w0*x1 + w1*x0 + y0*z1 - z0*y1,
     w0*y1 + w1*y0 + z0*x1 - x0*z1,
     w0*z1 + w1*z0 + x0*y1 - y0*x1]
end
function tilde(x)
    ret = Matrix{eltype(x)}(undef,3,3)
    ret[1,1] = 0.0
    ret[1,2] = -x[3]
    ret[1,3] = x[2]
    ret[2,1] = x[3]
    ret[2,2] = 0.0
    ret[2,3] = -x[1]
    ret[3,1] = -x[2]
    ret[3,2] = x[1]
    ret[3,3] = 0.0
    ret
end
