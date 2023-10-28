module URDF
using LinearAlgebra
using StaticArrays
using ColorTypes

struct Origin{T}
    xyz::SVector{3,T}
    rpy::SVector{3,T}
end

struct Inertial{T}
    origin::Origin{T}
    mass::T
    inertia::SMatrix{3,3,T}
end

struct Visual{T,gT,mT}
    name::Union{Nothing,String}
    origin::Origin{T}
    geometry::gT
    material::mT
end

abstract type Geometry{T} end
struct Box{T} <: Geometry{T}
    size::SVector{3,T}
end

struct Cylinder{T} <: Geometry{T}
    radius::T
    length::T
end

struct Sphere{T} <: Geometry{T}
    radius::T
end

struct MeshGeometry{T} <: Geometry{T}
    filename::String
    scale::SVector{3,T}
end

struct Material{colorType,textureType}
    name::String
    color::colorType
    texture::textureType
end

struct Collision{T,gT}
    name::Union{Nothing,String}
    origin::Origin{T}
    geometry::gT
end

mutable struct Link{T}
    name::String
    inertial::Union{Nothing,Inertial{T}}
    visual::Union{Nothing,Visual{T}}
    collision::Union{Nothing,Collision{T}}
end


struct Mobility{T,limitType}    
    axis::SVector{3,T}
    calibration::Union{Nothing,NamedTuple{(:rising, :falling), Tuple{T, T}}}
    dynamics::NamedTuple{(:damping, :friction), Tuple{T, T}}
    limit::Union{Nothing,limitType}
    # mimic::Mimic{T}
    # safety_controller::Safety{T}
end

mutable struct Joint{T}
    name::String
    type::String
    origin::Origin{T}
    parent::String
    child::String
    mobility::Union{Nothing,Mobility{T}}
end


end