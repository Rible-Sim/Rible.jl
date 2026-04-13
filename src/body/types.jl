"""
Rigid Body Types

This module defines the concrete types for rigid bodies.
"""

"""
Rigid Body Property Type 
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct RigidBodyProperty{N,T,L} <: AbstractRigidBodyProperty{N,T}
    "Is able to make contact with?"
    contactable::Bool
    "Is visible?"
    visible::Bool
    "id. Unique in a system"
    id::Int
    "Type or name."
    type::Symbol
    "Mass"
    mass::T
    "Inertia"
    inertia::SMatrix{N,N,T,L}
    "Centroid in local frame"
    mass_locus::Locus{N,T,L}
    "Anchor points in local frame"
    loci::Vector{Locus{N,T,L}}
end

"""
Rigid Body State Mutable Type.
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
mutable struct RigidBodyState{N,M,T,L} <: AbstractRigidBodyState{N,T}
    "Origin of local frame"
    origin_frame::CartesianFrame{N,M,T,L}
    "Position of mass center in global frame"
    mass_locus_state::LocusState{N,M,T,L}
    "Positions of anchor points in global frame"
    loci_states::Vector{LocusState{N,M,T,L}}
end

"""
Rigid Body Cache Type.
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct RigidBodyCache{CType,cacheType,jacType} <: AbstractBodyCache
    Co::CType
    Cg::CType
    Cps::Vector{CType}
    inertia_cache::cacheType
    jac_cache::jacType
end

"""
Rigid Body Type 
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct RigidBody{N,M,T,L,coordsType,cacheType,meshType} <: AbstractRigidBody{N,T}
    "Rigid Body Property "
    prop::RigidBodyProperty{N,T,L}
    "Rigid Body State "
    state::RigidBodyState{N,M,T,L}
    "Coordinates State"
    coords::coordsType
    "Cache"
    cache::cacheType
    "Rigid Body Mesh"
    mesh::meshType
end

struct FlexibleBodyProperty{N,T,L} <: AbstractFlexibleBodyProperty{N,T}
    "Is able to make contact with?"
    contactable::Bool
    "Is visible?"
    visible::Bool
    id::Int
    type::Symbol
    mass::T
    mass_locus::Locus{N,T,L}
    loci::Vector{Locus{N,T,L}}
end

struct FlexibleBodyCoordinatesCache{fT,eT,inertia_cacheType,ArrayT}
    funcs::fT
    e::eT
    ė::eT
    inertia_cache::inertia_cacheType
    So::ArrayT
    Sg::ArrayT
    Sps::Vector{ArrayT}
end

struct FlexibleBodyState{N,M,T,L} <: AbstractFlexibleBodyState{N,T}
    mass_locus_state::LocusState{N,M,T,L}
	"Positions of anchor points"
    loci_states::Vector{LocusState{N,M,T,L}}
end

struct FlexibleBody{
        N,T,
        propType <: AbstractFlexibleBodyProperty{N,T},
        stateType <: AbstractFlexibleBodyState{N,T},
        coordsType <: AbstractCoordinates{N,T},
        cacheType <: AbstractBodyCache,
        meshType
    } <: AbstractFlexibleBody{N,T}
    "Properties of a flexible body"
    prop::propType
    "State of a flexible body"
    state::stateType
    "Coordinates Info"
    coords::coordsType
    "Cache"
    cache::cacheType
    "Mesh"
    mesh::meshType
end


# Constructor implementations
"""
Rigid Body Property Constructor 
$(TYPEDSIGNATURES)
"""
function RigidBodyProperty(
        id::Integer,contactable::Bool,
        mass::T,inertia_input,
        mass_locus::Locus{N,T,L},
        loci::Vector{<:Locus} = Locus{N,T,L}[];
        visible = true,
        type = :generic
    ) where {T,N,L}
    mtype = StaticArrays.similar_type(inertia_input)
    inertia = mtype(inertia_input)
    loci_converted = Vector{Locus{N,T,L}}(loci)
    
    @debug "Logging RigidBodyProperty" id visible contactable mass inertia mass_locus loci_converted
    return RigidBodyProperty(
        contactable,
        visible,
        id,type,
        mass,
        inertia,
        mass_locus,
        loci_converted
    )
end

"""
Rigid Body State Constructor 
$(TYPEDSIGNATURES)
---
"""
function RigidBodyState(prop::RigidBodyProperty{N,T},
        origin_position_input,
        rotation_input,
        origin_velocity_input,
        angular_velocity_input
    ) where {N,T}
    (;mass_locus,loci) = prop
    num_of_loci = length(loci)
    origin_position = MVector{N}(origin_position_input)
    origin_velocity = MVector{N}(origin_velocity_input)
    if rotation_input isa Number
        rotation = rotation_matrix(rotation_input)
    else
        rotation = Matrix(rotation_input)
    end
    axes = Axes(SMatrix{N,N,T}(rotation))
    M = 2N-3
    ω = MVector{M}(angular_velocity_input)
    origin_frame = CartesianFrame(
        origin_position,
        origin_velocity,
        axes,
        ω
    )
    mass_locus_state = LocusState(
        mass_locus,
        origin_frame,
    )
    loci_states = [
        LocusState(
            lo,
            origin_frame,
        )
        for lo in loci
    ]
    RigidBodyState(
        origin_frame,
        mass_locus_state,
        loci_states,
    )
end

function RigidBody(prop,state,coords::AbstractCoordinates,mesh=nothing)
    cache = RigidBodyCache(prop,coords)
    RigidBody(prop,state,coords,cache,mesh)
end
