"""
Coordinates State Type.
$(TYPEDEF)
"""
mutable struct CoordinatesState{T,qT,sT,pT,cT} <: AbstractCoordinatesState
    t::T
    q::qT
    q̇::qT
    q̈::qT
    F::qT
    λ::qT
    s::sT
    p::pT
    c::cT
end

get_coords(coords::CoordinatesState) = coords.q
get_params(coords::CoordinatesState) = coords.c
get_auxiliary(coords::CoordinatesState) = coords.s
get_velocs(coords::CoordinatesState) = coords.q̇
get_accels(coords::CoordinatesState) = coords.q̈
get_free_coords(coords::CoordinatesState) = coords.q
get_free_velocs(coords::CoordinatesState) = coords.q̇
get_free_accels(coords::CoordinatesState) = coords.q̈
get_pres_coords(coords::CoordinatesState{T}) where T = T[]
get_pres_velocs(coords::CoordinatesState{T}) where T = T[]
get_pres_accels(coords::CoordinatesState{T}) where T = T[]
get_multipliers(coords::CoordinatesState) = coords.λ