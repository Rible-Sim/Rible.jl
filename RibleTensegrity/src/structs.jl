mutable struct LinearLaw{T}
    k::T
    F0::T
end

function (ll::LinearLaw)(Δl)
    (;F0, k) = ll
    F = F0 + k*Δl
end

"""
$(TYPEDSIGNATURES)
"""
mutable struct SMADistanceSpringDamperState{N,T}
    temp::T
    restlen::T
    length::T
    lengthdot::T
    tension::T
    direction::MArray{Tuple{N},T,1,N}
end

function SMADistanceSpringDamperState(temp,restlen,direction)
    SMADistanceSpringDamperState(temp,restlen,restlen,zero(restlen),zero(restlen),direction)
end


struct SMADistanceSpringDamper{N,T,F} <: AbstractForce
    id::Int
    law::F
    c::T
    state::SMADistanceSpringDamperState{N,T}
end

function SMADistanceSpringDamper3D(restlen::T,law::LawT,c::T,original_temp::T=0.0) where {LawT,T}
    direction = MVector{3}(one(T),zero(T),zero(T))
    state = SMADistanceSpringDamperState(original_temp,restlen,direction)
    SMADistanceSpringDamper(law,c,state)
end
