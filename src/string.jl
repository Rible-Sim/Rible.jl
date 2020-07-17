
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

struct SString{N,T}
    k::T
    c::T
    original_restlen::T
    state::SStringState{N,T}
end

function SString2D(origin_restlen::T,k::T,c=zero(k)) where T
    direction = MVector{2}(one(T),zero(T))
    state = SStringState(origin_restlen,direction)
    SString(k,c,origin_restlen,state)
end

function SString3D(origin_restlen::T,k::T,c=zero(k)) where T
    direction = MVector{3}(one(T),zero(T),zero(T))
    state = SStringState(origin_restlen,direction)
    SString(k,c,origin_restlen,state)
end

struct Actuator{N,NS,T}
    strings::SArray{Tuple{NS},SString{N,T},1,NS}
end
