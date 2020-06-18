
mutable struct SStringState{T}
    restlength::T
    length::T
    lengthdot::T
    tension::T
    direction::MArray{Tuple{2},T,1,2}
end

function SStringState(restlength,length,tension)
    SStringState(restlength,length,zero(length),tension,MVector(1.0,0.0))
end

struct SString{T}
    k::T
    c::T
    original_restlength::T
    state::SStringState{T}
end

function SString(k,origin_restlength,state)
    SString(k,zero(k),origin_restlength,state)
end
