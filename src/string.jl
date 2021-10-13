
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

function SHStringState(restlen, direction)
    SStringState(restlen,restlen,zero(restlen),1000*one(restlen),direction)
end

function SRStringState(restlen, direction)
    SStringState(restlen,restlen,zero(restlen),1628.58*one(restlen),direction)
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

struct SString{N,T}
    id::Int
    k::T
    c::T
    prestress::T
    state::SStringState{N,T}
end

function SString2D(id,origin_restlen::T,k::T,c=zero(k)) where T
    direction = MVector{2}(one(T),zero(T))
    state = SStringState(origin_restlen,direction)
    SString(id,k,c,state)
end

function SString3D(id,origin_restlen::T,k::T,c=zero(k)) where T
    direction = MVector{3}(one(T),zero(T),zero(T))
    state = SStringState(origin_restlen,direction)
    SString(id,k,c,state)
end

function SHString3D(id,origin_restlen::T,k::T,c=zero(k)) where T
    direction = MVector{3}(one(T),zero(T),zero(T))
    state = SHStringState(origin_restlen,direction)
    SString(id,k,c,state)
end

function SRString3D(id,origin_restlen::T,k::T,c=zero(k)) where T
    direction = MVector{3}(one(T),zero(T),zero(T))
    state = SRStringState(origin_restlen,direction)
    SString(id,k,c,state)
end

function SHRString3D(id,i,origin_restlen::T,k::T,prestress::T,c::T) where T
    direction = MVector{3}(one(T),zero(T),zero(T))
    if i%8 <= 4
        state = SHStringState(origin_restlen,direction)
    else
        state = SRStringState(origin_restlen,direction)
    end
    SString(id,k,c,prestress,state)
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

mutable struct ClusterSectionCablesState{N, T}
    restlen::T
    length::T 
    lengthdot::T
    tension::T
    direction::MArray{Tuple{N},T,1,N}
end

struct ClusterSectionCables{N, T}
    ID::Int
    state::ClusterSectionCablesState{N, T}
end

mutable struct ClusterCablesState{T}
    theta::Vector{T}
    s::Vector{T}
end

struct ClusterCables
    ID::Int
    state::ClusterCablesState
    section::Vector{ClusterSectionCables}
end

function ClusterSectionCables3D(id, restlen::T) where T
    direction = MVector{3}(one(T),zero(T),zero(T))
    section_state = ClusterSectionCablesState(restlen, restlen, 
        zero(restlen), zero(restlen), direction)
    return ClusterSectionCables(id, section_state)
end

function ClusterCables3D(Vector_Cluster_Section, section, id) 
    s = [Float64(0.0) for i in 1:section]
    theta = similar(s)
    state = ClusterCablesState(theta, s)
    return ClusterCables(id, state, Vector_Cluster_Section)
end