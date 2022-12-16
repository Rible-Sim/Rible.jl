struct SpatialFrame{T}
    n::SArray{Tuple{3},T,1,3}
    t1::SArray{Tuple{3},T,1,3}
    t2::SArray{Tuple{3},T,1,3}
end

function SpatialFrame(n)
    n /= norm(n)
    t1,t2 = NaturalCoordinates.HouseholderOrthogonalization(n)
    SpatialFrame(SVector{3}(n),SVector{3}(t1),SVector{3}(t2))
end

mutable struct FrictionalContactState{T}
    active::Bool
    persistent::Bool
    gap::T
    frame::SpatialFrame{T}
    v::SArray{Tuple{3},T,1,3}
    Λ::SArray{Tuple{3},T,1,3}
end

struct Contact{T}
    id::Int
    μ::T
    e::T
    state::FrictionalContactState{T}
end

function Contact(id,μ,e)
    active = false
    persistent = true
    i = one(μ)
    o = zero(μ)
    gap = i
    n = SVector(i,o,o)
    frame = SpatialFrame(n)
    v = SVector(o,o,o)
    Λ = SVector(o,o,o)
    state = FrictionalContactState(active,persistent,gap,frame,v,Λ)
    Contact(id,μ,e,state)
end

activate!(c::Contact,gap) = activate!(c.state,gap)

function activate!(state::FrictionalContactState,gap)
    state.gap = gap
    pre_active = state.active
    state.active = ifelse(gap<=0, true, false)
    state.persistent = (pre_active == state.active)
end

function update_contacts!(contacts, v, Λ)
    vs = split_by_lengths(v,3)
    Λs = split_by_lengths(Λ,3)
    for (ac,va,Λa) in zip(contacts,vs,Λs)
        ac.state.v = SVector{3}(va)
        ac.state.Λ = SVector{3}(Λa)
    end
end

mutable struct ApproxFrictionalImpulse{T}
    active::Bool
    Λn::T
    β::Vector{T}
end

ApproxFrictionalImpulse{T}(m) where T = ApproxFrictionalImpulse(false,zero(T),zeros(T,m))

struct ApproxFrictionalContact{T}
    ϵ::T
    μ::T
    frame::SpatialFrame{T}
    m::Int
    e::Vector{T}
    d::Vector{SArray{Tuple{3},T,1,3}}
    D::Array{T,2}
    impulse::ApproxFrictionalImpulse{T}
end

function ApproxFrictionalContact(ϵ::T,μ::T,m::Int) where T
    frame = SpatialFrame{T}()
    e = ones(m)
    d = [SVector{3,T}(cos(2π*i/m),sin(2π*i/m),0.0) for i = 0:m-1]
    D = Matrix{T}(undef,3,m)
    for i = 1:m
        D[:,i] .= d[i]
    end
    impulse = ApproxFrictionalImpulse{T}(m)
    ApproxFrictionalContact(
    ϵ,μ,frame,m,e,d,D,impulse
    )
end
