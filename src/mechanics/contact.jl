struct Contact{T}
    id::Int
    μ::T
    e::T
    state::FrictionalContactState{3,T}
end

function Contact(id,μ,e)
    active = false
    persistent = true
    i = one(μ)
    o = zero(μ)
    gap = i
    n = SVector(i,o,o)
    frame = spatial_frame(n)
    v = SVector(o,o,o)
    Λ = SVector(o,o,o)
    state = FrictionalContactState(active,persistent,gap,frame,v,Λ)
    Contact(id,μ,e,state)
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
    frame::Axes{3,T}
    m::Int
    e::Vector{T}
    d::Vector{SArray{Tuple{3},T,1,3}}
    D::Array{T,2}
    impulse::ApproxFrictionalImpulse{T}
end

function ApproxFrictionalContact(ϵ::T,μ::T,m::Int) where T
    frame = spatial_frame{T}()
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
