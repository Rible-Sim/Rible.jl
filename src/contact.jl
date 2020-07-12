struct ContactFrame{T}
    n::MArray{Tuple{3},T,1,3}
    t1::MArray{Tuple{3},T,1,3}
    t2::MArray{Tuple{3},T,1,3}
end

function ContactFrame{T}() where {T}
    i = one(T)
    o = zero(T)
    n = MVector(o,o,i)
    t1 = MVector(i,o,o)
    t2 = MVector(o,i,o)
    ContactFrame(n,t1,t2)
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
    frame::ContactFrame{T}
    m::Int
    e::Vector{T}
    d::Vector{SArray{Tuple{3},T,1,3}}
    D::Array{T,2}
    impulse::ApproxFrictionalImpulse{T}
end

function ApproxFrictionalContact(ϵ::T,μ::T,m::Int) where T
    frame = ContactFrame{T}()
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
