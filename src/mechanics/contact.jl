struct Contact{T}
    id::Int
    μ::T
    e::T
    state::ContactState{3,T}
end

function Contact(id,μ,e)
    active = false
    persistent = true
    i = one(μ)
    o = zero(μ)
    gap = i
    n = SVector(i,o,o)
    frame = spatial_axes(n)
    v = SVector(o,o,o)
    Λ = SVector(o,o,o)
    state = ContactState(active,persistent,gap,frame,v,Λ)
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
    frame = spatial_axes{T}()
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

abstract type AbstractContactEnvironment end
abstract type StaticContactEnvironment <: AbstractContactEnvironment end
struct StaticContactSurfaces{surfacesType} <: StaticContactEnvironment 
    surfaces::surfacesType
end

struct EmptySpace <: StaticContactEnvironment end
struct GravitySpace{T} <: StaticContactEnvironment 
    gravity::T
end

abstract type ContactGeometry end
abstract type ContactPrimitive <: ContactGeometry end
abstract type ConvexContactPrimitive <: ContactPrimitive end

struct Plane{T,N} <: ConvexContactPrimitive
    n::SArray{Tuple{N},T,1,N}
    d::T
    r::SArray{Tuple{N},T,1,N}
end

function Plane(n::AbstractVector{T},r::AbstractVector{T}) where T
    n /= norm(n)
    a = n[1]
    b = n[2]
    c = n[3]
    d = -(a*r[1]+b*r[2]+c*r[3])
    Plane(SVector{3}(n),d,SVector{3}(r))
end

function Plane(n::AbstractVector{T},d::T) where T
    n /= norm(n)
    a = n[1]
    b = n[2]
    c = n[3]
    o = zero(d)
    x = o
    y = o
    z = d/(-c)
    Plane(SVector{3}(n),d,SVector{3}(x,y,z))
end

function Plane(a::T,b::T,c::T,d::T) where T
    n = [a,b,c]
    n /= norm(n)
    Plane(n,d)
end

struct Sphere{T} <: ConvexContactPrimitive
    radius::T
end

struct HalfSpace{T,N} <: ConvexContactPrimitive
    n::SArray{Tuple{N},T,1,N}
    d::T
    r::SArray{Tuple{N},T,1,N}
end

function HalfSpace(n::AbstractVector{T},r::AbstractVector{T}) where T
    n /= norm(n)
    a = n[1]
    b = n[2]
    c = n[3]
    d = -(a*r[1]+b*r[2]+c*r[3])
    HalfSpace(SVector{3}(n),d,SVector{3}(r))
end

function HalfSpace(n::AbstractVector{T},d::T) where T
    n /= norm(n)
    a = n[1]
    b = n[2]
    c = n[3]
    o = zero(d)
    x = o
    y = o
    z = d/(-c)
    HalfSpace(SVector{3}(n),d,SVector{3}(x,y,z))
end

function HalfSpace(a::T,b::T,c::T,d::T) where T
    n = [a,b,c]
    n /= norm(n)
    HalfSpace(n,d)
end

struct Polytope{T,N,M} <: ConvexContactPrimitive
    vertices::SMatrix{N,M,T}
end

function signed_distance(x::AbstractVector{T},p::Union{Plane,HalfSpace}) where T
    (;n, d) = p
    transpose(n)*x + d
end

function distance(x::AbstractVector{T},p::Plane) where T
    abs(signed_distance(x,p))
end

function ison(x::AbstractVector{T},p::Plane;tol=eps(T)) where T
    distance(x,p) < tol
end

function contact_gap_and_normal(x::AbstractVector,p::Union{Plane,HalfSpace})
    signed_distance(x,p), p.n
end

function contact_gap_and_normal(x::AbstractVector,cp::Vector{<:Plane})
    gap_first, n_first = contact_gap_and_normal(x,first(cp))
    for p in cp[begin+1:end]
        gap, n = contact_gap_and_normal(x,p)
        if gap < 0
            return gap, n 
        end
    end
    gap_first, n_first
end

function contact_gap_and_normal(x::AbstractVector,halvspaces::Vector{<:HalfSpace{T,N}}) where {T,N}
    nhs = length(halvspaces)
    gaps = zeros(T,nhs)
    normals = Vector{SVector{N,T}}(undef,nhs)
    for (i,hs) in enumerate(halvspaces)
        gap, n = contact_gap_and_normal(x,hs)
        gaps[i] = gap
        normals[i] = n
    end
    if all(gaps .< 0) # penetration
        ihs = argmax(gaps)
        return gaps[ihs], normals[ihs]
    else
        idx = findall((x)->x>=0,gaps)
        postive_gaps = @view gaps[idx]
        postive_normals = @view normals[idx]
        ihs = argmin(postive_gaps)
        return postive_gaps[ihs], postive_normals[ihs]
    end
end
