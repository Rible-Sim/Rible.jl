
# Discrete point-set surfaces (e.g. for defining specific contact points)
function StaticContactSurfaces(surfaces, field=Gravity(9.81))
    StaticEnvironment(
        surfaces,
        field
    )
end

abstract type ContactGeometry <: AbstractGeometry end 
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

struct Sphere{T, N} <: ConvexContactPrimitive
    radius::T
    r::SArray{Tuple{N},T,1,N}
end

function Sphere(radius::T, r::AbstractVector{T}) where T
    Sphere(radius, SVector{3}(r))
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

struct Polytope{T,N,M,L} <: ConvexContactPrimitive
    vertices::SMatrix{N,M,T,L}
end

function signed_distance(x::StaticArray{Tuple{N}, T, 1},p::Union{Plane,HalfSpace}) where {T,N}
    (;n, d) = p
    transpose(n)*x + d
end

function signed_distance(x::StaticArray{Tuple{N}, T, 1},p::Sphere) where {T,N}
    norm(x-p.r) - p.radius
end

function distance(x::StaticArray{Tuple{N}, T, 1},p::Plane) where {T,N}
    abs(signed_distance(x,p))
end

function ison(x::StaticArray{Tuple{N}, T, 1},p::Plane;tol=eps(T)) where {T,N}
    distance(x,p) < tol
end

function contact_gap_and_normal(x::StaticArray{Tuple{N}, T, 1},p::Sphere) where {T,N}
    (; r, radius) = p
    gap = norm(x-r) - radius
    normal = (x-r)/norm(x-r)
    gap, normal
end

function contact_gap_and_normal(x::StaticArray{Tuple{N}, T, 1},spheres::Vector{<:Sphere{T,N}}) where {T,N}
    nhs = length(spheres)
    gaps = zeros(T,nhs)
    normals = Vector{SVector{N,T}}(undef,nhs)
    for (i,hs) in enumerate(spheres)
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

function contact_gap_and_normal(x::StaticArray{Tuple{N}, T, 1},p::Union{Plane,HalfSpace}) where {T,N}
    signed_distance(x,p), p.n
end

function contact_gap_and_normal(x::StaticArray{Tuple{N}, T, 1},cp::Vector{<:Plane}) where {T,N}
    gap_first, n_first = contact_gap_and_normal(x,first(cp))
    for p in cp[begin+1:end]
        gap, n = contact_gap_and_normal(x,p)
        if gap < 0
            return gap, n 
        end
    end
    gap_first, n_first
end

function contact_gap_and_normal(x::StaticArray{Tuple{N}, T, 1},nogeometry::NoGeometry) where {T,N}
    gap = typemax(T)
    normal = @SVector zeros(N)
    gap, normal
end

function contact_gap_and_normal(x::StaticArray{Tuple{N}, T, 1},halvspaces::Vector{<:HalfSpace{T,N}}) where {T,N}
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

