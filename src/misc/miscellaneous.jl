struct Plane{T}
    n::SArray{Tuple{3},T,1,3}
    d::T
    r::SArray{Tuple{3},T,1,3}
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

function signed_distance(x::AbstractVector{T},p::Plane) where T
    (;n, d) = p
    transpose(n)*x + d
end

function distance(x::AbstractVector{T},p::Plane) where T
    abs(signed_distance(x,p))
end

function ison(x::AbstractVector{T},p::Plane;tol=eps(T)) where T
    distance(x,p) < tol
end
