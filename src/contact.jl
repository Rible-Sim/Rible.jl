module CollisionDetection
using LinearAlgebra
using StaticArrays
abstract type CollisionShape end
abstract type CollisionPrimitive <: CollisionShape end
struct Sphere{T} <: CollisionPrimitive
    radius::T
end
struct Plane{T} <: CollisionPrimitive
    normal::SArray{Tuple{3},T,1,3}
    d::T
end
function Plane(normal,d)
    norm_normal = normal ./ norm(normal)
    Plane(SVector{3}(norm_normal),d)
end

struct Halfspace{T} <: CollisionPrimitive
    normal::SArray{Tuple{3},T,1,3}
    d::T
end
function Halfspace(normal,d)
    norm_normal = normal ./ norm(normal)
    Halfspace(SVector{3}(norm_normal),d)
end

struct CollisionObject{A,T}
    shape::A
    translation::MArray{Tuple{3},T,1,3}
    rotation::MArray{Tuple{3,3},T,2,9}
end
function sphere_object(radius,center)
    shape = Sphere(radius)
    translation = MVector{3}(center)
    rotation = MMatrix{3,3}(Matrix(one(radius)*I,3,3))
    CollisionObject(shape,translation,rotation)
end
struct CapsuleObject{T}
    c1::MArray{Tuple{3},T,1,3}
    c2::MArray{Tuple{3},T,1,3}
    radius::T
end
function halfspace_object(normal,d)
    shape = Halfspace(normal,d)
    translation = @MVector zeros(typeof(d),3)
    rotation = MMatrix{3,3}(Matrix(one(d)*I,3,3))
    CollisionObject(shape,translation,rotation)
end
function plane_object(normal,d)
    shape = Plane(normal,d)
    translation = @MVector zeros(typeof(d),3)
    rotation = MMatrix{3,3}(Matrix(one(d)*I,3,3))
    CollisionObject(shape,translation,rotation)
end
mutable struct Contact{T}
    status::Symbol
    position::MArray{Tuple{3},T,1,3}
    penetration::T
    normal::MArray{Tuple{3},T,1,3}
    contact_points::Array{MArray{Tuple{3},T,1,3},1}
end
# closest_point??
function isintersect(sphere::CollisionObject{Sphere{T},T},plane::CollisionObject{Plane{T},T}) where T
    dist = sphere.translation ⋅ plane.shape.normal - plane.shape.d
    return abs(dist) <= sphere.shape.radius
end
function isintersect(sphere::CollisionObject{Sphere{T},T},hs::CollisionObject{Halfspace{T},T}) where T
    dist = sphere.translation ⋅ hs.shape.normal - hs.shape.d
    return dist <= sphere.shape.radius
end
function isinside(sphere::CollisionObject{Sphere{T},T},hs::CollisionObject{Halfspace{T},T}) where T
    dist = sphere.translation ⋅ hs.shape.normal - hs.shape.d
    return dist < -sphere.shape.radius
end
function isintersect(s1::CollisionObject{Sphere,T},s2::CollisionObject{Sphere,T}) where T
    vec12 = s2.translation - s1.translation
    return norm(vec12) <= s1.shape.radius + s2.shape.radius
end
const ZeroMVector3Float64 = MVector(0.0,0.0,0.0)
const emptycontact = Contact{Float64}(:nocontact,
                                    ZeroMVector3Float64,
                                    0.0,
                                    ZeroMVector3Float64,
                                    [ZeroMVector3Float64])
function contact(s1::CollisionObject{Sphere{T},T},s2::CollisionObject{Sphere{T},T}) where T
    if s2.translation ≈ s1.translation
        return emptycontact
    else
        vec12 = s2.translation - s1.translation
        normvec12 = norm(vec12)
        normal = vec12./normvec12
        penetration =  normvec12 - s1.shape.radius - s2.shape.radius
        if penetration > 0.0
            status = :nocontact
        else
            status = :contact
        end
        contact_position = s1.translation + normal*(s1.shape.radius + penetration/2)
        return Contact(status,
                        contact_position,
                        penetration,
                        normal,
                        [contact_position])
    end
end
function contact(hs::CollisionObject{Halfspace,T},sphere::CollisionObject{Sphere,T}) where T
    dist = sphere.translation ⋅ hs.shape.normal - hs.shape.d
    penetration = dist - sphere.radius
    if penetration > 0.0
        status = :nocontact
    else
        status = :contact
    end
    normal = hs.shape.normal
    contact_position = sphere.translation - normal * (sphere.shape.radius + penetration/2)
    return Contact(status,
                    contact_position,
                    penetration,
                    normal,
                    [contact_position])
end

end  # module CollisionDetection
