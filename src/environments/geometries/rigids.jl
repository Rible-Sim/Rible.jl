
struct ContactRigidBodies{bodiesType,fieldType} <: RigidBodyContactEnvironment
    contact_bodies::bodiesType
    field::fieldType
end

function ContactRigidBodies(contact_bodies;gravity=Gravity())
    ContactRigidBodies(
        contact_bodies,
        gravity
    )
end

abstract type ContactRigid end 

mutable struct RigidSphere{T} <: ContactRigid
    bid::Int
    pid::Int
    radius::T
end

mutable struct RigidCylinder{T} <: ContactRigid
    bid::Int
    pid::Int
    radius::T
    halflen::T
end

mutable struct RigidCube{T} <: ContactRigid
    bid::Int
    pid::Int
    length::T
    width::T
    height::T
end

function contact_gap_and_normal(x::StaticArray{Tuple{N}, T, 1},contact::RigidSphere,bodies) where {T,N}
    (;bid) = contact_body
    foreach(bodies) do body
        if body.prop.id == bid
            radius = body.prop.loci[1].position |> norm
            r = body.state.origin_frame.position
            gap = norm(x-r) - radius
            normal = (x-r)/norm(x-r)
            return gap, normal
        end
    end
end

function contact_gap_and_normal!(gaps, normals, local_positions, contact_ids, x::StaticArray{Tuple{N}, T, 1}, i, contact_body::RigidSphere, bodies, 
    ) where {T,N}
    (; bid, radius) = contact_body
    foreach(bodies) do body
        if (body.prop.id == bid)
            (;X) = body.state.origin_frame.axes
            r = body.state.origin_frame.position
            gap = norm(x - r) - radius
            normal = (x - r) / norm(x - r)
            gaps[i] = gap
            normals[i] = normal
            contact_ids[i] = i
            local_positions[i] = X\(x - r)
        end
    end
end

function contact_gap_and_normal!(gaps, normals, local_positions, contact_ids, x::StaticArray{Tuple{N}, T, 1}, i, contact_body::RigidCylinder, bodies, 
    ) where {T,N}
    (; bid, radius, halflen) = contact_body
    foreach(bodies) do body
        if (body.prop.id == bid)
            (;X) = body.state.origin_frame.axes
            rc = body.state.origin_frame.position
            p = x - rc
            p̄ = X\p
            px, py, pz = p̄
            if abs(px) > halflen
                gap = 1.0
                normal = SVector(1.0, 0.0, 0.0)
            elseif norm(SVector(py, pz)) > radius
                gap = 1.0
                normal = SVector(1.0, 0.0, 0.0)
            else
                gap = norm(SVector(py, pz)) - radius
                normal = X*SVector(0.0, py, pz) 
            end
            if gap < 0 
                ## @show gap
            end
            gaps[i] = gap
            normals[i] = normal / norm(normal)
            contact_ids[i] = i
            local_positions[i] = p̄
        end
    end
end

function contact_gap_and_normal!(gaps, normals, local_positions, contact_ids, x::StaticArray{Tuple{N}, T, 1}, i, contact_body::RigidCube, bodies, 
    ) where {T,N}
    (; bid, length, width, height) = contact_body
    hx, hy, hz = length, width, height
    foreach(bodies) do body
        if (body.prop.id == bid)
            (;X) = body.state.origin_frame.axes
            rc = body.state.origin_frame.position
            p = x - rc
            p̄ = X\p
            px, py, pz = p̄
            if (abs(px) > hx) || (abs(py) > hy) || (abs(pz) > hz)
                gap = 1.0
                normal = SVector(1.0, 0.0, 0.0)
            else
                normal_vector = [SVector(1.0, 0.0, 0.0),
                                    SVector(-1.0, 0.0, 0.0),
                                    SVector(0.0, 1.0, 0.0),
                                    SVector(0.0, -1.0, 0.0),
                                    SVector(0.0, 0.0, 1.0),
                                    SVector(0.0, 0.0, -1.0)]
                gap_vector = [hx-px, px+hx, hy-py, py+hy, hz-pz, pz+hz]
                min_idx = argmin(gap_vector)
                gap = - gap_vector[min_idx]
                normal = X*normal_vector[min_idx]
                # @show px, py, pz
                # @show gap
            end
            gaps[i] = gap
            normals[i] = normal
            contact_ids[i] = i
            local_positions[i] = p̄
        end
    end
end

function contact_gap_and_normal(x::StaticArray{Tuple{N}, T, 1}, contact_bodies::Vector{<:ContactRigid}, bodies) where {T,N}
    ncb = length(contact_bodies)
    gaps = zeros(T, ncb)
    normals = Vector{SVector{3,T}}(undef, ncb)
    local_positions = Vector{SVector{3,T}}(undef, ncb)
    contact_ids = Vector{Int}(undef, ncb)
    for (i, contact_body) in enumerate(contact_bodies)
        contact_gap_and_normal!(
            gaps, normals, local_positions, contact_ids,
            x, i,contact_body, bodies, 
        )
    end
    if all(gaps .<= 0)
        icb = argmax(gaps)
        contact_bodies[icb].pid += 1
        return gaps[icb], normals[icb], contact_ids[icb], local_positions[icb]
    else
        idx = findall((x) -> x > 0, gaps)
        positive_gaps = @view gaps[idx]
        positive_normals = @view normals[idx]
        positive_contact_body_ids = @view contact_ids[idx]
        icb = argmin(positive_gaps)
        return positive_gaps[icb], positive_normals[icb], positive_contact_body_ids[icb], local_positions[icb]
    end
end