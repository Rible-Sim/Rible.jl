#todo parameterization of joints
abstract type AbstractJoint end
struct NoJoint <: AbstractJoint end

"""
Linear constraint joint.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct LinearJoint{T,bodyType} <: AbstractJoint
    body::bodyType
    num_of_cstr::Int
    A::Matrix{T}
    violations::Vector{T}
end

"""
$(TYPEDSIGNATURES)
"""
function LinearJoint(body,A,violations)
    num_of_cstr = size(A,1)
    LinearJoint(body,num_of_cstr,A,violations)
end

function get_joint_idx(joint::LinearJoint)
    (;body) = joint
    coords = body.coords
    ncoords = get_num_of_coords(coords)
    full_idx = collect(1:ncoords)
    full_idx
end

"""
Joint that represents the body's own degrees of freedom (unconstrained).
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct BodyJoint{T,bodyType} <: AbstractJoint
    num_of_cstr::Int
    body::bodyType
end

"""
$(TYPEDSIGNATURES)
"""
BodyJoint(body) = BodyJoint{get_numbertype(body),typeof(body)}(0,body)

function get_joint_idx(joint::BodyJoint)
    (;body) = joint
    coords = body.coords
    ncoords = get_num_of_coords(coords)
    full_idx = collect(1:ncoords)
    full_idx
end

get_numbertype(::BodyJoint{T}) where{T} = T

get_num_of_dims(joint::BodyJoint) = get_num_of_dims(joint.body)

get_num_of_cstr(joint::AbstractJoint) = joint.num_of_cstr

get_num_of_cstr(joint::BodyJoint) = 0

"""
Joint that fixes a point on a body in space.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct FixedPointJoint{T,N,bodyType} <: AbstractJoint
    body::bodyType
    num_of_cstr::Int
    A::Matrix{T}
    r̄o::SVector{N,T} 
    ro::SVector{N,T}
end

function FixedPointJoint(body::AbstractBody{N,T}, r̄o) where {N,T}
    (;prop,coords) = body
    q,_ = body_state2coords_state(body)
    c = to_local_coords(coords,r̄o)
    A = to_position_jacobian(coords,q,c) |> Matrix
    ro = to_position(coords,q,c)
    FixedPointJoint(
        body,
        N,
        A,
        SVector{N}(r̄o),
        SVector{N}(ro),
    )
end

function get_joint_info(joint_type::Symbol)
    (joint_type == :FloatingSpherical)  && (return (ntrl = 3, nrot = 3, num_of_dof = 6, num_of_cstr = 0, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[]  ))
    (joint_type == :OrbitalSpherical)   && (return (ntrl = 2, nrot = 3, num_of_dof = 5, num_of_cstr = 1, mask_1st = Int[]  , mask_2nd = [1]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[]  ))
    (joint_type == :PlanarSpherical)    && (return (ntrl = 2, nrot = 3, num_of_dof = 5, num_of_cstr = 1, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [3]  , mask_3rd_egg = Int[], mask_4th = Int[]  ))
    (joint_type == :PrismaticSpherical) && (return (ntrl = 1, nrot = 3, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [1,2], mask_3rd_egg = Int[], mask_4th = Int[]  ))
    (joint_type == :Spherical)          && (return (ntrl = 0, nrot = 3, num_of_dof = 3, num_of_cstr = 3, mask_1st = [1,2,3], mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[]  ))
    (joint_type == :FloatingUniversal)  && (return (ntrl = 3, nrot = 2, num_of_dof = 5, num_of_cstr = 1, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [3]    )) #t'*b
    (joint_type == :OrbitalUniversal)   && (return (ntrl = 2, nrot = 2, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[]  , mask_2nd = [1]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [3]    )) #t'*b
    (joint_type == :PlanarUniversal)    && (return (ntrl = 2, nrot = 2, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [3]  , mask_3rd_egg = Int[], mask_4th = [3]    )) #t'*b
    (joint_type == :PrismaticUniversal) && (return (ntrl = 1, nrot = 2, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [1,2], mask_3rd_egg = Int[], mask_4th = [3]    )) #t'*b
    (joint_type == :UniversalPrismatic) && (return (ntrl = 1, nrot = 2, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = [1,2], mask_4th = [3]    )) #t'*b
    (joint_type == :Universal)          && (return (ntrl = 0, nrot = 2, num_of_dof = 2, num_of_cstr = 4, mask_1st = [1,2,3], mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [3]    )) #t'*b
    (joint_type == :FloatingRevolute)   && (return (ntrl = 3, nrot = 1, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2]  )) #t'*n, b'*n
    (joint_type == :OrbitalRevolute)    && (return (ntrl = 2, nrot = 1, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = [1]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2]  )) #t'*n, b'*n
    (joint_type == :PlanarRevolute)     && (return (ntrl = 2, nrot = 1, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [3]  , mask_3rd_egg = Int[], mask_4th = [1,2]  )) #t'*n, b'*n
    (joint_type == :Cylindrical)        && (return (ntrl = 1, nrot = 1, num_of_dof = 2, num_of_cstr = 4, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [1,2], mask_3rd_egg = Int[], mask_4th = [1,2]  )) #t'*n, b'*n
    (joint_type == :Revolute)           && (return (ntrl = 0, nrot = 1, num_of_dof = 1, num_of_cstr = 5, mask_1st = [1,2,3], mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2]  )) #t'*n, b'*n
    (joint_type == :Floating)           && (return (ntrl = 3, nrot = 0, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t'*n, b'*n, t'*b
    (joint_type == :Orbital)            && (return (ntrl = 2, nrot = 0, num_of_dof = 2, num_of_cstr = 4, mask_1st = Int[]  , mask_2nd = [1]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t'*n, b'*n, t'*b
    (joint_type == :Planar)             && (return (ntrl = 2, nrot = 0, num_of_dof = 2, num_of_cstr = 4, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [3]  , mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t'*n, b'*n, t'*b
    (joint_type == :Prismatic)          && (return (ntrl = 1, nrot = 0, num_of_dof = 1, num_of_cstr = 5, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [1,2], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t'*n, b'*n, t'*b
    (joint_type == :Fixed)              && (return (ntrl = 0, nrot = 0, num_of_dof = 0, num_of_cstr = 6, mask_1st = [1,2,3], mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t'*n, b'*n, t'*b
    error("Unknown joint type: $joint_type")
end

# (joint_type == :Revolute)           && (return (ntrl = 0, nrot = 1, num_of_dof = 1, num_of_cstr = 5, mask_1st = [1,2,3], mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2]  )) #t'*n, b'*n
# (joint_type == :FloatingUniversal)  && (return (ntrl = 3, nrot = 2, num_of_dof = 5, num_of_cstr = 1, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [3]    )) #t'*b

# (joint_type == :Prismatic)          && (return (ntrl = 1, nrot = 0, num_of_dof = 1, num_of_cstr = 5, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [1,2], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t'*n, b'*n, t'*b
# (joint_type == :PlanarSpherical)    && (return (ntrl = 2, nrot = 3, num_of_dof = 5, num_of_cstr = 1, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [3]  , mask_3rd_egg = Int[], mask_4th = Int[]  ))

"""
Generic joint prototype used for standard joint types.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct PrototypeJoint{hen2eggType,maskType,valueType,cacheType,posType} <: AbstractJoint
    hen2egg::hen2eggType
    num_of_cstr::Int
    num_of_dof::Int
    mask_1st::maskType
    mask_2nd::maskType
    mask_3rd::maskType
    mask_4th::maskType
    violations::valueType
    cache::cacheType
    loci_position_hen::posType
    loci_position_egg::posType
end


function ProtoJoint(hen2egg,joint_type::Symbol)
    (;hen,egg) = hen2egg
    joint_info = get_joint_info(joint_type)
    (;  ntrl, nrot, 
        num_of_dof, num_of_cstr, 
        mask_1st, 
        mask_2nd,
        mask_3rd_hen, 
        mask_3rd_egg, 
        mask_4th, 
    ) = joint_info
    mask_3rd = vcat(mask_3rd_hen,mask_3rd_egg.+3)
    coords_hen = hen.body.coords
    coords_egg = egg.body.coords

    state_hen = hen.body.state
    state_egg = egg.body.state
    q_hen,_ = cartesian_frame2coords(coords_hen,state_hen.origin_frame)
    q_egg,_ = cartesian_frame2coords(coords_egg,state_egg.origin_frame)
    
    cache, violations = build_joint_cache(
        coords_hen,coords_egg,
        hen.position,
        egg.position,
        hen.trl_axes,
        egg.trl_axes,
        hen.rot_axes,
        egg.rot_axes,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg
    )

    # @show mask_1st,mask_2nd,mask_3rd,mask_4th
    # @show values
    PrototypeJoint(
        hen2egg,
        num_of_cstr,
        num_of_dof,
        mask_1st,
        mask_2nd,
        mask_3rd,
        mask_4th,
        violations,
        cache,
        hen.position,
        egg.position
    )
end

get_numbertype(::FixedPointJoint{T}) where{T} = T


get_numbertype(::LinearJoint{T}) where{T} = T


function get_numbertype(joint::PrototypeJoint)
    get_numbertype(joint.hen2egg.hen.body)
end

function get_joint_idx(joint::PrototypeJoint)
    (;hen,egg) = joint.hen2egg
    coords_hen = hen.body.coords
    coords_egg = egg.body.coords
    ncoords_hen = get_num_of_coords(coords_hen)
    ncoords_egg = get_num_of_coords(coords_egg)
    full_idx = vcat(
        collect(1:ncoords_hen),
        collect(1:ncoords_egg) .+ ncoords_hen
    )
    full_idx
end

get_num_of_dims(joint::PrototypeJoint) = get_num_of_dims(joint.hen2egg.hen.body)
