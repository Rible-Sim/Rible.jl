#todo parameterization of joints
abstract type AbstractJoint end

"""
$(TYPEDEF)
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
function LinearJoint(id,body,A,violations)
    num_of_cstr = size(A,1)
    LinearJoint(id,body,num_of_cstr,A,violations)
end

function get_joint_idx(joint::LinearJoint)
    (;body) = joint
    nmcs = body.coords.nmcs
    free_idx = body.coords.free_idx
    ncoords = get_num_of_coords(nmcs)
    full_idx = collect(1:ncoords)
    free_idx = free_idx
    full_idx, free_idx
end

function get_appar_idx(appar::Apparatus{<:LinearJoint},bodyid2sys_free_coords)
    (;body) = appar.joint
    bid = body.prop.id
    bodyid2sys_free_coords[bid]
end

function FixedBodyConstraint(id::Int,body::AbstractBody)
    nmcs = body.coords.nmcs
    state = body.state
    q,_ = cartesian_frame2coords(
        nmcs,
        state.origin_frame
    )
    independent_free_idx = find_independent_free_idx(nmcs,q)
    violations = q[independent_free_idx]
    FixedIndicesConstraint(id,body,independent_free_idx,violations)
end

function FixedIndicesConstraint(id::Int,body,idx,violations)
    num_of_cstr = length(idx)
    num_of_coords = get_num_of_coords(body)
    T = get_numbertype(body)
    A = zeros(T,num_of_cstr,num_of_coords)
    for (i,j) in enumerate(idx)
        A[i,j] = 1
    end
    joint = LinearJoint(body,num_of_cstr,A,violations)
    full_idx, free_idx = get_joint_idx(joint)
    force = nothing
    Apparatus(
        id,
        joint,
        force,
        0,
        full_idx, free_idx
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
    (joint_type == :Floating)           && (return (ntrl = 3, nrot = 0, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t'*b, t'*n, b'*n
    (joint_type == :Orbital)            && (return (ntrl = 2, nrot = 0, num_of_dof = 2, num_of_cstr = 4, mask_1st = Int[]  , mask_2nd = [1]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t'*b, t'*n, b'*n
    (joint_type == :Planar)             && (return (ntrl = 2, nrot = 0, num_of_dof = 2, num_of_cstr = 4, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [3]  , mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t'*b, t'*n, b'*n
    (joint_type == :Prismatic)          && (return (ntrl = 1, nrot = 0, num_of_dof = 1, num_of_cstr = 5, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [1,2], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t'*b, t'*n, b'*n
    (joint_type == :Fixed)              && (return (ntrl = 0, nrot = 0, num_of_dof = 0, num_of_cstr = 6, mask_1st = [1,2,3], mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t'*b, t'*n, b'*n
end

FloatingSphericalJoint(id,hen2egg,force=nothing)  = PrototypeJoint(id,hen2egg,force,:FloatingSpherical)
OrbitalSphericalJoint(id,hen2egg,force=nothing)   = PrototypeJoint(id,hen2egg,force,:OrbitalSpherical)
PlanarSphericalJoint(id,hen2egg,force=nothing)    = PrototypeJoint(id,hen2egg,force,:PlanarSpherical)
PrismaticSphericalJoint(id,hen2egg,force=nothing) = PrototypeJoint(id,hen2egg,force,:PrismaticSpherical)
SphericalJoint(id,hen2egg,force=nothing)          = PrototypeJoint(id,hen2egg,force,:Spherical)
FloatingUniversalJoint(id,hen2egg,force=nothing)  = PrototypeJoint(id,hen2egg,force,:FloatingUniversal)
OrbitalUniversalJoint(id,hen2egg,force=nothing)   = PrototypeJoint(id,hen2egg,force,:OrbitalUniversal)
PlanarUniversalJoint(id,hen2egg,force=nothing)    = PrototypeJoint(id,hen2egg,force,:PlanarUniversal)
PrismaticUniversalJoint(id,hen2egg,force=nothing) = PrototypeJoint(id,hen2egg,force,:PrismaticUniversal)
UniversalPrismaticJoint(id,hen2egg,force=nothing) = PrototypeJoint(id,hen2egg,force,:UniversalPrismatic)
UniversalJoint(id,hen2egg,force=nothing)          = PrototypeJoint(id,hen2egg,force,:Universal)
FloatingRevoluteJoint(id,hen2egg,force=nothing)   = PrototypeJoint(id,hen2egg,force,:FloatingRevolute)
OrbitalRevoluteJoint(id,hen2egg,force=nothing)    = PrototypeJoint(id,hen2egg,force,:OrbitalRevolute)
PlanarRevoluteJoint(id,hen2egg,force=nothing)     = PrototypeJoint(id,hen2egg,force,:PlanarRevolute)
CylindricalJoint(id,hen2egg,force=nothing)        = PrototypeJoint(id,hen2egg,force,:Cylindrical)
RevoluteJoint(id,hen2egg,force=nothing)           = PrototypeJoint(id,hen2egg,force,:Revolute)
FloatingJoint(id,hen2egg,force=nothing)           = PrototypeJoint(id,hen2egg,force,:Floating)
OrbitalJoint(id,hen2egg,force=nothing)            = PrototypeJoint(id,hen2egg,force,:Orbital)
PlanarJoint(id,hen2egg,force=nothing)             = PrototypeJoint(id,hen2egg,force,:Planar)
PrismaticJoint(id,hen2egg,force=nothing)          = PrototypeJoint(id,hen2egg,force,:Prismatic)
FixedJoint(id,hen2egg,force=nothing)              = PrototypeJoint(id,hen2egg,force,:Fixed)

const PinJoint = SphericalJoint

"""
$(TYPEDEF)
"""
struct PrototypeJoint{hen2eggType,maskType,valueType,cacheType} <: AbstractJoint
    hen2egg::hen2eggType
    num_of_cstr::Int
    num_of_dof::Int
    mask_1st::maskType
    mask_2nd::maskType
    mask_3rd::maskType
    mask_4th::maskType
    violations::valueType
    cache::cacheType
end

"""
A prototype joint is defined for a `Hen2Egg`.

Translation uses either the axis on Hen or on Egg.

Rotation uses the axis on Egg, because the local angular velocity is defined relative to Egg.

The order of translation and rotation matters!

$(TYPEDSIGNATURES)
"""
function PrototypeJoint(id,hen2egg,force,joint_type::Symbol) 
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
    nmcs_hen = hen.body.coords.nmcs
    nmcs_egg = egg.body.coords.nmcs

    state_hen = hen.body.state
    state_egg = egg.body.state
    q_hen,_ = cartesian_frame2coords(nmcs_hen,state_hen.origin_frame)
    q_egg,_ = cartesian_frame2coords(nmcs_egg,state_egg.origin_frame)
    
    cache, values = build_joint_cache(
        nmcs_hen,nmcs_egg,
        hen.body.prop.loci[hen.pid].position,
        egg.body.prop.loci[egg.pid].position,
        hen.body.prop.loci[hen.trlid].axes,
        egg.body.prop.loci[egg.trlid].axes,
        hen.body.prop.loci[hen.rotid].axes,
        egg.body.prop.loci[egg.rotid].axes,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg
    )
    # @show mask_1st,mask_2nd,mask_3rd,mask_4th
    # @show values
    joint = PrototypeJoint(
        hen2egg,
        num_of_cstr,
        num_of_dof,
        mask_1st,
        mask_2nd,
        mask_3rd,
        mask_4th,
        values,
        cache
    )
    full_idx, free_idx = get_joint_idx(joint)
    Apparatus(
        id,
        joint,
        force,
        0,
        full_idx, free_idx
    )
end

"""
$(TYPEDEF)
"""
struct CableJoint{hen2eggType} <: AbstractJoint
    hen2egg::hen2eggType
    num_of_cstr::Int
end


function get_joint_idx(joint::Union{PrototypeJoint,CableJoint})
    (;hen,egg) = joint.hen2egg
    nmcs_hen = hen.body.coords.nmcs
    nmcs_egg = egg.body.coords.nmcs
    free_idx_hen = hen.body.coords.free_idx
    free_idx_egg = egg.body.coords.free_idx
    ncoords_hen = get_num_of_coords(nmcs_hen)
    ncoords_egg = get_num_of_coords(nmcs_egg)
    full_idx = vcat(
        collect(1:ncoords_hen),
        collect(1:ncoords_egg) .+ ncoords_hen
    )
    free_idx = vcat(
        free_idx_hen,
        free_idx_egg .+ ncoords_hen
    )
    full_idx, free_idx
end

function get_appar_idx(appar::Apparatus{<:jointType},bodyid2sys_free_coords) where {jointType<:Union{PrototypeJoint,CableJoint}}
    (;joint) = appar
    (;hen,egg) = joint.hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    free_hen = bodyid2sys_free_coords[id_hen]
    free_egg = bodyid2sys_free_coords[id_egg]
    sys_free_coords_idx = vcat(
        free_hen,
        free_egg
    )
    sys_free_coords_idx
end