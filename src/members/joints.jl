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

function get_joint_idx(cst::LinearJoint,bodyid2sys_free_coords)
    (;body) = cst
    bid = body.prop.id
    nmcs = body.coords.nmcs
    free_idx = body.coords.free_idx
    ncoords = get_num_of_coords(nmcs)
    full_idx = collect(1:ncoords)
    free_idx = free_idx
    sys_free_coords_idx = bodyid2sys_free_coords[bid]
    # full_idx, sys_full_idx
    full_idx, free_idx, sys_free_coords_idx
end

function FixedBodyConstraint(id::Int,indexed,body::AbstractBody)
    (;bodyid2sys_free_coords,num_of_free_coords) = indexed
    bodyid = body.prop.id
    nmcs = body.coords.nmcs
    state = body.state
    q,_ = cartesian_frame2coords(
        nmcs,
        state.origin_frame
    )
    independent_free_idx = find_independent_free_idx(nmcs,q)
    jointed_sys_free_idx = bodyid2sys_free_coords[bodyid][independent_free_idx]
    violations = q[independent_free_idx]
    num_of_cstr = num_of_free_idx = length(jointed_sys_free_idx)
    A = zeros(eltype(q),num_of_cstr,num_of_free_idx)
    for i = 1:num_of_cstr
        A[i,i] = 1
    end
    LinearJoint(
        id,
        num_of_cstr,
        independent_free_idx,
        jointed_sys_free_idx,
        A,
        violations
    )
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
    force = nothing
    Apparatus(
        id,
        Val(true),
        Val(false),
        joint,
        force,
    )
end

function get_joint_info(joint_type::Symbol)
    (joint_type == :FloatingSpherical)  && (return (ntrl = 3, nrot = 3, num_of_dof = 6, num_of_cstr = 0, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[]  ))
    (joint_type == :OrbitalSpherical)   && (return (ntrl = 2, nrot = 3, num_of_dof = 5, num_of_cstr = 1, mask_1st = Int[]  , mask_2nd = [1]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[]  ))
    (joint_type == :PlanarSpherical)    && (return (ntrl = 2, nrot = 3, num_of_dof = 5, num_of_cstr = 1, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [3]  , mask_3rd_egg = Int[], mask_4th = Int[]  ))
    (joint_type == :PrismaticSpherical) && (return (ntrl = 1, nrot = 3, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [1,2], mask_3rd_egg = Int[], mask_4th = Int[]  ))
    (joint_type == :Spherical)          && (return (ntrl = 0, nrot = 3, num_of_dof = 3, num_of_cstr = 3, mask_1st = [1,2,3], mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[]  ))
    (joint_type == :FloatingUniversal)  && (return (ntrl = 3, nrot = 2, num_of_dof = 5, num_of_cstr = 1, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [3]    )) #t1'*t2
    (joint_type == :OrbitalUniversal)   && (return (ntrl = 2, nrot = 2, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[]  , mask_2nd = [1]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [3]    )) #t1'*t2
    (joint_type == :PlanarUniversal)    && (return (ntrl = 2, nrot = 2, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [3]  , mask_3rd_egg = Int[], mask_4th = [3]    )) #t1'*t2
    (joint_type == :PrismaticUniversal) && (return (ntrl = 1, nrot = 2, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [1,2], mask_3rd_egg = Int[], mask_4th = [3]    )) #t1'*t2
    (joint_type == :UniversalPrismatic) && (return (ntrl = 1, nrot = 2, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = [1,2], mask_4th = [3]    )) #t1'*t2
    (joint_type == :Universal)          && (return (ntrl = 0, nrot = 2, num_of_dof = 2, num_of_cstr = 4, mask_1st = [1,2,3], mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [3]    )) #t1'*t2
    (joint_type == :FloatingRevolute)   && (return (ntrl = 3, nrot = 1, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2]  )) #t1'*n, t2'*n
    (joint_type == :OrbitalRevolute)    && (return (ntrl = 2, nrot = 1, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = [1]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2]  )) #t1'*n, t2'*n
    (joint_type == :PlanarRevolute)     && (return (ntrl = 2, nrot = 1, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [3]  , mask_3rd_egg = Int[], mask_4th = [1,2]  )) #t1'*n, t2'*n
    (joint_type == :Cylindrical)        && (return (ntrl = 1, nrot = 1, num_of_dof = 2, num_of_cstr = 4, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [1,2], mask_3rd_egg = Int[], mask_4th = [1,2]  )) #t1'*n, t2'*n
    (joint_type == :Revolute)           && (return (ntrl = 0, nrot = 1, num_of_dof = 1, num_of_cstr = 5, mask_1st = [1,2,3], mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2]  )) #t1'*n, t2'*n
    (joint_type == :Floating)           && (return (ntrl = 3, nrot = 0, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Orbital)            && (return (ntrl = 2, nrot = 0, num_of_dof = 2, num_of_cstr = 4, mask_1st = Int[]  , mask_2nd = [1]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Planar)             && (return (ntrl = 2, nrot = 0, num_of_dof = 2, num_of_cstr = 4, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [3]  , mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Prismatic)          && (return (ntrl = 1, nrot = 0, num_of_dof = 1, num_of_cstr = 5, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [1,2], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Fixed)              && (return (ntrl = 0, nrot = 0, num_of_dof = 0, num_of_cstr = 6, mask_1st = [1,2,3], mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t1'*t2, t1'*n, t2'*n
end

FloatingSphericalJoint(id,hen2egg)  = PrototypeJoint(id,hen2egg,:FloatingSpherical)
OrbitalSphericalJoint(id,hen2egg)   = PrototypeJoint(id,hen2egg,:OrbitalSpherical)
PlanarSphericalJoint(id,hen2egg)    = PrototypeJoint(id,hen2egg,:PlanarSpherical)
PrismaticSphericalJoint(id,hen2egg) = PrototypeJoint(id,hen2egg,:PrismaticSpherical)
SphericalJoint(id,hen2egg)          = PrototypeJoint(id,hen2egg,:Spherical)
FloatingUniversalJoint(id,hen2egg)  = PrototypeJoint(id,hen2egg,:FloatingUniversal)
OrbitalUniversalJoint(id,hen2egg)   = PrototypeJoint(id,hen2egg,:OrbitalUniversal)
PlanarUniversalJoint(id,hen2egg)    = PrototypeJoint(id,hen2egg,:PlanarUniversal)
PrismaticUniversalJoint(id,hen2egg) = PrototypeJoint(id,hen2egg,:PrismaticUniversal)
UniversalPrismaticJoint(id,hen2egg) = PrototypeJoint(id,hen2egg,:UniversalPrismatic)
UniversalJoint(id,hen2egg)          = PrototypeJoint(id,hen2egg,:Universal)
FloatingRevoluteJoint(id,hen2egg)   = PrototypeJoint(id,hen2egg,:FloatingRevolute)
OrbitalRevoluteJoint(id,hen2egg)    = PrototypeJoint(id,hen2egg,:OrbitalRevolute)
PlanarRevoluteJoint(id,hen2egg)     = PrototypeJoint(id,hen2egg,:PlanarRevolute)
CylindricalJoint(id,hen2egg)        = PrototypeJoint(id,hen2egg,:Cylindrical)
RevoluteJoint(id,hen2egg)           = PrototypeJoint(id,hen2egg,:Revolute)
FloatingJoint(id,hen2egg)           = PrototypeJoint(id,hen2egg,:Floating)
OrbitalJoint(id,hen2egg)            = PrototypeJoint(id,hen2egg,:Orbital)
PlanarJoint(id,hen2egg)             = PrototypeJoint(id,hen2egg,:Planar)
PrismaticJoint(id,hen2egg)          = PrototypeJoint(id,hen2egg,:Prismatic)
FixedJoint(id,hen2egg)              = PrototypeJoint(id,hen2egg,:Fixed)

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
$(TYPEDSIGNATURES)
"""
function PrototypeJoint(id,hen2egg,joint_type::Symbol) 
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
    nmcs_hen = hen.bodysig.coords.nmcs
    nmcs_egg = egg.bodysig.coords.nmcs

    state_hen = hen.bodysig.state
    state_egg = egg.bodysig.state
    q_hen,_ = cartesian_frame2coords(nmcs_hen,state_hen.origin_frame)
    q_egg,_ = cartesian_frame2coords(nmcs_egg,state_egg.origin_frame)
    
    cache, values = build_joint_cache(
        nmcs_hen,nmcs_egg,
        hen.bodysig.prop.loci[hen.pid].position,
        egg.bodysig.prop.loci[egg.pid].position,
        hen.bodysig.prop.loci[hen.trlid].axes,
        egg.bodysig.prop.loci[egg.trlid].axes,
        hen.bodysig.prop.loci[hen.rotid].axes,
        egg.bodysig.prop.loci[egg.rotid].axes,
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
    force = nothing
    Apparatus(
        id,
        Val(true),
        Val(false),
        joint,
        force,
    )
end

"""
$(TYPEDEF)
"""
struct CableJoint{hen2eggType} <: AbstractJoint
    hen2egg::hen2eggType
    num_of_cstr::Int
end


function get_joint_idx(cst::Union{PrototypeJoint,CableJoint},bodyid2sys_free_coords)
    (;hen,egg) = cst.hen2egg
    id_hen = hen.bodysig.prop.id
    id_egg = egg.bodysig.prop.id
    nmcs_hen = hen.bodysig.coords.nmcs
    nmcs_egg = egg.bodysig.coords.nmcs
    free_idx_hen = hen.bodysig.coords.free_idx
    free_idx_egg = egg.bodysig.coords.free_idx
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
    free_hen = bodyid2sys_free_coords[id_hen]
    free_egg = bodyid2sys_free_coords[id_egg]
    sys_free_coords_idx = vcat(
        free_hen,
        free_egg
    )
    # full_idx, sys_full_idx
    full_idx, free_idx, sys_free_coords_idx
end