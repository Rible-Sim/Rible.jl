#done use the basic types of Constraints
#todo parameterization of joints
#note can full rotation cstr be linear?
abstract type ExtrinsicConstraints{T} end


get_numbertype(cst::ExtrinsicConstraints{<:AbstractArray{T}}) where T = T


"""
刚体通用线性约束类。
$(TYPEDEF)
"""
struct LinearJoint{T} <: ExtrinsicConstraints{T}
    id::Int
    num_of_cstr::Int
    free_idx::Vector{Int64}
    sys_free_idx::Vector{Int64}
    A::Matrix{T}
    violations::Vector{T}
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function LinearJoint(id,indexed,A,violations)
    num_of_cstr = size(A,1)
    free_idx = collect(1:size(cst.A,2))
    sys_free_idx = collect(1:indexed.num_of_free_coords)
    LinearJoint(id,num_of_cstr,free_idx,sys_free_idx,A,violations)
end

function cstr_function(cst::LinearJoint,structure,q,c)
    (;A,sys_free_idx,violations) = cst
    A*q[sys_free_idx] .- violations
end

function cstr_jacobian(cst::LinearJoint,structure,q)
    cst.A
end

function cstr_forces_jacobian(cst::LinearJoint,st,q,λ)
    (;sys_free_idx) = cst
    n = length(sys_free_idx)
    zeros(eltype(λ),n,n)
end

function cstr_velocity_jacobian(cst::LinearJoint,st,q,q̇)
    (;num_of_cstr,sys_free_idx) = cst
    n = length(sys_free_idx)
    zeros(eltype(q̇),num_of_cstr,n)
end

"""
固定约束构造子。
$(TYPEDSIGNATURES)
"""
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

FloatingSphericalJoint(id,indexed,hen2egg)  = PrototypeJoint(id,indexed,hen2egg,:FloatingSpherical)
OrbitalSphericalJoint(id,indexed,hen2egg)   = PrototypeJoint(id,indexed,hen2egg,:OrbitalSpherical)
PlanarSphericalJoint(id,indexed,hen2egg)    = PrototypeJoint(id,indexed,hen2egg,:PlanarSpherical)
PrismaticSphericalJoint(id,indexed,hen2egg) = PrototypeJoint(id,indexed,hen2egg,:PrismaticSpherical)
SphericalJoint(id,indexed,hen2egg)          = PrototypeJoint(id,indexed,hen2egg,:Spherical)
FloatingUniversalJoint(id,indexed,hen2egg)  = PrototypeJoint(id,indexed,hen2egg,:FloatingUniversal)
OrbitalUniversalJoint(id,indexed,hen2egg)   = PrototypeJoint(id,indexed,hen2egg,:OrbitalUniversal)
PlanarUniversalJoint(id,indexed,hen2egg)    = PrototypeJoint(id,indexed,hen2egg,:PlanarUniversal)
PrismaticUniversalJoint(id,indexed,hen2egg) = PrototypeJoint(id,indexed,hen2egg,:PrismaticUniversal)
UniversalPrismaticJoint(id,indexed,hen2egg) = PrototypeJoint(id,indexed,hen2egg,:UniversalPrismatic)
UniversalJoint(id,indexed,hen2egg)          = PrototypeJoint(id,indexed,hen2egg,:Universal)
FloatingRevoluteJoint(id,indexed,hen2egg)   = PrototypeJoint(id,indexed,hen2egg,:FloatingRevolute)
OrbitalRevoluteJoint(id,indexed,hen2egg)    = PrototypeJoint(id,indexed,hen2egg,:OrbitalRevolute)
PlanarRevoluteJoint(id,indexed,hen2egg)     = PrototypeJoint(id,indexed,hen2egg,:PlanarRevolute)
CylindricalJoint(id,indexed,hen2egg)        = PrototypeJoint(id,indexed,hen2egg,:Cylindrical)
RevoluteJoint(id,indexed,hen2egg)           = PrototypeJoint(id,indexed,hen2egg,:Revolute)
FloatingJoint(id,indexed,hen2egg)           = PrototypeJoint(id,indexed,hen2egg,:Floating)
OrbitalJoint(id,indexed,hen2egg)            = PrototypeJoint(id,indexed,hen2egg,:Orbital)
PlanarJoint(id,indexed,hen2egg)             = PrototypeJoint(id,indexed,hen2egg,:Planar)
PrismaticJoint(id,indexed,hen2egg)          = PrototypeJoint(id,indexed,hen2egg,:Prismatic)
FixedJoint(id,indexed,hen2egg)              = PrototypeJoint(id,indexed,hen2egg,:Fixed)

const PinJoint = SphericalJoint


"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct PrototypeJoint{hen2eggType,maskType,valueType,cacheType} <: ExtrinsicConstraints{valueType}
    id::Int
    hen2egg::hen2eggType
    num_of_cstr::Int
    num_of_dof::Int
    free_idx::Vector{Int}
    sys_free_idx::Vector{Int}
    mask_1st::maskType
    mask_2nd::maskType
    mask_3rd::maskType
    mask_4th::maskType
    violations::valueType
    cache::cacheType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function PrototypeJoint(id,indexed,hen2egg,joint_type::Symbol) 
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
    (;bodyid2sys_free_coords,bodyid2sys_full_coords,num_of_free_coords) = indexed
    mask_3rd = vcat(mask_3rd_hen,mask_3rd_egg.+3)
    nmcs_hen = hen.rbsig.coords.nmcs
    nmcs_egg = egg.rbsig.coords.nmcs
    free_idx_hen = hen.rbsig.coords.free_idx
    free_idx_egg = egg.rbsig.coords.free_idx
    ncoords_hen = NCF.get_num_of_coords(nmcs_hen)
    ncoords_egg = NCF.get_num_of_coords(nmcs_egg)
    free_idx = vcat(
        free_idx_hen,
        free_idx_egg .+ ncoords_hen
    )
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    free_hen = bodyid2sys_free_coords[id_hen]
    free_egg = bodyid2sys_free_coords[id_egg]
    sys_free_idx = vcat(
        free_hen,
        free_egg
    )

    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    q_hen,_ = cartesian_frame2coords(nmcs_hen,state_hen.origin_frame)
    q_egg,_ = cartesian_frame2coords(nmcs_egg,state_egg.origin_frame)
    
    cache, values = build_joint_cache(
        nmcs_hen,nmcs_egg,
        hen.rbsig.prop.loci[hen.pid].position,
        egg.rbsig.prop.loci[egg.pid].position,
        hen.rbsig.prop.loci[hen.trlid].axes,
        egg.rbsig.prop.loci[egg.trlid].axes,
        hen.rbsig.prop.loci[hen.rotid].axes,
        egg.rbsig.prop.loci[egg.rotid].axes,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg
    )
    # @show mask_1st,mask_2nd,mask_3rd,mask_4th
    # @show values
    PrototypeJoint(
        id,hen2egg,
        num_of_cstr,
        num_of_dof,
        free_idx,
        sys_free_idx,
        mask_1st,
        mask_2nd,
        mask_3rd,
        mask_4th,
        values,
        cache
    )
end

function cstr_function(cst::PrototypeJoint,structure::Structure,q, c = get_local_coords(structure))
    (;indexed,numbered) = structure.connectivity
    (;bodyid2sys_full_coords) = indexed
    (;sys_loci2coords_idx,bodyid2sys_loci_idx) = numbered
    (;
        num_of_cstr,hen2egg,
        mask_1st,
        violations,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = cst
    (;hen,egg) = hen2egg
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    nmcs_hen = hen.rbsig.coords.nmcs
    nmcs_egg = egg.rbsig.coords.nmcs
    T = eltype(q)
    ret = zeros(T,num_of_cstr)
    q_hen = @view q[bodyid2sys_full_coords[id_hen]]
    q_egg = @view q[bodyid2sys_full_coords[id_egg]]
    q = vcat(q_hen,q_egg)
    # cstr violations
    # @show cst.id
    get_joint_violations!(
        ret,
        nmcs_hen, nmcs_egg,
        hen.rbsig.prop.loci[hen.pid].position,
        egg.rbsig.prop.loci[egg.pid].position,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg,
        violations
    )
    ret
end

function cstr_jacobian(cst::PrototypeJoint,structure::Structure,q,c = get_local_coords(structure))
    (;indexed,numbered) = structure.connectivity
    (;bodyid2sys_free_coords,
      bodyid2sys_full_coords,
      num_of_free_coords,
      num_of_full_coords,
    ) = indexed
    (;sys_loci2coords_idx,bodyid2sys_loci_idx) = numbered
    (;
        num_of_cstr,
        free_idx,
        sys_free_idx,
        hen2egg,
        mask_1st,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = cst
    (;hen,egg) = hen2egg
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    nmcs_hen = hen.rbsig.coords.nmcs
    nmcs_egg = egg.rbsig.coords.nmcs
    T = eltype(q)
    ret = zeros(T,num_of_cstr,length(sys_free_idx))
    q_hen = @view q[bodyid2sys_full_coords[id_hen]]
    q_egg = @view q[bodyid2sys_full_coords[id_egg]]
    # translate
    get_joint_jacobian!(
        ret,
        nmcs_hen, nmcs_egg,
        hen.rbsig.prop.loci[hen.pid].position,
        egg.rbsig.prop.loci[egg.pid].position,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        free_idx,
        q_hen,q_egg
    )
    ret
end

function cstr_forces_jacobian(cst::PrototypeJoint,structure,q,λ)
    (;indexed,numbered) = structure.connectivity
    (;bodyid2sys_free_coords,
      bodyid2sys_full_coords,
      num_of_free_coords,
      num_of_full_coords,
    ) = indexed
    (;
        num_of_cstr,
        hen2egg,
        free_idx,
        mask_1st,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = cst
    (;hen,egg) = hen2egg
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    q_hen = @view q[bodyid2sys_full_coords[id_hen]]
    q_egg = @view q[bodyid2sys_full_coords[id_egg]]
    nmcs_hen = hen.rbsig.coords.nmcs
    nmcs_egg = egg.rbsig.coords.nmcs
    T = eltype(λ)
    num_of_free_idx = length(free_idx)
    ret = zeros(T,num_of_free_idx,num_of_free_idx)
    get_joint_forces_jacobian!(
        ret,
        num_of_cstr,
        nmcs_hen, nmcs_egg,
        hen.rbsig.prop.loci[hen.pid].position,
        egg.rbsig.prop.loci[egg.pid].position,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        free_idx,
        q_hen,q_egg,
        λ
    )
    ret
end

function cstr_velocity_jacobian(cst::PrototypeJoint,structure,q,q̇)
    (;indexed,numbered) = structure.connectivity
    (;bodyid2sys_free_coords,
      bodyid2sys_full_coords,
      num_of_free_coords,
      num_of_full_coords,
    ) = indexed
    (;
        num_of_cstr,
        hen2egg,
        free_idx,
        mask_1st,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = cst
    (;hen,egg) = hen2egg
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    q_hen = @view q[bodyid2sys_full_coords[id_hen]]
    q_egg = @view q[bodyid2sys_full_coords[id_egg]]
    q̇_hen = @view q̇[bodyid2sys_full_coords[id_hen]]
    q̇_egg = @view q̇[bodyid2sys_full_coords[id_egg]]
    nmcs_hen = hen.rbsig.coords.nmcs
    nmcs_egg = egg.rbsig.coords.nmcs
    T = eltype(q̇)
    num_of_free_idx = length(free_idx)
    ret = zeros(T,num_of_cstr,num_of_free_idx)
    get_joint_velocity_jacobian!(
        ret,
        num_of_cstr,
        nmcs_hen, nmcs_egg,
        hen.rbsig.prop.loci[hen.pid].position,
        egg.rbsig.prop.loci[egg.pid].position,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        free_idx,
        q_hen,q_egg,
        q̇_hen,q̇_egg,
    )
    ret
end
