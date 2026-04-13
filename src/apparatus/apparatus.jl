mutable struct ApparatusCache{φType,AType,FJType,FType,QType,KType,CType,cacheType}
    # cstr
    φ::φType
    # cstr Jacobian
    A::AType
    # cstr force Jacobian, symmetric
    ∂Aᵀλ∂q::FJType
    # cstr velocity Jacobian
    ∂Aq̇∂q::FType
    # generalized force
    Q::QType
    #∂Q∂q
    K::KType
    #∂Q∂q̇
    C::CType
    named::cacheType
end

@inline function zero!(A::SparseMatrixCSC)
    fill!(A.nzval, zero(eltype(A)))
end

@inline function zero!(A::AbstractArray)
    fill!(A, zero(eltype(A)))
end

function ApparatusCache(joint::AbstractJoint,force::AbstractForce,full_coords_idx, )
    T = get_numbertype(joint)
    num_of_dims = get_num_of_dims(joint)
    num_of_cstr = get_num_of_cstr(joint)
    num_appar_full_idx = length(full_coords_idx)
    
    φ = zeros(T,num_of_cstr)
    A = spzeros(T,num_of_cstr,num_appar_full_idx)
    ∂Aᵀλ∂q = spzeros(T,num_appar_full_idx,num_appar_full_idx)
    ∂Aq̇∂q = spzeros(T,num_of_cstr,num_appar_full_idx)
    Q = zeros(T,num_appar_full_idx)
    K = spzeros(T,num_appar_full_idx,num_appar_full_idx)
    C = spzeros(T,num_appar_full_idx,num_appar_full_idx)

    J = zeros(T,num_of_dims,num_appar_full_idx)
    DJ = zeros(T,num_of_dims,num_appar_full_idx)
    named = @eponymtuple(
        J,DJ
    )
    
    cache = ApparatusCache(
        φ,
        A,
        ∂Aᵀλ∂q,
        ∂Aq̇∂q,
        Q,
        K,
        C,
        named
    )
end


struct Apparatus{jointType,forceType,cacheType}
    id::Int
    joint::jointType
    force::forceType
    num_of_aux_var::Int
    full_coords_idx::Vector{Int}
    free_coords_idx::Vector{Int}
    cache::cacheType
end

function Base.isless(a::Apparatus,b::Apparatus)
    isless(a.id,b.id)
end

get_id(appar::Apparatus) = appar.id

function get_appar_idx(appar::Apparatus{<:LinearJoint},bodyid2sys_full_coords)
    (;body) = appar.joint
    bid = body.prop.id
    bodyid2sys_full_coords[bid]
end

function get_appar_idx(appar::Apparatus{<:FixedPointJoint},bodyid2sys_full_coords)
    (;body) = appar.joint
    bid = body.prop.id
    bodyid2sys_full_coords[bid]
end


function get_appar_idx(appar::Apparatus,bodyid2sys_full_coords)
    Int[]
end

function get_appar_idx(appar::Apparatus{<:jointType},bodyid2sys_full_coords) where {jointType<:PrototypeJoint}
    (;joint) = appar
    (;hen,egg) = joint.hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    full_hen = bodyid2sys_full_coords[id_hen]
    full_egg = bodyid2sys_full_coords[id_egg]
    appar_full_coords_idx = vcat(
        full_hen,
        full_egg
    )
end

function FixedIndicesApparatus(id::Int,body,idx,violations)
    num_of_cstr = length(idx)
    num_of_coords = get_num_of_coords(body)
    T = get_numbertype(body)
    A = zeros(T,num_of_cstr,num_of_coords)
    for (i,j) in enumerate(idx)
        A[i,j] = 1
    end
    joint = LinearJoint(body,num_of_cstr,A,violations)
    full_idx = get_joint_idx(joint)
    force = NoForce()
    Apparatus(
        id,
        joint,
        force,
        0,
        full_idx,full_idx,
        ApparatusCache(joint,force,full_idx,)
    )
end

function FixedBodyApparatus(id::Int,body::AbstractBody)
    state = body.state
    q,_ = cartesian_frame2coords(
        body.coords,
        state.origin_frame
    )
    independent_free_idx = find_independent_free_idx(body.coords,q)
    violations = q[independent_free_idx]
    FixedIndicesApparatus(id,body,independent_free_idx,violations)
end


"""
A prototype joint is defined for a `Hen2Egg`.

Translation uses either the axis on Hen or on Egg.

Rotation uses the axis on Egg, because the local angular velocity is defined relative to Egg.

The order of translation and rotation matters!

$(TYPEDSIGNATURES)
"""
function proto_joint_apparatus(id,hen2egg,force,joint_type::Symbol;
    ) 
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
    joint = PrototypeJoint(
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

    full_idx = get_joint_idx(joint)

    num_of_aux_var = get_num_of_aux_var(force)
    Apparatus(
        id,
        joint,
        force,
        num_of_aux_var,
        full_idx, 
        full_idx,
        ApparatusCache(joint,force,full_idx,)
    )
end

FloatingSphericalJoint(id,hen2egg,force=NoForce())  = proto_joint_apparatus(id,hen2egg,force,:FloatingSpherical)
OrbitalSphericalJoint(id,hen2egg,force=NoForce())   = proto_joint_apparatus(id,hen2egg,force,:OrbitalSpherical)
PlanarSphericalJoint(id,hen2egg,force=NoForce())    = proto_joint_apparatus(id,hen2egg,force,:PlanarSpherical)
PrismaticSphericalJoint(id,hen2egg,force=NoForce()) = proto_joint_apparatus(id,hen2egg,force,:PrismaticSpherical)
SphericalJoint(id,hen2egg,force=NoForce())          = proto_joint_apparatus(id,hen2egg,force,:Spherical)
FloatingUniversalJoint(id,hen2egg,force=NoForce())  = proto_joint_apparatus(id,hen2egg,force,:FloatingUniversal)
OrbitalUniversalJoint(id,hen2egg,force=NoForce())   = proto_joint_apparatus(id,hen2egg,force,:OrbitalUniversal)
PlanarUniversalJoint(id,hen2egg,force=NoForce())    = proto_joint_apparatus(id,hen2egg,force,:PlanarUniversal)
PrismaticUniversalJoint(id,hen2egg,force=NoForce()) = proto_joint_apparatus(id,hen2egg,force,:PrismaticUniversal)
UniversalPrismaticJoint(id,hen2egg,force=NoForce()) = proto_joint_apparatus(id,hen2egg,force,:UniversalPrismatic)
UniversalJoint(id,hen2egg,force=NoForce())          = proto_joint_apparatus(id,hen2egg,force,:Universal)
FloatingRevoluteJoint(id,hen2egg,force=NoForce())   = proto_joint_apparatus(id,hen2egg,force,:FloatingRevolute)
OrbitalRevoluteJoint(id,hen2egg,force=NoForce())    = proto_joint_apparatus(id,hen2egg,force,:OrbitalRevolute)
PlanarRevoluteJoint(id,hen2egg,force=NoForce())     = proto_joint_apparatus(id,hen2egg,force,:PlanarRevolute)
CylindricalJoint(id,hen2egg,force=NoForce())        = proto_joint_apparatus(id,hen2egg,force,:Cylindrical)
RevoluteJoint(id,hen2egg,force=NoForce())           = proto_joint_apparatus(id,hen2egg,force,:Revolute)
FloatingJoint(id,hen2egg,force=NoForce())           = proto_joint_apparatus(id,hen2egg,force,:Floating)
OrbitalJoint(id,hen2egg,force=NoForce())            = proto_joint_apparatus(id,hen2egg,force,:Orbital)
PlanarJoint(id,hen2egg,force=NoForce())             = proto_joint_apparatus(id,hen2egg,force,:Planar)
PrismaticJoint(id,hen2egg,force=NoForce())          = proto_joint_apparatus(id,hen2egg,force,:Prismatic)
FixedJoint(id,hen2egg,force=NoForce())              = proto_joint_apparatus(id,hen2egg,force,:Fixed)

const PinJoint = SphericalJoint

function stretch!(appar::Apparatus,c) end


function get_params!(params,appar::Apparatus,)
end


get_num_of_params(appar::Apparatus) = 0

get_num_of_dims(appar::Apparatus{<:AbstractJoint,<:AbstractForce}) = get_num_of_dims(appar.joint)
get_num_of_dims(appar::Apparatus{NoJoint,<:AbstractForce}) = get_num_of_dims(appar.force)

get_numbertype(appar::Apparatus{<:AbstractJoint,<:AbstractForce}) = get_numbertype(appar.joint)
get_numbertype(appar::Apparatus{NoJoint,<:AbstractForce}) = get_numbertype(appar.force)
