#done use the basic types of Constraints
#todo parameterization of joints
#note can full rotation constraints be linear?
"""
ID
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct ID{sigType,pidType,aidType}
    "Signifier of body"
    rbsig::sigType
    "Index of the anchor point"
    pid::pidType
    "Index of the translational axis"
    trlid::aidType
    "Index of the rotational axis"
    rotid::aidType
end

function ID(rbsig,pid,trlid)
    ID(rbsig,pid,trlid,trlid)
end

function ID(rbsig,pid)
    ID(rbsig,pid,1,1)
end

"""
Hen 2 Egg
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct Hen2Egg{henType<:ID,eggType<:ID}
    "id"
    id::Int
    "hen/parent/predecessor"
    hen::henType
    "egg/child/successor"
    egg::eggType
end

"""
空约束类。
$(TYPEDEF)
"""
struct EmptyConstraint{T} <: ExternalConstraints{T}
    nconstraints::Int64
    indices::Vector{Int64}
    values::T
end

get_numbertype(cst::ExternalConstraints{<:AbstractArray{T}}) where T = T

"""
空约束构造子。
$(TYPEDEF)
"""
function EmptyConstraint(values=Vector{Float64}())
    EmptyConstraint(0,Vector{Int64}(),values)
end

function make_constraints_function(::EmptyConstraint)
    inner_constraints_function(q) = Vector{eltype(q)}()
    inner_constraints_function(q,d) = Vector{eltype(q)}()
    inner_constraints_function
end

function make_constraints_jacobian(::EmptyConstraint)
    inner_constraints_jacobian(q) = Array{eltype(q)}(undef,0,length(q))
end

"""
刚体固定约束类。
$(TYPEDEF)
"""
struct FixedBodyConstraint{T} <: ExternalConstraints{T}
    nconstraints::Int64
    indices::Vector{Int64}
    values::T
end

"""
固定约束构造子。
$(TYPEDSIGNATURES)
"""
function FixedBodyConstraint(rbs,mem2sysfull,bodyid)
    body = rbs[bodyid]
    nmcs = body.state.cache.funcs.nmcs
    q_rb = body.state.coords.q
    pres_idx = find_full_pres_indices(nmcs,q_rb)
    indices = body2q[bodyid][pres_idx]
    values = q_rb[pres_idx]
    FixedBodyConstraint(length(indices),indices,values)
end

function make_constraints_function(cst::FixedBodyConstraint)
    (;indices, values) = cst
    @inline @inbounds inner_constraints_function(q)   = q[indices]-values
    @inline @inbounds inner_constraints_function(q,d) = q[indices]-d
    inner_constraints_function
end

function make_constraints_jacobian(cst::FixedBodyConstraint)
    num_of_constraints = cst.nconstraints
    indices = cst.indices
    @inline @inbounds function inner_constraints_jacobian(q)
        nq = length(q)
        ret = zeros(eltype(q),num_of_constraints,nq)
        for (iΦ,i) in enumerate(indices)
            ret[iΦ,i] = 1
        end
        ret
    end
end

"""
刚体坐标固定约束类，适用于单个坐标。
$(TYPEDEF)
"""
struct FixedIndicesConstraint{T} <: ExternalConstraints{T}
    id::Int64
    nconstraints::Int64
    indices::Vector{Int64}
    values::T
end

"""
刚体坐标固定约束构造子。
$(TYPEDSIGNATURES)
"""
function FixedIndicesConstraint(id,indices,values)
    FixedIndicesConstraint(id,length(indices),indices,values)
end

function make_constraints_function(cst::FixedIndicesConstraint,st)
    @unpack indices, values = cst
    @inline @inbounds inner_constraints_function(q)   = q[indices]-values
    @inline @inbounds inner_constraints_function(q,d) = q[indices]-d
    inner_constraints_function
end

function make_constraints_jacobian(cst::FixedIndicesConstraint,st)
    (;indexed,numbered) = st.connectivity
    num_of_constraints = cst.nconstraints
    indices = cst.indices
    (;sysfree,nfree) = indexed
    @inline @inbounds function inner_constraints_jacobian(q)
        nq = length(q)
        ret = zeros(eltype(q),num_of_constraints,nq)
        for (iΦ,i) in enumerate(indices)
            ret[iΦ,i] = 1
        end
        ret[:,sysfree]
    end
end



"""
刚体通用线性约束类。
$(TYPEDEF)
"""
struct LinearJoint{valueType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::Vector{valueType}
    A::Matrix{valueType}
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function LinearJoint(A,values)
    num_of_constraints = size(A,1)
    LinearJoint(id,num_of_constraints,values,A)
end

function make_constraints_function(cst::LinearJoint,indexed,numbered)
    (;mem2sysfull) = indexed
    (;nconstraints,values,A) = cst
    function _inner_constraints_function(q,d)
        ret = zeros(eltype(q),nconstraints)
        ret .= A*q .- d
        ret
    end
    inner_constraints_function(q)   = _inner_constraints_function(q,values)
    inner_constraints_function(q,d) = _inner_constraints_function(q,d)
    inner_constraints_function
end

function make_constraints_jacobian(cst::LinearJoint,indexed,numbered)
    (;sysfree,nfree) = indexed
    (;nconstraints,values,A) = cst
    function _inner_constraints_jacobian(q,c)
        q̌ = @view q[sysfree]
        ret = zeros(eltype(q̌),nconstraints,nfree)
        ret .= A
        ret
    end
    inner_constraints_jacobian(q)   = _inner_constraints_jacobian(q,0)
    inner_constraints_jacobian(q,c) = _inner_constraints_jacobian(q,c)
    inner_constraints_jacobian
end

"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct PrototypeJoint{valueType,hen2eggType,axesType,maskType} <: ExternalConstraints{valueType}
    id::Int
    hen2egg::hen2eggType
    nconstraints::Int
    ndof::Int
    axes_trl_hen::axesType
    axes_trl_egg::axesType
    axes_rot_hen::axesType
    axes_rot_egg::axesType
    mask_1st::maskType
    mask_2nd::maskType
    mask_3rd_hen::maskType
    mask_3rd_egg::maskType
    mask_4th::maskType
    values::valueType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function PrototypeJoint(id,hen2egg,joint_type::Symbol) 
    (;hen,egg) = hen2egg
    joint_info = get_joint_info(joint_type)
    (;ntrl, nrot, ndof, ncsts, mask_1st, mask_2nd, mask_3rd_hen, mask_3rd_egg, mask_4th) = joint_info
    T = get_numbertype(hen.rbsig)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    q_hen = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.origin_position,state_hen.R)
    q_egg = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.origin_position,state_egg.R)
    R_hen = NCF.find_rotation(nmcs_hen,q_hen)
    R_egg = NCF.find_rotation(nmcs_egg,q_egg)
    # translate     
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    r_hen2egg = C_egg*q_egg .- C_hen*q_hen
    # translate on hen
    axes_trl_hen = Axes(inv(R_hen)*get_orthonormal_axes(R_hen*hen.rbsig.prop.loci[hen.rotid].axes.normal)) 
    trl_hen_n = R_hen*axes_trl_hen.normal
    trl_hen_t = R_hen*axes_trl_hen.tangent
    trl_hen_b = R_hen*axes_trl_hen.bitangent
    # translate on egg
    axes_trl_egg = Axes(inv(R_egg)*get_orthonormal_axes(R_egg*egg.rbsig.prop.loci[egg.trlid].axes.normal)) 
    trl_egg_n = R_egg*axes_trl_egg.normal
    trl_egg_t = R_egg*axes_trl_egg.tangent
    trl_egg_b = R_egg*axes_trl_egg.bitangent
    # rotate of egg
    axes_rot_glb = get_orthonormal_axes(R_egg*egg.rbsig.prop.loci[egg.rotid].axes.normal)
    axes_rot_egg = Axes(inv(R_egg)*axes_rot_glb)
    axes_rot_hen = Axes(inv(R_hen)*axes_rot_glb)
    rot_egg_n = R_egg*axes_rot_egg.normal
    rot_egg_t = R_egg*axes_rot_egg.tangent
    rot_egg_b = R_egg*axes_rot_egg.bitangent
    rot_hen_n = R_hen*axes_rot_hen.normal
    rot_hen_t = R_hen*axes_rot_hen.tangent
    rot_hen_b = R_hen*axes_rot_hen.bitangent
    # first    
    Φ_1st = r_hen2egg
    # second
    Φ_2nd = [
        rot_hen_t'*rot_egg_b,
        rot_hen_t'*rot_egg_n,
        rot_hen_b'*rot_egg_n
    ]
    # third
    Φ_3rd_hen = [
        trl_hen_n'*r_hen2egg,
        trl_hen_t'*r_hen2egg,
        trl_hen_b'*r_hen2egg,
    ]
    Φ_3rd_egg = [
        trl_egg_n'*r_hen2egg,
        trl_egg_t'*r_hen2egg,
        trl_egg_b'*r_hen2egg,
    ]
    # fourth    
    Φ_4th = [r_hen2egg'*r_hen2egg]
    # constraint values
    values = vcat(
        Φ_4th[mask_4th],
        Φ_3rd_hen[mask_3rd_hen], # translate prior to 
        Φ_3rd_egg[mask_3rd_egg], 
        Φ_1st[mask_1st],
        Φ_2nd[mask_2nd], # rotate
    )
    # @show joint_info
    # @show values
    PrototypeJoint(
        id,hen2egg,
        ncsts,ndof,
        axes_trl_hen,axes_trl_egg,axes_rot_hen,axes_rot_egg,
        mask_1st, mask_2nd, mask_3rd_hen, mask_3rd_egg, mask_4th,
        values
    )
end

function get_joint_info(joint_type::Symbol)
    (joint_type == :FloatingSpherical)  && (return (ntrl = 3, nrot = 3, ndof = 6, ncsts = 0, mask_1st = Int[],   mask_2nd = Int[]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] ))
    (joint_type == :OrbitalSpherical)   && (return (ntrl = 2, nrot = 3, ndof = 5, ncsts = 1, mask_1st = Int[],   mask_2nd = Int[]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1]   ))
    (joint_type == :PlanarSpherical)    && (return (ntrl = 2, nrot = 3, ndof = 5, ncsts = 1, mask_1st = Int[],   mask_2nd = Int[]  , mask_3rd_hen = [1],   mask_3rd_egg = Int[], mask_4th = Int[] ))
    (joint_type == :PrismaticSpherical) && (return (ntrl = 1, nrot = 3, ndof = 4, ncsts = 2, mask_1st = Int[],   mask_2nd = Int[]  , mask_3rd_hen = [2,3], mask_3rd_egg = Int[], mask_4th = Int[] ))
    (joint_type == :Spherical)          && (return (ntrl = 0, nrot = 3, ndof = 3, ncsts = 3, mask_1st = [1,2,3], mask_2nd = Int[]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] ))
    (joint_type == :FloatingUniversal)  && (return (ntrl = 3, nrot = 2, ndof = 5, ncsts = 1, mask_1st = Int[],   mask_2nd = [1]    , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2
    (joint_type == :OrbitalUniversal)   && (return (ntrl = 2, nrot = 2, ndof = 4, ncsts = 2, mask_1st = Int[],   mask_2nd = [1]    , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1]   )) #t1'*t2
    (joint_type == :PlanarUniversal)    && (return (ntrl = 2, nrot = 2, ndof = 4, ncsts = 2, mask_1st = Int[],   mask_2nd = [1]    , mask_3rd_hen = [1],   mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2
    (joint_type == :PrismaticUniversal) && (return (ntrl = 1, nrot = 2, ndof = 3, ncsts = 3, mask_1st = Int[],   mask_2nd = [1]    , mask_3rd_hen = [2,3], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2
    (joint_type == :UniversalPrismatic) && (return (ntrl = 1, nrot = 2, ndof = 3, ncsts = 3, mask_1st = Int[],   mask_2nd = [1]    , mask_3rd_hen = Int[], mask_3rd_egg = [2,3], mask_4th = Int[] )) #t1'*t2
    (joint_type == :Universal)          && (return (ntrl = 0, nrot = 2, ndof = 2, ncsts = 4, mask_1st = [1,2,3], mask_2nd = [1]    , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2
    (joint_type == :FloatingRevolute)   && (return (ntrl = 3, nrot = 1, ndof = 4, ncsts = 2, mask_1st = Int[],   mask_2nd = [2,3]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*n, t2'*n
    (joint_type == :OrbitalRevolute)    && (return (ntrl = 2, nrot = 1, ndof = 3, ncsts = 3, mask_1st = Int[],   mask_2nd = [2,3]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1]   )) #t1'*n, t2'*n
    (joint_type == :PlanarRevolute)     && (return (ntrl = 2, nrot = 1, ndof = 3, ncsts = 3, mask_1st = Int[],   mask_2nd = [2,3]  , mask_3rd_hen = [1],   mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*n, t2'*n
    (joint_type == :Cylindrical)        && (return (ntrl = 1, nrot = 1, ndof = 2, ncsts = 4, mask_1st = Int[],   mask_2nd = [2,3]  , mask_3rd_hen = [2,3], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*n, t2'*n
    (joint_type == :Revolute)           && (return (ntrl = 0, nrot = 1, ndof = 1, ncsts = 5, mask_1st = [1,2,3], mask_2nd = [2,3]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*n, t2'*n
    (joint_type == :Floating)           && (return (ntrl = 3, nrot = 0, ndof = 3, ncsts = 3, mask_1st = Int[],   mask_2nd = [1,2,3], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Orbital)            && (return (ntrl = 2, nrot = 0, ndof = 2, ncsts = 4, mask_1st = Int[],   mask_2nd = [1,2,3], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1]   )) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Planar)             && (return (ntrl = 2, nrot = 0, ndof = 2, ncsts = 4, mask_1st = Int[],   mask_2nd = [1,2,3], mask_3rd_hen = [1],   mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Prismatic)          && (return (ntrl = 1, nrot = 0, ndof = 1, ncsts = 5, mask_1st = Int[],   mask_2nd = [1,2,3], mask_3rd_hen = [2,3], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Fixed)              && (return (ntrl = 0, nrot = 0, ndof = 0, ncsts = 6, mask_1st = [1,2,3], mask_2nd = [1,2,3], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2, t1'*n, t2'*n
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