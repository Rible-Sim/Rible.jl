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
    "No. of the anchor point"
    pid::pidType
    "No. of the translational axis"
    trlid::aidType
    "No. of the rotational axis"
    rotid::aidType
end

function ID(rbsig,pid,trlid)
    ID(rbsig,pid,trlid,trlid)
end

function ID(rbsig,pid)
    ID(rbsig,pid,1,1)
end

"""
点对点类。
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct End2End{henType<:ID,eggType<:ID}
    "编号"
    id::Int
    "起始点"
    hen::henType
    "终止点"
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

function make_Φ(::EmptyConstraint)
    inner_Φ(q) = Vector{eltype(q)}()
    inner_Φ(q,d) = Vector{eltype(q)}()
    inner_Φ
end

function make_A(::EmptyConstraint)
    inner_A(q) = Array{eltype(q)}(undef,0,length(q))
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
function FixedBodyConstraint(rbs,mem2sysfull,rbid)
    rb = rbs[rbid]
    lncs = rb.state.cache.funcs.lncs
    q_rb = rb.state.coords.q
    pres_idx = find_full_pres_indices(lncs,q_rb)
    indices = body2q[rbid][pres_idx]
    values = q_rb[pres_idx]
    FixedBodyConstraint(length(indices),indices,values)
end

function make_Φ(cst::FixedBodyConstraint)
    (;indices, values) = cst
    @inline @inbounds inner_Φ(q)   = q[indices]-values
    @inline @inbounds inner_Φ(q,d) = q[indices]-d
    inner_Φ
end

function make_A(cst::FixedBodyConstraint)
    nΦ = cst.nconstraints
    indices = cst.indices
    @inline @inbounds function inner_A(q)
        nq = length(q)
        ret = zeros(eltype(q),nΦ,nq)
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

function make_Φ(cst::FixedIndicesConstraint,st)
    @unpack indices, values = cst
    @inline @inbounds inner_Φ(q)   = q[indices]-values
    @inline @inbounds inner_Φ(q,d) = q[indices]-d
    inner_Φ
end

function make_A(cst::FixedIndicesConstraint,st)
    (;indexed,numbered) = st.connectivity
    nΦ = cst.nconstraints
    indices = cst.indices
    (;sysfree,nfree) = indexed
    @inline @inbounds function inner_A(q)
        nq = length(q)
        ret = zeros(eltype(q),nΦ,nq)
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
    nΦ = size(A,1)
    LinearJoint(id,nΦ,values,A)
end

function make_Φ(cst::LinearJoint,indexed,numbered)
    (;mem2sysfull) = indexed
    (;nconstraints,values,A) = cst
    function _inner_Φ(q,d)
        ret = zeros(eltype(q),nconstraints)
        ret .= A*q .- d
        ret
    end
    inner_Φ(q)   = _inner_Φ(q,values)
    inner_Φ(q,d) = _inner_Φ(q,d)
    inner_Φ
end

function make_A(cst::LinearJoint,indexed,numbered)
    (;sysfree,nfree) = indexed
    (;nconstraints,values,A) = cst
    function _inner_A(q,c)
        q̌ = @view q[sysfree]
        ret = zeros(eltype(q̌),nconstraints,nfree)
        ret .= A
        ret
    end
    inner_A(q)   = _inner_A(q,0)
    inner_A(q,c) = _inner_A(q,c)
    inner_A
end

"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct PrototypeJoint{valueType,e2eType,axesType,maskType} <: ExternalConstraints{valueType}
    id::Int
    e2e::e2eType
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
function PrototypeJoint(id,e2e,joint_type::Symbol) 
    (;hen,egg) = e2e
    joint_info = get_joint_info(joint_type)
    (;ntrl, nrot, ndof, ncsts, mask_1st, mask_2nd, mask_3rd_hen, mask_3rd_egg, mask_4th) = joint_info
    T = get_numbertype(hen.rbsig)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    q_hen = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R)
    q_egg = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R)
    R_hen = NCF.find_R(nmcs_hen,q_hen)
    R_egg = NCF.find_R(nmcs_egg,q_egg)
    # translate     
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    r_hen2egg = C_egg*q_egg .- C_hen*q_hen
    # translate on hen
    axes_trl_hen = Axes(inv(R_hen)*get_orthonormal_axes(R_hen*hen.rbsig.prop.ās[hen.rotid])) 
    trl_hen_n = R_hen*axes_trl_hen.n
    trl_hen_t1 = R_hen*axes_trl_hen.t1
    trl_hen_t2 = R_hen*axes_trl_hen.t2
    # translate on egg
    axes_trl_egg = Axes(inv(R_egg)*get_orthonormal_axes(R_egg*egg.rbsig.prop.ās[egg.trlid])) 
    trl_egg_n = R_egg*axes_trl_egg.n
    trl_egg_t1 = R_egg*axes_trl_egg.t1
    trl_egg_t2 = R_egg*axes_trl_egg.t2
    # rotate of egg
    axes_rot_glb = get_orthonormal_axes(R_egg*egg.rbsig.prop.ās[egg.rotid])
    axes_rot_egg = Axes(inv(R_egg)*axes_rot_glb)
    axes_rot_hen = Axes(inv(R_hen)*axes_rot_glb)
    rot_egg_n = R_egg*axes_rot_egg.n
    rot_egg_t1 = R_egg*axes_rot_egg.t1
    rot_egg_t2 = R_egg*axes_rot_egg.t2
    rot_hen_n = R_hen*axes_rot_hen.n
    rot_hen_t1 = R_hen*axes_rot_hen.t1
    rot_hen_t2 = R_hen*axes_rot_hen.t2
    # first    
    Φ_1st = r_hen2egg
    # second
    Φ_2nd = [
        rot_hen_t1'*rot_egg_t2,
        rot_hen_t1'*rot_egg_n,
        rot_hen_t2'*rot_egg_n
    ]
    # third
    Φ_3rd_hen = [
        trl_hen_n'*r_hen2egg,
        trl_hen_t1'*r_hen2egg,
        trl_hen_t2'*r_hen2egg,
    ]
    Φ_3rd_egg = [
        trl_egg_n'*r_hen2egg,
        trl_egg_t1'*r_hen2egg,
        trl_egg_t2'*r_hen2egg,
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
        id,e2e,
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

FloatingSphericalJoint(id,e2e)  = PrototypeJoint(id,e2e,:FloatingSpherical)
OrbitalSphericalJoint(id,e2e)   = PrototypeJoint(id,e2e,:OrbitalSpherical)
PlanarSphericalJoint(id,e2e)    = PrototypeJoint(id,e2e,:PlanarSpherical)
PrismaticSphericalJoint(id,e2e) = PrototypeJoint(id,e2e,:PrismaticSpherical)
SphericalJoint(id,e2e)          = PrototypeJoint(id,e2e,:Spherical)
FloatingUniversalJoint(id,e2e)  = PrototypeJoint(id,e2e,:FloatingUniversal)
OrbitalUniversalJoint(id,e2e)   = PrototypeJoint(id,e2e,:OrbitalUniversal)
PlanarUniversalJoint(id,e2e)    = PrototypeJoint(id,e2e,:PlanarUniversal)
PrismaticUniversalJoint(id,e2e) = PrototypeJoint(id,e2e,:PrismaticUniversal)
UniversalPrismaticJoint(id,e2e) = PrototypeJoint(id,e2e,:UniversalPrismatic)
UniversalJoint(id,e2e)          = PrototypeJoint(id,e2e,:Universal)
FloatingRevoluteJoint(id,e2e)   = PrototypeJoint(id,e2e,:FloatingRevolute)
OrbitalRevoluteJoint(id,e2e)    = PrototypeJoint(id,e2e,:OrbitalRevolute)
PlanarRevoluteJoint(id,e2e)     = PrototypeJoint(id,e2e,:PlanarRevolute)
CylindricalJoint(id,e2e)        = PrototypeJoint(id,e2e,:Cylindrical)
RevoluteJoint(id,e2e)           = PrototypeJoint(id,e2e,:Revolute)
FloatingJoint(id,e2e)           = PrototypeJoint(id,e2e,:Floating)
OrbitalJoint(id,e2e)            = PrototypeJoint(id,e2e,:Orbital)
PlanarJoint(id,e2e)             = PrototypeJoint(id,e2e,:Planar)
PrismaticJoint(id,e2e)          = PrototypeJoint(id,e2e,:Prismatic)
FixedJoint(id,e2e)              = PrototypeJoint(id,e2e,:Fixed)

const PinJoint = SphericalJoint