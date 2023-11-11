#done use the basic types of Constraints
#todo parameterization of joints
#note can full rotation cstr be linear?
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
struct EmptyConstraint{T} <: ExtrinsicConstraints{T}
    num_of_cstr::Int64
    idx::Vector{Int64}
    values::T
end

get_numbertype(cst::ExtrinsicConstraints{<:AbstractArray{T}}) where T = T

"""
空约束构造子。
$(TYPEDEF)
"""
function EmptyConstraint(values=Vector{Float64}())
    EmptyConstraint(0,Vector{Int64}(),values)
end

function make_cstr_function(::EmptyConstraint)
    inner_cstr_function(q) = Vector{eltype(q)}()
    inner_cstr_function(q,d) = Vector{eltype(q)}()
    inner_cstr_function
end

function make_cstr_jacobian(::EmptyConstraint)
    inner_cstr_jacobian(q) = Array{eltype(q)}(undef,0,length(q))
end

"""
刚体固定约束类。
$(TYPEDEF)
"""
struct FixedBodyConstraint{T} <: ExtrinsicConstraints{T}
    num_of_cstr::Int64
    idx::Vector{Int64}
    values::T
end

"""
固定约束构造子。
$(TYPEDSIGNATURES)
"""
function FixedBodyConstraint(rbs,bodyid2sys_full_coords,bodyid)
    body = rbs[bodyid]
    nmcs = body.coords.nmcs
    q_rb = body.state.coords.q
    pres_idx = find_full_pres_idx(nmcs,q_rb)
    idx = bodyid2q[bodyid][pres_idx]
    values = q_rb[pres_idx]
    FixedBodyConstraint(length(idx),idx,values)
end

function make_cstr_function(cst::FixedBodyConstraint)
    (;idx, values) = cst
    @inline @inbounds inner_cstr_function(q)   = q[idx]-values
    @inline @inbounds inner_cstr_function(q,d) = q[idx]-d
    inner_cstr_function
end

function make_cstr_jacobian(cst::FixedBodyConstraint)
    num_of_cstr = cst.num_of_cstr
    idx = cst.idx
    @inline @inbounds function inner_cstr_jacobian(q)
        nq = length(q)
        ret = zeros(eltype(q),num_of_cstr,nq)
        for (iΦ,i) in enumerate(idx)
            ret[iΦ,i] = 1
        end
        ret
    end
end

"""
刚体坐标固定约束类，适用于单个坐标。
$(TYPEDEF)
"""
struct FixedIndicesConstraint{T} <: ExtrinsicConstraints{T}
    id::Int64
    num_of_cstr::Int64
    idx::Vector{Int64}
    values::T
end

"""
刚体坐标固定约束构造子。
$(TYPEDSIGNATURES)
"""
function FixedIndicesConstraint(id,idx,values)
    FixedIndicesConstraint(id,length(idx),idx,values)
end

function make_cstr_function(cst::FixedIndicesConstraint,st)
    (;idx, values) = cst
    @inline @inbounds inner_cstr_function(q)   = q[idx]-values
    @inline @inbounds inner_cstr_function(q,d) = q[idx]-d
    inner_cstr_function
end

function make_cstr_jacobian(cst::FixedIndicesConstraint,st)
    (;indexed,numbered) = st.connectivity
    num_of_cstr = cst.num_of_cstr
    idx = cst.idx
    (;sys_free_idx,num_of_free_coords) = indexed
    @inline @inbounds function inner_cstr_jacobian(q)
        nq = length(q)
        ret = zeros(eltype(q),num_of_cstr,nq)
        for (iΦ,i) in enumerate(idx)
            ret[iΦ,i] = 1
        end
        ret[:,sys_free_idx]
    end
end



"""
刚体通用线性约束类。
$(TYPEDEF)
"""
struct LinearJoint{valueType} <: ExtrinsicConstraints{valueType}
    id::Int
    num_of_cstr::Int
    values::Vector{valueType}
    A::Matrix{valueType}
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function LinearJoint(A,values)
    num_of_cstr = size(A,1)
    LinearJoint(id,num_of_cstr,values,A)
end

function make_cstr_function(cst::LinearJoint,indexed,numbered)
    (;bodyid2sys_full_coords) = indexed
    (;num_of_cstr,values,A) = cst
    function _inner_cstr_function(q,d)
        ret = zeros(eltype(q),num_of_cstr)
        ret .= A*q .- d
        ret
    end
    inner_cstr_function(q)   = _inner_cstr_function(q,values)
    inner_cstr_function(q,d) = _inner_cstr_function(q,d)
    inner_cstr_function
end

function make_cstr_jacobian(cst::LinearJoint,indexed,numbered)
    (;sys_free_idx,num_of_free_coords) = indexed
    (;num_of_cstr,values,A) = cst
    function _inner_cstr_jacobian(q,c)
        q̌ = @view q[sys_free_idx]
        ret = zeros(eltype(q̌),num_of_cstr,num_of_free_coords)
        ret .= A
        ret
    end
    inner_cstr_jacobian(q)   = _inner_cstr_jacobian(q,0)
    inner_cstr_jacobian(q,c) = _inner_cstr_jacobian(q,c)
    inner_cstr_jacobian
end

function get_joint_info(joint_type::Symbol)
    (joint_type == :FloatingSpherical)  && (return (ntrl = 3, nrot = 3, num_of_dof = 6, num_of_cstr = 0, mask_1st = Int[],   mask_2nd = Int[]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] ))
    (joint_type == :OrbitalSpherical)   && (return (ntrl = 2, nrot = 3, num_of_dof = 5, num_of_cstr = 1, mask_1st = Int[],   mask_2nd = Int[]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1]   ))
    (joint_type == :PlanarSpherical)    && (return (ntrl = 2, nrot = 3, num_of_dof = 5, num_of_cstr = 1, mask_1st = Int[],   mask_2nd = Int[]  , mask_3rd_hen = [1],   mask_3rd_egg = Int[], mask_4th = Int[] ))
    (joint_type == :PrismaticSpherical) && (return (ntrl = 1, nrot = 3, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[],   mask_2nd = Int[]  , mask_3rd_hen = [2,3], mask_3rd_egg = Int[], mask_4th = Int[] ))
    (joint_type == :Spherical)          && (return (ntrl = 0, nrot = 3, num_of_dof = 3, num_of_cstr = 3, mask_1st = [1,2,3], mask_2nd = Int[]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] ))
    (joint_type == :FloatingUniversal)  && (return (ntrl = 3, nrot = 2, num_of_dof = 5, num_of_cstr = 1, mask_1st = Int[],   mask_2nd = [1]    , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2
    (joint_type == :OrbitalUniversal)   && (return (ntrl = 2, nrot = 2, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[],   mask_2nd = [1]    , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1]   )) #t1'*t2
    (joint_type == :PlanarUniversal)    && (return (ntrl = 2, nrot = 2, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[],   mask_2nd = [1]    , mask_3rd_hen = [1],   mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2
    (joint_type == :PrismaticUniversal) && (return (ntrl = 1, nrot = 2, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[],   mask_2nd = [1]    , mask_3rd_hen = [2,3], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2
    (joint_type == :UniversalPrismatic) && (return (ntrl = 1, nrot = 2, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[],   mask_2nd = [1]    , mask_3rd_hen = Int[], mask_3rd_egg = [2,3], mask_4th = Int[] )) #t1'*t2
    (joint_type == :Universal)          && (return (ntrl = 0, nrot = 2, num_of_dof = 2, num_of_cstr = 4, mask_1st = [1,2,3], mask_2nd = [1]    , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2
    (joint_type == :FloatingRevolute)   && (return (ntrl = 3, nrot = 1, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[],   mask_2nd = [2,3]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*n, t2'*n
    (joint_type == :OrbitalRevolute)    && (return (ntrl = 2, nrot = 1, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[],   mask_2nd = [2,3]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1]   )) #t1'*n, t2'*n
    (joint_type == :PlanarRevolute)     && (return (ntrl = 2, nrot = 1, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[],   mask_2nd = [2,3]  , mask_3rd_hen = [1],   mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*n, t2'*n
    (joint_type == :Cylindrical)        && (return (ntrl = 1, nrot = 1, num_of_dof = 2, num_of_cstr = 4, mask_1st = Int[],   mask_2nd = [2,3]  , mask_3rd_hen = [2,3], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*n, t2'*n
    (joint_type == :Revolute)           && (return (ntrl = 0, nrot = 1, num_of_dof = 1, num_of_cstr = 5, mask_1st = [1,2,3], mask_2nd = [2,3]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*n, t2'*n
    (joint_type == :Floating)           && (return (ntrl = 3, nrot = 0, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[],   mask_2nd = [1,2,3], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Orbital)            && (return (ntrl = 2, nrot = 0, num_of_dof = 2, num_of_cstr = 4, mask_1st = Int[],   mask_2nd = [1,2,3], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1]   )) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Planar)             && (return (ntrl = 2, nrot = 0, num_of_dof = 2, num_of_cstr = 4, mask_1st = Int[],   mask_2nd = [1,2,3], mask_3rd_hen = [1],   mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Prismatic)          && (return (ntrl = 1, nrot = 0, num_of_dof = 1, num_of_cstr = 5, mask_1st = Int[],   mask_2nd = [1,2,3], mask_3rd_hen = [2,3], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Fixed)              && (return (ntrl = 0, nrot = 0, num_of_dof = 0, num_of_cstr = 6, mask_1st = [1,2,3], mask_2nd = [1,2,3], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[] )) #t1'*t2, t1'*n, t2'*n
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
刚体铰接约束类。
$(TYPEDEF)
"""
struct PrototypeJoint{hen2eggType,axesType,maskType,halfType,hessType,valueType} <: ExtrinsicConstraints{valueType}
    id::Int
    hen2egg::hen2eggType
    num_of_cstr::Int
    num_of_dof::Int
    axes_trl_hen::axesType
    axes_trl_egg::axesType
    axes_rot_hen::axesType
    axes_rot_egg::axesType
    mask_1st::maskType
    mask_2nd::maskType
    mask_3rd::maskType
    mask_4th::maskType
    transformations::halfType
    halves::Vector{halfType}
    hessians::Vector{hessType}
    violations::valueType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function PrototypeJoint(id,hen2egg,joint_type::Symbol) 
    (;hen,egg) = hen2egg
    joint_info = get_joint_info(joint_type)
    (;  ntrl, nrot, 
        num_of_dof, num_of_cstr, 
        mask_1st, mask_2nd, 
        mask_3rd_hen, 
        mask_3rd_egg, 
        mask_4th
    ) = joint_info
    mask_3rd = vcat(mask_3rd_hen,mask_3rd_egg.+3)
    # bits_1st = in.([1,2,3],Ref(mask_1st))
    # bits_2nd = in.([1,2,3],Ref(mask_2nd))
    # bits_3rd = in.([1,2,3,4,5,6],Ref(mask_3rd))
    # bits_4th = in.([1,],Ref(mask_4th))
    T = get_numbertype(hen.rbsig)
    nmcs_hen = hen.rbsig.coords.nmcs
    nmcs_egg = egg.rbsig.coords.nmcs
    num_of_coords_hen = get_num_of_coords(nmcs_hen)
    num_of_coords_egg = get_num_of_coords(nmcs_egg)
    num_of_jointed_coords = num_of_coords_hen + num_of_coords_egg
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    q_hen = cartesian_frame2coords(nmcs_hen,state_hen.origin_position,state_hen.R)
    q_egg = cartesian_frame2coords(nmcs_egg,state_egg.origin_position,state_egg.R)
    X_hen = NCF.get_X(nmcs_hen,q_hen)
    X_egg = NCF.get_X(nmcs_egg,q_egg)
    invX̄_hen = nmcs_hen.data.invX̄
    invX̄_egg = nmcs_egg.data.invX̄
    # translate
    c_hen = to_local_coords(nmcs_hen,hen.rbsig.prop.loci[hen.pid].position)
    C_hen = to_transformation(nmcs_hen,c_hen)
    c_egg = to_local_coords(nmcs_egg,egg.rbsig.prop.loci[egg.pid].position)
    C_egg = to_transformation(nmcs_egg,c_egg)
    J = [-C_hen C_egg;] |> sparse 
    q = vcat(q_hen,q_egg)
    axes_trl_hen = invX̄_hen*hen.rbsig.prop.loci[hen.rotid].axes
    axes_trl_egg = invX̄_egg*egg.rbsig.prop.loci[egg.trlid].axes
    axes_rot_egg = invX̄_egg*egg.rbsig.prop.loci[egg.rotid].axes
    select_uvw_hen = BlockDiagonal([nmcs_hen.conversion_to_X,zero(nmcs_egg.conversion_to_X)])
    select_uvw_egg = BlockDiagonal([zero(nmcs_hen.conversion_to_X),nmcs_egg.conversion_to_X])
    I3_Bool = I(3)
    # translate
    heros = spzeros(T,num_of_jointed_coords,num_of_jointed_coords)
    # hes 1st
    half_1st = fill(heros,3)
    # hes 4th
    half_4th = [(J'*J) |> sparse]
    # hes 3rd
    half_3rd = fill(heros,6)
    half_2nd = fill(heros,3)
    if (nmcs_hen isa NCF.NC3D12C) && (nmcs_egg isa NCF.NC3D12C)
        axes_rot_hen = inv(X_hen)*X_egg*axes_rot_egg
        for i = 1:3
            axis_hen = axes_trl_hen.X[:,i]
            axis_zero = zero(axis_hen)
            half_3rd[i] = select_uvw_hen'*kron(vcat(0,axis_hen,0,axis_zero),I3_Bool)*J |> sparse
        end
        # hes 3rd on egg
        # translate on egg
        for i = 1:3
            axis_egg = axes_trl_egg.X[:,i]
            axis_zero = zero(axis_egg)
            half_3rd[3+i] = select_uvw_egg'*kron(vcat(0,axis_zero,0,axis_egg),I3_Bool)*J |> sparse
        end
        # hes 2nd
        # rotate of egg
        axes_idx = [
            (2,3),
            (2,1),
            (3,1)
        ]
        for (i,(id_axis_hen,id_axis_egg)) in enumerate(axes_idx)
            axis_hen = axes_rot_hen.X[:,id_axis_hen]
            axis_egg = axes_rot_egg.X[:,id_axis_egg]
            axis_zero = zero(axis_egg)
            half_2nd[i] = select_uvw_hen'*
                kron(vcat(0,axis_hen,0,axis_zero),I3_Bool)*
                kron(vcat(0,axis_zero,0,axis_egg),I3_Bool)'*
                select_uvw_egg |> sparse
        end
    else
        # not to be used
        axes_rot_hen = axes_rot_egg
    end
    hess_1st = [(H .+ H') |> Symmetric for H in half_1st]
    hess_4th = [(H .+ H') |> Symmetric for H in half_4th]
    hess_3rd = [(H .+ H') |> Symmetric for H in half_3rd]
    hess_2nd = [(H .+ H') |> Symmetric for H in half_2nd]
    # cstr cache
    halves = vcat(
        half_1st[mask_1st],
        half_4th[mask_4th],
        half_3rd[mask_3rd],
        half_2nd[mask_2nd];
    )
    hessians = vcat(
        hess_1st[mask_1st],
        hess_4th[mask_4th],
        hess_3rd[mask_3rd],
        hess_2nd[mask_2nd];
    )
    # cstr values
    Refq = Ref(q)
    RefqT = Ref(q')
    transformations = J[mask_1st,:]
    values = RefqT.*halves.*Refq
    values[mask_1st] .+= transformations*q
    # @show joint_info
    # @show values
    PrototypeJoint(
        id,hen2egg,
        num_of_cstr,
        num_of_dof,
        axes_trl_hen,
        axes_trl_egg,
        axes_rot_hen,
        axes_rot_egg,
        mask_1st,
        mask_2nd,
        mask_3rd,
        mask_4th,
        transformations,
        halves,
        hessians,
        values
    )
end
