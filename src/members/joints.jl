#done use the basic types of Constraints
#todo parameterization of joints
#note can full rotation cstr be linear?

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
    id::Int
    num_of_cstr::Int64
    idx::Vector{Int64}
    values::T
end

"""
固定约束构造子。
$(TYPEDSIGNATURES)
"""
function FixedBodyConstraint(id::Int,body::AbstractBody,indexed)
    (;bodyid2sys_full_coords) = indexed
    bodyid = body.prop.id
    nmcs = body.coords.nmcs
    state = body.state
    q = cartesian_frame2coords(
        nmcs,
        state.origin_position,state.R
    )
    pres_idx = find_independent_idx(nmcs,q)
    idx = bodyid2sys_full_coords[bodyid][pres_idx]
    values = q[pres_idx]
    FixedBodyConstraint(id,length(idx),idx,values)
end

function make_cstr_function(cst::FixedBodyConstraint,st)
    (;idx, values) = cst
    @inline @inbounds inner_cstr_function(q)   = q[idx]-values
    @inline @inbounds inner_cstr_function(q,d) = q[idx]-d
    inner_cstr_function
end

function make_cstr_jacobian(cst::FixedBodyConstraint,st)
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
    (joint_type == :FloatingSpherical)  && (return (ntrl = 3, nrot = 3, num_of_dof = 6, num_of_cstr = 0, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[]  ))
    (joint_type == :OrbitalSpherical)   && (return (ntrl = 2, nrot = 3, num_of_dof = 5, num_of_cstr = 1, mask_1st = Int[]  , mask_2nd = [1]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[]  ))
    (joint_type == :PlanarSpherical)    && (return (ntrl = 2, nrot = 3, num_of_dof = 5, num_of_cstr = 1, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [1]  , mask_3rd_egg = Int[], mask_4th = Int[]  ))
    (joint_type == :PrismaticSpherical) && (return (ntrl = 1, nrot = 3, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [2,3], mask_3rd_egg = Int[], mask_4th = Int[]  ))
    (joint_type == :Spherical)          && (return (ntrl = 0, nrot = 3, num_of_dof = 3, num_of_cstr = 3, mask_1st = [1,2,3], mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = Int[]  ))
    (joint_type == :FloatingUniversal)  && (return (ntrl = 3, nrot = 2, num_of_dof = 5, num_of_cstr = 1, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1]    )) #t1'*t2
    (joint_type == :OrbitalUniversal)   && (return (ntrl = 2, nrot = 2, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[]  , mask_2nd = [1]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1]    )) #t1'*t2
    (joint_type == :PlanarUniversal)    && (return (ntrl = 2, nrot = 2, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [1]  , mask_3rd_egg = Int[], mask_4th = [1]    )) #t1'*t2
    (joint_type == :PrismaticUniversal) && (return (ntrl = 1, nrot = 2, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [2,3], mask_3rd_egg = Int[], mask_4th = [1]    )) #t1'*t2
    (joint_type == :UniversalPrismatic) && (return (ntrl = 1, nrot = 2, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = [2,3], mask_4th = [1]    )) #t1'*t2
    (joint_type == :Universal)          && (return (ntrl = 0, nrot = 2, num_of_dof = 2, num_of_cstr = 4, mask_1st = [1,2,3], mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1]    )) #t1'*t2
    (joint_type == :FloatingRevolute)   && (return (ntrl = 3, nrot = 1, num_of_dof = 4, num_of_cstr = 2, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [2,3]  )) #t1'*n, t2'*n
    (joint_type == :OrbitalRevolute)    && (return (ntrl = 2, nrot = 1, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = [1]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [2,3]  )) #t1'*n, t2'*n
    (joint_type == :PlanarRevolute)     && (return (ntrl = 2, nrot = 1, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [1]  , mask_3rd_egg = Int[], mask_4th = [2,3]  )) #t1'*n, t2'*n
    (joint_type == :Cylindrical)        && (return (ntrl = 1, nrot = 1, num_of_dof = 2, num_of_cstr = 4, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [2,3], mask_3rd_egg = Int[], mask_4th = [2,3]  )) #t1'*n, t2'*n
    (joint_type == :Revolute)           && (return (ntrl = 0, nrot = 1, num_of_dof = 1, num_of_cstr = 5, mask_1st = [1,2,3], mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [2,3]  )) #t1'*n, t2'*n
    (joint_type == :Floating)           && (return (ntrl = 3, nrot = 0, num_of_dof = 3, num_of_cstr = 3, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Orbital)            && (return (ntrl = 2, nrot = 0, num_of_dof = 2, num_of_cstr = 4, mask_1st = Int[]  , mask_2nd = [1]  , mask_3rd_hen = Int[], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Planar)             && (return (ntrl = 2, nrot = 0, num_of_dof = 2, num_of_cstr = 4, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [1]  , mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Prismatic)          && (return (ntrl = 1, nrot = 0, num_of_dof = 1, num_of_cstr = 5, mask_1st = Int[]  , mask_2nd = Int[], mask_3rd_hen = [2,3], mask_3rd_egg = Int[], mask_4th = [1,2,3])) #t1'*t2, t1'*n, t2'*n
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
刚体铰接约束类。
$(TYPEDEF)
"""
struct PrototypeJoint{hen2eggType,maskType,valueType,cacheType} <: ExtrinsicConstraints{valueType}
    id::Int
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
铰接约束构造子。
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
    nmcs_hen = hen.rbsig.coords.nmcs
    nmcs_egg = egg.rbsig.coords.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    q_hen = cartesian_frame2coords(nmcs_hen,state_hen.origin_position,state_hen.R)
    q_egg = cartesian_frame2coords(nmcs_egg,state_egg.origin_position,state_egg.R)
    if (nmcs_hen isa QCF.QC) && (nmcs_egg isa QCF.QC)
        if mask_3rd == [2,3]
            mask_3rd = [1,2]
        elseif mask_3rd == [1]
            mask_3rd = [3]
        end
        if mask_3rd == [2,3] .+ 3
            mask_3rd = [1,2] .+ 3
        elseif mask_3rd == [1] .+ 3
            mask_3rd = [3] .+ 3
        end
        if mask_4th == [2,3]
            mask_4th = [1,2]
        elseif mask_4th == [1]
            mask_4th = [3]
        end
    end
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
    @show mask_1st,mask_2nd,mask_3rd,mask_4th
    @show values
    PrototypeJoint(
        id,hen2egg,
        num_of_cstr,
        num_of_dof,
        mask_1st,
        mask_2nd,
        mask_3rd,
        mask_4th,
        values,
        cache
    )
end


function make_cstr_function(cst::PrototypeJoint,st::Structure)
    (;indexed,numbered) = st.connectivity
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
    function _inner_cstr_function(q,c)
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
    function inner_cstr_function(q)
        c = get_local_coords(st)
        _inner_cstr_function(q,c)
    end
    inner_cstr_function
end

function make_cstr_jacobian(cst::PrototypeJoint,st::Structure)
    (;indexed,numbered) = st.connectivity
    (;bodyid2sys_free_coords,
      bodyid2sys_full_coords,
      num_of_free_coords,
      num_of_full_coords,
    ) = indexed
    (;sys_loci2coords_idx,bodyid2sys_loci_idx) = numbered
    (;
        num_of_cstr,hen2egg,
        mask_1st,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = cst
    (;hen,egg) = hen2egg
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    free_hen = bodyid2sys_free_coords[id_hen]
    free_egg = bodyid2sys_free_coords[id_egg]
    nmcs_hen = hen.rbsig.coords.nmcs
    nmcs_egg = egg.rbsig.coords.nmcs
    free_idx_hen =  hen.rbsig.coords.free_idx
    free_idx_egg =  egg.rbsig.coords.free_idx
    function _inner_cstr_jacobian(q,c)
        T = eltype(q)
        ret = zeros(T,num_of_cstr,num_of_free_coords)
        q_hen = @view q[bodyid2sys_full_coords[id_hen]]
        q_egg = @view q[bodyid2sys_full_coords[id_egg]]
        # translate
        get_joint_jacobian!(ret,
            nmcs_hen, nmcs_egg,
            hen.rbsig.prop.loci[hen.pid].position,
            egg.rbsig.prop.loci[egg.pid].position,
            cache,
            mask_1st,mask_2nd,mask_3rd,mask_4th,
            free_idx_hen,free_idx_egg,
            free_hen,free_egg,
            q_hen,q_egg
        )
        ret
    end
    function inner_cstr_jacobian(q,c)
        _inner_cstr_jacobian(q,c)
    end
    function inner_cstr_jacobian(q)
        c = get_local_coords(st)
        _inner_cstr_jacobian(q,c)
    end
    inner_cstr_jacobian
end

function get_jointed_free_idx(cst)
    (;num_of_cstr,hen2egg,) = cst
    (;hen,egg) = hen2egg
    free_idx_hen = hen.rbsig.coords.free_idx
    free_idx_egg = egg.rbsig.coords.free_idx
    ncoords_hen = NCF.get_num_of_coords(hen.rbsig.coords.nmcs)
    # ncoords_egg = NCF.get_num_of_coords(egg.rbsig.coords.nmcs)
    free_idx = vcat(
        free_idx_hen,
        free_idx_egg .+ ncoords_hen
    )
end

function get_jointed_free(cst,indexed)
    (;num_of_cstr,hen2egg,) = cst
    (;bodyid2sys_free_coords,bodyid2sys_full_coords,num_of_free_coords) = indexed
    (;hen,egg) = hen2egg
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    mem2sysfree1 = bodyid2sys_free_coords[id_hen]
    mem2sysfree2 = bodyid2sys_free_coords[id_egg]
    cst2sysfree = vcat(
        mem2sysfree1,
        mem2sysfree2
    )
end

function make_cstr_forces_jacobian(cst::PrototypeJoint,st)
    (;
        num_of_cstr,
        cache
    ) = cst
    (;hessians) = cache
    free_idx = get_jointed_free_idx(cst)
    function cstr_forces_jacobian(λ)
        ret = [
            begin
                a = -λ[i] .* hessians[i][free_idx,free_idx]
                # display(a)
                a 
            end
            for i = 1:num_of_cstr
        ]
        sum(ret)
    end
end

function get_jointed_free_idx(cst::LinearJoint)
    free_idx = collect(1:size(cst.A,2))
end

function get_jointed_free(cst::LinearJoint,indexed)
    cst2sysfree = collect(1:indexed.num_of_free_coords)
end

function make_cstr_forces_jacobian(cst::LinearJoint,st)
    free_idx = get_jointed_free_idx(cst)
    function cstr_forces_jacobian(λ)
        zeros(eltype(λ),length(free_idx),length(free_idx))
    end
end
