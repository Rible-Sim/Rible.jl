#todo use the basic types of Constraints
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
    "No. of the axis"
    aid::aidType
end

function ID(rbsig,pid)
    ID(rbsig,pid,1)
end

"""
点对点类。
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct End2End{end1Type<:ID,end2Type<:ID}
    "编号"
    id::Int
    "起始点"
    end1::end1Type
    "终止点"
    end2::end2Type
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
    axes_rot_hen::axesType
    axes_rot_egg::axesType
    mask_1st::maskType
    mask_2nd::maskType
    mask_3rd::maskType
    mask_4th::maskType
    values::valueType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function PrototypeJoint(id,e2e,joint_type::Symbol) 
    hen = e2e.end1
    egg = e2e.end2
    joint_info = get_joint_info(joint_type)
    (;ntrl, nrot, ndof, ncsts, mask_1st, mask_2nd, mask_3rd, mask_4th) = joint_info
    T = get_numbertype(hen.rbsig)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    q_hen = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R)
    q_egg = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R)
    X_hen = NCF.make_X(nmcs_hen,q_hen)
    X_egg = NCF.make_X(nmcs_egg,q_egg)
    # translate on hen
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    axes_trl_hen = Axes(inv(X_hen)*get_orthonormal_axes(X_hen*ā_hen)) 
    trl_hen_n = X_hen*axes_trl_hen.n
    trl_hen_t1 = X_hen*axes_trl_hen.t1
    trl_hen_t2 = X_hen*axes_trl_hen.t2
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    r_hen2egg = C_egg*q_egg .- C_hen*q_hen
    # rotate of egg
    ā_egg = egg.rbsig.prop.ās[egg.aid]
    axes_rot_glb = get_orthonormal_axes(X_egg*ā_egg)
    axes_rot_egg = Axes(inv(X_egg)*axes_rot_glb)
    axes_rot_hen = Axes(inv(X_hen)*axes_rot_glb)
    rot_egg_n = X_egg*axes_rot_egg.n
    rot_egg_t1 = X_egg*axes_rot_egg.t1
    rot_egg_t2 = X_egg*axes_rot_egg.t2
    rot_hen_n = X_hen*axes_rot_hen.n
    rot_hen_t1 = X_hen*axes_rot_hen.t1
    rot_hen_t2 = X_hen*axes_rot_hen.t2
    # first    
    Φ_1st = r_hen2egg
    # second
    Φ_2nd = [
        rot_hen_t1'*rot_egg_t2,
        rot_hen_t1'*rot_egg_n,
        rot_hen_t2'*rot_egg_n
    ]
    # third
    Φ_3rd = [
        trl_hen_n'*r_hen2egg,
        trl_hen_t1'*r_hen2egg,
        trl_hen_t2'*r_hen2egg,
    ]
    # fourth    
    Φ_4th = [r_hen2egg'*r_hen2egg]
    # constraint values
    values = vcat(
        Φ_4th[mask_4th],
        Φ_3rd[mask_3rd], # translate prior to 
        Φ_1st[mask_1st],
        Φ_2nd[mask_2nd], # rotate
    )
    @show joint_info
    @show values
    PrototypeJoint(
        id,e2e,
        ncsts,ndof,
        axes_trl_hen,axes_rot_hen,axes_rot_egg,
        mask_1st, mask_2nd, mask_3rd, mask_4th,
        values
    )
end

function make_Φ(cst::PrototypeJoint,indexed,numbered)
    (;nconstraints,e2e,values,axes_trl_hen,axes_rot_hen,axes_rot_egg,mask_1st, mask_2nd, mask_3rd, mask_4th) = cst
    (;num2sys,mem2num) = numbered
    (;mem2sysfull) = indexed
    hen = e2e.end1
    egg = e2e.end2
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    function inner_Φ(q)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        # translate on hen
        trl_hen_n = X_hen*axes_trl_hen.n
        trl_hen_t1 = X_hen*axes_trl_hen.t1
        trl_hen_t2 = X_hen*axes_trl_hen.t2
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        r_hen2egg = C_egg*q_egg .- C_hen*q_hen
        # rotate of egg
        rot_hen_t1 = X_hen*axes_rot_hen.t1
        rot_hen_t2 = X_hen*axes_rot_hen.t2
        rot_egg_n = X_egg*axes_rot_egg.n
        rot_egg_t2 = X_egg*axes_rot_egg.t2
        # first
        Φ_1st = r_hen2egg
        # second
        Φ_2nd = [
            rot_hen_t1'*rot_egg_t2,
            rot_hen_t1'*rot_egg_n,
            rot_hen_t2'*rot_egg_n
        ]
        # third
        Φ_3rd = [
            trl_hen_n'*r_hen2egg,
            trl_hen_t1'*r_hen2egg,
            trl_hen_t2'*r_hen2egg,
        ]
        # fourth    
        Φ_4th = [r_hen2egg'*r_hen2egg]
        # constraint values
        ret = vcat(
            Φ_4th[mask_4th],
            Φ_3rd[mask_3rd], # translate prior to 
            Φ_1st[mask_1st],
            Φ_2nd[mask_2nd], # rotate
        ) .- values
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        # translate on hen
        trl_hen_n = X_hen*axes_trl_hen.n
        trl_hen_t1 = X_hen*axes_trl_hen.t1
        trl_hen_t2 = X_hen*axes_trl_hen.t2
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)*q_hen
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)*q_egg
        r_hen2egg = C_egg*q_egg .- C_hen*q_hen
        # rotate of egg
        rot_hen_t1 = X_hen*axes_rot_hen.t1
        rot_hen_t2 = X_hen*axes_rot_hen.t2
        rot_egg_n = X_egg*axes_rot_egg.n
        rot_egg_t2 = X_egg*axes_rot_egg.t2
        # first    
        Φ_1st = r_hen2egg
        # second
        Φ_2nd = [
            rot_hen_t1'*rot_egg_t2,
            rot_hen_t1'*rot_egg_n,
            rot_hen_t2'*rot_egg_n
        ]
        # third
        Φ_3rd = [
            trl_hen_n'*r_hen2egg,
            trl_hen_t1'*r_hen2egg,
            trl_hen_t2'*r_hen2egg,
        ]
        # fourth    
        Φ_4th = [r_hen2egg'*r_hen2egg]
        # constraint values
        ret = vcat(
            Φ_4th[mask_4th],
            Φ_3rd[mask_3rd], # translate prior to 
            Φ_1st[mask_1st],
            Φ_2nd[mask_2nd], # rotate
        ) .- values
        ret
    end
    inner_Φ
end

function make_A(cst::PrototypeJoint,indexed,numbered)
    (;nconstraints,e2e,axes_trl_hen,axes_rot_hen,axes_rot_egg, mask_1st, mask_2nd, mask_3rd, mask_4th) = cst
    (;mem2sysfree,mem2sysfull,nfree) = indexed
    (;num2sys,mem2num) = numbered
    hen = e2e.end1
    egg = e2e.end2
    uci_hen =  hen.rbsig.state.cache.free_idx
    uci_egg =  egg.rbsig.state.cache.free_idx
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    free_hen = mem2sysfree[id_hen]
    free_egg = mem2sysfree[id_egg]
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    Y_hen = NCF.get_conversion(nmcs_hen)
    Y_egg = NCF.get_conversion(nmcs_egg)
    function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        # translate on hen
        trl_hen_n = X_hen*axes_trl_hen.n
        trl_hen_t1 = X_hen*axes_trl_hen.t1
        trl_hen_t2 = X_hen*axes_trl_hen.t2
        trl_hen = [trl_hen_n trl_hen_t1 trl_hen_t2]
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        r_hen2egg = C_egg*q_egg .- C_hen*q_hen
        # rotate of egg
        rot_hen_t1 = X_hen*axes_rot_hen.t1
        rot_hen_t2 = X_hen*axes_rot_hen.t2
        rot_egg_n = X_egg*axes_rot_egg.n
        rot_egg_t2 = X_egg*axes_rot_egg.t2
        # jac
        o3 = zeros(eltype(q),3)
        # jac first
        A_1st = zeros(eltype(q),3,nfree)
        A_1st[:,free_hen] .= -C_hen[:,uci_hen]
        A_1st[:,free_egg] .=  C_egg[:,uci_egg]
        # jac third
        A_3rd = zeros(eltype(q),3,nfree)
        A_3rd[:,free_hen] .= kron(
                hcat(
                    o3, transpose(trl_hen)
                ),
                transpose(r_hen2egg)
            )[:,uci_hen]
        A_3rd[:,free_hen] .-= transpose(trl_hen)*C_hen[:,uci_hen]
        A_3rd[:,free_egg] .+= transpose(trl_hen)*C_egg[:,uci_egg]
        # jac 2nd
        A_2nd = zeros(eltype(q),3,nfree)
        A_2nd[1,free_hen] = (
                transpose(
                    kron(vcat(0,axes_rot_hen.t1),rot_egg_t2)
                )*Y_hen
            )[:,uci_hen]
        A_2nd[1,free_egg] = (
                transpose(
                    kron(vcat(0,axes_rot_egg.t2),rot_hen_t1)
                )*Y_egg
            )[:,uci_egg]
        A_2nd[2,free_hen] = (
                transpose(
                    kron(vcat(0,axes_rot_hen.t1),rot_egg_n)
                )*Y_hen
            )[:,uci_hen]
        A_2nd[2,free_egg] = (
                transpose(
                    kron(vcat(0,axes_rot_egg.n),rot_hen_t1)
                )*Y_egg
            )[:,uci_egg]
        A_2nd[3,free_hen] = (
                transpose(
                    kron(vcat(0,axes_rot_hen.t2),rot_egg_n)
                )*Y_hen
            )[:,uci_hen]
        A_2nd[3,free_egg] = (
                transpose(
                    kron(vcat(0,axes_rot_egg.n),rot_hen_t2)
                )*Y_egg
            )[:,uci_egg]
        # jac 4th
        # to be implemented
        A_4th = zeros(eltype(q),1,nfree)
        A_4th[1,free_hen] = -2r_hen2egg'*C_hen[:,uci_hen]
        A_4th[1,free_egg] =  2r_hen2egg'*C_egg[:,uci_egg]
        # constraint values
        ret = vcat(
            A_4th[mask_4th,:],
            A_3rd[mask_3rd,:], # translate prior to 
            A_1st[mask_1st,:],
            A_2nd[mask_2nd,:], # rotate
        )
        ret
    end
    function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        # translate on hen
        trl_hen_n = X_hen*axes_trl_hen.n
        trl_hen_t1 = X_hen*axes_trl_hen.t1
        trl_hen_t2 = X_hen*axes_trl_hen.t2
        trl_hen = [trl_hen_n trl_hen_t1 trl_hen_t2]
        c_hen = c[num2sys[mem2num[id_hen][hen.pid]]]
        c_egg = c[num2sys[mem2num[id_egg][egg.pid]]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)
        r_hen2egg = C_egg*q_egg .- C_hen*q_hen
        # rotate of egg
        rot_hen_t1 = X_hen*axes_rot_hen.t1
        rot_hen_t2 = X_hen*axes_rot_hen.t2
        rot_egg_n = X_egg*axes_rot_egg.n
        rot_egg_t2 = X_egg*axes_rot_egg.t2
        # jac
        o3 = zeros(eltype(q),3)
        # jac first
        A_1st = zeros(eltype(q),3,nfree)
        A_1st[:,free_hen] .= -C_hen[:,uci_hen]
        A_1st[:,free_egg] .=  C_egg[:,uci_egg]
        # jac third
        A_3rd = zeros(eltype(q),3,nfree)
        A_3rd[:,free_hen] .= kron(
                hcat(
                    o3, transpose(trl_hen)
                ),
                transpose(r_hen2egg)
            )[:,uci_hen]
        A_3rd[:,free_hen] .-= transpose(trl_hen)*C_hen[:,uci_hen]
        A_3rd[:,free_egg] .+= transpose(trl_hen)*C_egg[:,uci_egg]
        # jac 2nd
        A_2nd = zeros(eltype(q),3,nfree)
        A_2nd[1,free_hen] = (
                transpose(
                    kron(vcat(0,axes_rot_hen.t1),rot_egg_t2)
                )*Y_hen
            )[:,uci_hen]
        A_2nd[1,free_egg] = (
                transpose(
                    kron(vcat(0,axes_rot_egg.t2),rot_hen_t1)
                )*Y_egg
            )[:,uci_egg]
        A_2nd[2,free_hen] = (
                transpose(
                    kron(vcat(0,axes_rot_hen.t1),rot_egg_n)
                )*Y_hen
            )[:,uci_hen]
        A_2nd[2,free_egg] = (
                transpose(
                    kron(vcat(0,axes_rot_egg.n),rot_hen_t1)
                )*Y_egg
            )[:,uci_egg]
        A_2nd[3,free_hen] = (
                transpose(
                    kron(vcat(0,axes_rot_hen.t2),rot_egg_n)
                )*Y_hen
            )[:,uci_hen]
        A_2nd[3,free_egg] = (
                transpose(
                    kron(vcat(0,axes_rot_egg.n),rot_hen_t2)
                )*Y_egg
            )[:,uci_egg]
        # jac 4th
        A_4th = zeros(eltype(q),1,nfree)
        A_4th[1,free_hen] = -2r_hen2egg'*C_hen[:,uci_hen]
        A_4th[1,free_egg] =  2r_hen2egg'*C_egg[:,uci_egg]
        # constraint values
        ret = vcat(
            A_4th[mask_4th,:],
            A_3rd[mask_3rd,:], # translate prior to 
            A_1st[mask_1st,:],
            A_2nd[mask_2nd,:], # rotate
        )
        ret
    end
    inner_A
end

function make_Φqᵀq(cst::PrototypeJoint,numbered,c)
    (;e2e,axes_trl_hen,axes_rot_hen,axes_rot_egg,mask_1st,mask_2nd,mask_3rd,mask_4th,) = cst
    (;mem2num,num2sys) = numbered
    hen = e2e.end1
    egg = e2e.end2
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    nld_hen = NCF.get_nlocaldim(hen.rbsig.state.cache.funcs.nmcs)
    nld_egg = NCF.get_nlocaldim(egg.rbsig.state.cache.funcs.nmcs)
    Y_hen = NCF.get_conversion(hen.rbsig.state.cache.funcs.nmcs)
    Y_egg = NCF.get_conversion(egg.rbsig.state.cache.funcs.nmcs)
    cv = BlockDiagonal([Y_hen,Y_egg])
    ndim = NCF.get_ndim(hen.rbsig.state.cache.funcs.nmcs)
    T = get_numbertype(cst)
    I_Int = NCF.make_I(Int,ndim)
    # translate on hen
    c_hen = c[num2sys[mem2num[id_hen][hen.pid]]]
    c_egg = c[num2sys[mem2num[id_egg][egg.pid]]]
    C_hen = hen.rbsig.state.cache.funcs.C(c_hen)
    C_egg = egg.rbsig.state.cache.funcs.C(c_egg)
    # first
    Φqᵀq_1st =[
        begin
            ret = kron(
                zeros(
                    T,
                    (1+nld_hen)+(1+nld_egg),
                    (1+nld_hen)+(1+nld_egg),
                ),
                I_Int
            )
        end
        for i = 1:3
    ]
    Φqᵀq_3rd = [
        begin
        d̄ = axes_trl_hen.X[:,i]
            ret_raw = zeros(
                T,
                (1+nld_hen)+(1+nld_egg),
                (1+nld_hen)+(1+nld_egg),
            )
            ret_raw[1,2:1+nld_hen] = -d̄
            ret_raw[2:1+nld_hen,1] = -d̄
            ret_raw[2:1+nld_hen,2:1+nld_hen] = -c_hen*d̄'-d̄*c_hen'
            ret_raw[(1+nld_hen)+1,2:1+nld_egg] = d̄
            ret_raw[2:1+nld_egg,(1+nld_hen)+1] = d̄
            ret_raw[(1+nld_hen)+2:end,2:1+nld_egg] = c_egg*d̄'
            ret_raw[2:1+nld_egg,(1+nld_hen)+2:end] = d̄*c_egg'
            transpose(cv)*kron(ret_raw,I_Int)*cv
        end
        for i = 1:3
    ]
    Φqᵀq_2nd = [
        begin
            ret_raw = zeros(
                T,
                (1+nld_hen)+(1+nld_egg),
                (1+nld_hen)+(1+nld_egg),
            )
            ret_raw[2:1+nld_hen,(1+nld_hen)+2:(1+nld_hen)+1+nld_egg] = axes_rot_hen.t1*axes_rot_egg.t2'
            ret_raw[(1+nld_hen)+2:(1+nld_hen)+1+nld_egg,2:1+nld_hen] = axes_rot_egg.t2*axes_rot_hen.t1'
            transpose(cv)*kron(ret_raw,I_Int)*cv
        end,
        begin
            ret_raw = zeros(
                T,
                (1+nld_hen)+(1+nld_egg),
                (1+nld_hen)+(1+nld_egg),
            )
            ret_raw[2:1+nld_hen,(1+nld_hen)+2:(1+nld_hen)+1+nld_egg] = axes_rot_hen.t1*axes_rot_egg.n'
            ret_raw[(1+nld_hen)+2:(1+nld_hen)+1+nld_egg,2:1+nld_hen] = axes_rot_egg.n*axes_rot_hen.t1'
            transpose(cv)*kron(ret_raw,I_Int)*cv
        end,
        begin
            ret_raw = zeros(
                T,
                (1+nld_hen)+(1+nld_egg),
                (1+nld_hen)+(1+nld_egg),
            )
            ret_raw[2:1+nld_hen,(1+nld_hen)+2:(1+nld_hen)+1+nld_egg] = axes_rot_hen.t2*axes_rot_egg.n'
            ret_raw[(1+nld_hen)+2:(1+nld_hen)+1+nld_egg,2:1+nld_hen] = axes_rot_egg.n*axes_rot_hen.t2'
            transpose(cv)*kron(ret_raw,I_Int)*cv
        end,
    ]    
    Φqᵀq_4th = [
        begin
            ret = kron(
                zeros(
                    T,
                    (1+nld_hen)+(1+nld_egg),
                    (1+nld_hen)+(1+nld_egg),
                ),
                I_Int
            )
            ret[1:(1+nld_hen)*ndim,1:(1+nld_hen)*ndim] = 2C_hen'C_hen
            ret[(1+nld_hen)*ndim+1:(1+nld_hen)*ndim+(1+nld_egg)*ndim,1:(1+nld_hen)*ndim] = -2C_egg'C_hen
            ret[1:(1+nld_hen)*ndim,(1+nld_hen)*ndim+1:(1+nld_hen)*ndim+(1+nld_egg)*ndim] = -2C_hen'C_egg
            ret[(1+nld_hen)*ndim+1:(1+nld_hen)*ndim+(1+nld_egg)*ndim,(1+nld_hen)*ndim+1:(1+nld_hen)*ndim+(1+nld_egg)*ndim] = 2C_egg'C_egg
            ret
        end
    ]
    Φqᵀq = vcat(
        Φqᵀq_4th[mask_4th],
        Φqᵀq_3rd[mask_3rd],
        Φqᵀq_1st[mask_1st],
        Φqᵀq_2nd[mask_2nd],
    )
    Φqᵀq
end

function get_jointed_free_idx(cst)
    (;nconstraints,e2e,) = cst
    hen = e2e.end1
    egg = e2e.end2
    uci_hen = hen.rbsig.state.cache.free_idx
    uci_egg = egg.rbsig.state.cache.free_idx
    ncoords1 = NCF.get_ncoords(hen.rbsig.state.cache.funcs.nmcs)
    # ncoords2 = NCF.get_ncoords(egg.rbsig.state.cache.nmcs)
    uci = vcat(
        uci_hen,
        uci_egg .+ ncoords1
    )
end

function get_jointed_free(cst,indexed)
    (;nconstraints,e2e,) = cst
    (;mem2sysfree,mem2sysfull,nfree) = indexed
    hen = e2e.end1
    egg = e2e.end2
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    mem2sysfree1 = mem2sysfree[id_hen]
    mem2sysfree2 = mem2sysfree[id_egg]
    cst2sysfree = vcat(
        mem2sysfree1,
        mem2sysfree2
    )
end

function make_∂Aᵀλ∂q(cst,numbered,c)
    (;nconstraints) = cst
    uci = get_jointed_free_idx(cst)
    Φqᵀq = make_Φqᵀq(cst,numbered,c)
    function ∂Aᵀλ∂q(λ)
        ret = [
            begin
                a = Φqᵀq[i][uci,uci] .* λ[i]
                # display(a)
                a 
            end
            for i = 1:nconstraints
        ]
        sum(ret)
    end
end

function get_jointed_free_idx(cst::LinearJoint)
    uci = collect(1:size(cst.A,2))
end

function get_jointed_free(cst::LinearJoint,indexed)
    cst2sysfree = collect(1:indexed.nfree)
end

function make_∂Aᵀλ∂q(cst::LinearJoint,numbered,c)
    uci = get_jointed_free_idx(cst)
    function ∂Aᵀλ∂q(λ)
        zeros(eltype(λ),length(uci),length(uci))
    end
end

function get_joint_info(joint_type::Symbol)
    (joint_type == :FloatingSpherical)  && (return (ntrl = 3, nrot = 3, ndof = 6, ncsts = 0, mask_1st = Int[],   mask_2nd = Int[]     , mask_3rd = Int[], mask_4th = Int[] ))
    (joint_type == :OrbitalSpherical)   && (return (ntrl = 2, nrot = 3, ndof = 5, ncsts = 1, mask_1st = Int[],   mask_2nd = Int[]     , mask_3rd = Int[], mask_4th = [1]   ))
    (joint_type == :PlanarSpherical)    && (return (ntrl = 2, nrot = 3, ndof = 5, ncsts = 1, mask_1st = Int[],   mask_2nd = Int[]     , mask_3rd = [1],   mask_4th = Int[] ))
    (joint_type == :PrismaticSpherical) && (return (ntrl = 1, nrot = 3, ndof = 4, ncsts = 2, mask_1st = Int[],   mask_2nd = Int[]     , mask_3rd = [2,3], mask_4th = Int[] ))
    (joint_type == :Spherical)          && (return (ntrl = 0, nrot = 3, ndof = 3, ncsts = 3, mask_1st = [1,2,3], mask_2nd = Int[]     , mask_3rd = Int[], mask_4th = Int[] ))
    (joint_type == :FloatingUniversal)  && (return (ntrl = 3, nrot = 2, ndof = 5, ncsts = 1, mask_1st = Int[],   mask_2nd = [1]       , mask_3rd = Int[], mask_4th = Int[] )) #t1'*t2
    (joint_type == :OrbitalUniversal)   && (return (ntrl = 2, nrot = 2, ndof = 4, ncsts = 2, mask_1st = Int[],   mask_2nd = [1]       , mask_3rd = Int[], mask_4th = [1]   )) #t1'*t2
    (joint_type == :PlanarUniversal)    && (return (ntrl = 2, nrot = 2, ndof = 4, ncsts = 2, mask_1st = Int[],   mask_2nd = [1]       , mask_3rd = [1],   mask_4th = Int[] )) #t1'*t2
    (joint_type == :PrismaticUniversal) && (return (ntrl = 1, nrot = 2, ndof = 3, ncsts = 3, mask_1st = Int[],   mask_2nd = [1]       , mask_3rd = [2,3], mask_4th = Int[] )) #t1'*t2
    (joint_type == :Universal)          && (return (ntrl = 0, nrot = 2, ndof = 2, ncsts = 4, mask_1st = [1,2,3], mask_2nd = [1]       , mask_3rd = Int[], mask_4th = Int[] )) #t1'*t2
    (joint_type == :FloatingRevolute)   && (return (ntrl = 3, nrot = 1, ndof = 4, ncsts = 2, mask_1st = Int[],   mask_2nd = [2,3]     , mask_3rd = Int[], mask_4th = Int[] )) #t1'*n, t2'*n
    (joint_type == :OrbitalRevolute)    && (return (ntrl = 2, nrot = 1, ndof = 3, ncsts = 3, mask_1st = Int[],   mask_2nd = [2,3]     , mask_3rd = Int[], mask_4th = [1]   )) #t1'*n, t2'*n
    (joint_type == :PlanarRevolute)     && (return (ntrl = 2, nrot = 1, ndof = 3, ncsts = 3, mask_1st = Int[],   mask_2nd = [2,3]     , mask_3rd = [1],   mask_4th = Int[] )) #t1'*n, t2'*n
    (joint_type == :Cylindrical)        && (return (ntrl = 1, nrot = 1, ndof = 2, ncsts = 4, mask_1st = Int[],   mask_2nd = [2,3]     , mask_3rd = [2,3], mask_4th = Int[] )) #t1'*n, t2'*n
    (joint_type == :Revolute)           && (return (ntrl = 0, nrot = 1, ndof = 1, ncsts = 5, mask_1st = [1,2,3], mask_2nd = [2,3]     , mask_3rd = Int[], mask_4th = Int[] )) #t1'*n, t2'*n
    (joint_type == :Floating)           && (return (ntrl = 3, nrot = 0, ndof = 3, ncsts = 3, mask_1st = Int[],   mask_2nd = [1,2,3]   , mask_3rd = Int[], mask_4th = Int[] )) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Orbital)            && (return (ntrl = 2, nrot = 0, ndof = 2, ncsts = 4, mask_1st = Int[],   mask_2nd = [1,2,3]   , mask_3rd = Int[], mask_4th = [1]   )) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Planar)             && (return (ntrl = 2, nrot = 0, ndof = 2, ncsts = 4, mask_1st = Int[],   mask_2nd = [1,2,3]   , mask_3rd = [1],   mask_4th = Int[] )) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Prismatic)          && (return (ntrl = 1, nrot = 0, ndof = 1, ncsts = 5, mask_1st = Int[],   mask_2nd = [1,2,3]   , mask_3rd = [2,3], mask_4th = Int[] )) #t1'*t2, t1'*n, t2'*n
    (joint_type == :Fixed)              && (return (ntrl = 0, nrot = 0, ndof = 0, ncsts = 6, mask_1st = [1,2,3], mask_2nd = [1,2,3]   , mask_3rd = Int[], mask_4th = Int[] )) #t1'*t2, t1'*n, t2'*n
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