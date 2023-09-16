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

function find_order(nmcs_hen,nmcs_egg,q_hen,q_egg,X_hen,X_egg,nΦ1,nΦ2)
    T = eltype(q_hen)
    nq1 = length(q_hen)
    nq2 = length(q_egg)
    A = zeros(T,nΦ1+nΦ2+9,nq1+nq2)
    A[    1:nΦ1,          1:nq1    ] = NCF.make_Φq(nmcs_hen,collect(1:nq1),collect(1:nΦ1))(q_hen)
    A[nΦ1+1:nΦ1+nΦ2,  nq1+1:nq1+nq2] = NCF.make_Φq(nmcs_egg,collect(1:nq2),collect(1:nΦ2))(q_egg)
    k =     nΦ1+nΦ2    
    A_rest = @view A[k+1:end,:]
    rot_jac!(A_rest,collect(1:9),
        q_hen,q_egg,
        X_hen,X_egg,
        collect(1:nq1),
        collect(nq1+1:nq1+nq2),
        collect(1:nq1),
        collect(1:nq2),
    )
    _,pidx = rref_with_pivots!(Matrix(transpose(A)),min(size(A)...)*eps(real(float(one(eltype(A))))))
    @assert length(pidx) >= 15
    order = pidx[k+1:k+3] .- 12
    order
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
刚体坐标固定约束类，适用于单个坐标。
$(TYPEDEF)
"""
struct FixedIndicesConstraint{T} <: ExternalConstraints{T}
    nconstraints::Int64
    indices::Vector{Int64}
    values::T
end

"""
刚体坐标固定约束构造子。
$(TYPEDSIGNATURES)
"""
function FixedIndicesConstraint(indices,values)
    FixedIndicesConstraint(length(indices),indices,values)
end

function make_Φ(cst::FixedIndicesConstraint,indexed,numbered)
    @unpack indices, values = cst
    @inline @inbounds inner_Φ(q)   = q[indices]-values
    @inline @inbounds inner_Φ(q,d) = q[indices]-d
    inner_Φ
end

function make_A(cst::FixedIndicesConstraint,indexed,numbered)
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
struct PinJoint{valueType,e2eType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::valueType
    e2e::e2eType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function PinJoint(id,e2e)
    rb1 = e2e.hen.rbsig
    nΦ = get_ndim(rb1)
    T = get_numbertype(rb1)
    values = zeros(T,nΦ)
    PinJoint(id,nΦ,values,e2e)
end

function make_Φ(cst::PinJoint,indexed,numbered)
    (;nconstraints,values,e2e) = cst
    (;mem2num,num2sys) = numbered
    (;mem2sysfull) = indexed
    hen = e2e.end1
    egg = e2e.end2
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    function inner_Φ(q)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[hen.rbsig.prop.id]]
        q_egg = @view q[mem2sysfull[egg.rbsig.prop.id]]
        ret .= C_hen*q_hen.-C_egg*q_egg
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        id_hen = hen.rbsig.prop.id
        id_egg = egg.rbsig.prop.id
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        c_hen = c[num2sys[mem2num[id_hen][hen.pid]]]
        c_egg = c[num2sys[mem2num[id_egg][egg.pid]]]
        ret .= hen.rbsig.state.cache.funcs.C(c_hen)*q_hen .-
               egg.rbsig.state.cache.funcs.C(c_egg)*q_egg
        ret
    end
    inner_Φ
end

function make_A(cst::PinJoint,indexed,numbered)
    (;nconstraints,values,e2e) = cst
    (;mem2sysfree,nfree) = indexed
    (;mem2num,num2sys) = numbered
    hen = e2e.end1
    egg = e2e.end2
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    uci_hen =  hen.rbsig.state.cache.free_idx
    uci_egg =  egg.rbsig.state.cache.free_idx
    function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nfree)
        ret[:,mem2sysfree[hen.rbsig.prop.id]] =  C_hen[:,uci_hen]
        ret[:,mem2sysfree[egg.rbsig.prop.id]] = -C_egg[:,uci_egg]
        ret
    end
    function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
        id_hen = hen.rbsig.prop.id
        id_egg = egg.rbsig.prop.id
        c_hen = c[num2sys[mem2num[id_hen][hen.pid]]]
        c_egg = c[num2sys[mem2num[id_egg][egg.pid]]]
        ret[:,mem2sysfree[hen.rbsig.prop.id]] =  hen.rbsig.state.cache.funcs.C(c_hen)[:,uci_hen]
        ret[:,mem2sysfree[egg.rbsig.prop.id]] = -egg.rbsig.state.cache.funcs.C(c_egg)[:,uci_egg]
        ret
    end
    inner_A
end

"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct RevoluteJoint{valueType,e2eType,maskType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::valueType
    e2e::e2eType
    mask::maskType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function RevoluteJoint(id,e2e)
    hen = e2e.end1
    egg = e2e.end2
    rb1 = hen.rbsig
    nΦ = 5
    T = get_numbertype(rb1)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    q_hen,_ = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R,state_hen.ṙo,state_hen.ω)
    q_egg,_ = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R,state_egg.ṙo,state_egg.ω)
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    values = zeros(T,nΦ)
    values[1:3] .= C_hen*q_hen .- C_egg*q_egg
    X_hen = NCF.make_X(nmcs_hen,q_hen)
    X_egg = NCF.make_X(nmcs_egg,q_egg)
    a_hen = X_hen*ā_hen
    inprods = transpose(X_egg)*a_hen
    mask = 1:3 .!== argmax(abs.(inprods))
    values[4:5] = inprods[mask]
    RevoluteJoint(id,nΦ,values,e2e,mask)
end

function make_Φ(cst::RevoluteJoint,indexed,numbered)
    (;nconstraints,values,e2e,mask) = cst
    (;mem2num,num2sys) = numbered
    (;mem2sysfull) = indexed
    hen = e2e.end1
    egg = e2e.end2
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs


    
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    function inner_Φ(q)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[hen.rbsig.prop.id]]
        q_egg = @view q[mem2sysfull[egg.rbsig.prop.id]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        ret[1:3] .= C_hen*q_hen.-C_egg*q_egg - values[1:3]
        ret[4:5] = transpose(X_egg[:,mask])*a_hen - values[4:5]
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        id_hen = hen.rbsig.prop.id
        id_egg = egg.rbsig.prop.id
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        c_hen = c[num2sys[mem2num[id_hen][hen.pid]]]
        c_egg = c[num2sys[mem2num[id_egg][egg.pid]]]
        ret[1:3] .= hen.rbsig.state.cache.funcs.C(c_hen)*q_hen .-
                    egg.rbsig.state.cache.funcs.C(c_egg)*q_egg
        ret[1:3] .= C_hen*q_hen.-C_egg*q_egg .- values[1:3]
        ret[4:5] = transpose(X_egg[:,mask])*a_hen - values[4:5]
        ret
    end
    inner_Φ
end

function make_A(cst::RevoluteJoint,indexed,numbered)
    (;nconstraints,e2e,mask) = cst
    (;mem2sysfree,mem2sysfull,nfree) = indexed
    (;mem2num,num2sys) = numbered
    hen = e2e.end1
    egg = e2e.end2
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    free_hen = mem2sysfree[id_hen]
    free_egg = mem2sysfree[id_egg]
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs

    

    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    uci_hen =  hen.rbsig.state.cache.free_idx
    uci_egg =  egg.rbsig.state.cache.free_idx
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    Y_hen = NCF.get_conversion(hen.rbsig.state.cache.funcs.nmcs)
    Y_egg = NCF.get_conversion(egg.rbsig.state.cache.funcs.nmcs)
    function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        ret[1:3,free_hen] .=  C_hen[:,uci_hen]
        ret[1:3,free_egg] .= -C_egg[:,uci_egg]
        o3 = zero(a_hen)
        ret[4:5,free_hen] .= (transpose(
            kron(vcat(0,ā_hen),X_egg)[:,mask])*Y_hen)[:,uci_hen]
        ret[4:5,free_egg] .= (transpose(
            [
                o3 o3 o3;
                a_hen o3 o3; 
                o3 a_hen o3;  
                o3 o3 a_hen; 
            ][:,mask])*Y_egg)[:,uci_egg]
        ret
    end
    function inner_A(q,c)
        T = eltype(q)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        c_hen = c[num2sys[mem2num[id_hen][hen.pid]]]
        c_egg = c[num2sys[mem2num[id_egg][egg.pid]]]
        ret[1:3,free_hen] .=  hen.rbsig.state.cache.funcs.C(c_hen)[:,uci_hen]
        ret[1:3,free_egg] .= -egg.rbsig.state.cache.funcs.C(c_egg)[:,uci_egg]
        o3 = zero(a_hen)
        ret[4:5,free_hen] .= (transpose(
            kron(vcat(0,ā_hen),X_egg)[:,mask])*Y_hen)[:,uci_hen]
        ret[4:5,free_egg] .= (transpose(
            [
                o3 o3 o3;
                a_hen o3 o3;
                o3 a_hen o3;  
                o3 o3 a_hen;
            ][:,mask])*Y_egg)[:,uci_egg]
        ret
    end
    inner_A
end


"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct FloatingRevoluteJoint{valueType,e2eType,maskType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::valueType
    e2e::e2eType
    mask::maskType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function FloatingRevoluteJoint(id,e2e)
    hen = e2e.end1
    egg = e2e.end2
    rb1 = hen.rbsig
    nΦ = 2
    T = get_numbertype(rb1)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    q_hen,_ = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R,state_hen.ṙo,state_hen.ω)
    q_egg,_ = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R,state_egg.ṙo,state_egg.ω)
    values = zeros(T,nΦ)
    X_hen = NCF.make_X(nmcs_hen,q_hen)
    X_egg = NCF.make_X(nmcs_egg,q_egg)
    a_hen = X_hen*ā_hen
    inprods = transpose(X_egg)*a_hen
    mask = 1:3 .!== argmax(abs.(inprods))
    values[1:2] = inprods[mask]
    FloatingRevoluteJoint(id,nΦ,values,e2e,mask)
end

function make_Φ(cst::FloatingRevoluteJoint,indexed,numbered)
    (;nconstraints,values,e2e,mask) = cst
    (;mem2num,num2sys) = numbered
    (;mem2sysfull) = indexed
    hen = e2e.end1
    egg = e2e.end2
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    function inner_Φ(q)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[hen.rbsig.prop.id]]
        q_egg = @view q[mem2sysfull[egg.rbsig.prop.id]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        ret[1:2] = transpose(X_egg[:,mask])*a_hen - values[1:2]
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        id_hen = hen.rbsig.prop.id
        id_egg = egg.rbsig.prop.id
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        ret[1:2] = transpose(X_egg[:,mask])*a_hen - values[1:2]
        ret
    end
    inner_Φ
end

function make_A(cst::FloatingRevoluteJoint,indexed,numbered)
    (;nconstraints,e2e,mask) = cst
    (;mem2sysfree,mem2sysfull,nfree) = indexed
    (;mem2num,num2sys) = numbered
    hen = e2e.end1
    egg = e2e.end2
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    free_hen = mem2sysfree[id_hen]
    free_egg = mem2sysfree[id_egg]
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    
    uci_hen =  hen.rbsig.state.cache.free_idx
    uci_egg =  egg.rbsig.state.cache.free_idx
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    Y_hen = NCF.get_conversion(hen.rbsig.state.cache.funcs.nmcs)
    Y_egg = NCF.get_conversion(egg.rbsig.state.cache.funcs.nmcs)
    function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        o3 = zero(a_hen)
        ret[1:2,free_hen] .= (transpose(
            kron(vcat(0,ā_hen),X_egg)[:,mask])*Y_hen)[:,uci_hen]
        ret[1:2,free_egg] .= (transpose(
            [
                o3 o3 o3;
                a_hen o3 o3; 
                o3 a_hen o3;  
                o3 o3 a_hen; 
            ][:,mask])*Y_egg)[:,uci_egg]
        ret
    end
    function inner_A(q,c)
        T = eltype(q)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        o3 = zero(a_hen)
        ret[1:2,free_hen] .= (transpose(
            kron(vcat(0,ā_hen),X_egg)[:,mask])*Y_hen)[:,uci_hen]
        ret[1:2,free_egg] .= (transpose(
            [
                o3 o3 o3;
                a_hen o3 o3;
                o3 a_hen o3;  
                o3 o3 a_hen;
            ][:,mask])*Y_egg)[:,uci_egg]
        ret
    end
    inner_A
end


"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct PrismaticSphericalJoint{valueType,e2eType,axesType,orderType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::valueType
    e2e::e2eType
    axes::axesType
    order::orderType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function PrismaticSphericalJoint(id,e2e)
    hen = e2e.end1
    egg = e2e.end2
    nΦ = 2
    T = get_numbertype(hen.rbsig)


    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    axes = SpatialFrame(ā_hen)
    q_hen = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R)
    q_egg = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R)
    X_hen = NCF.make_X(nmcs_hen,q_hen)
    a1t1 = X_hen*axes.t1
    a1t2 = X_hen*axes.t2
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    values = zeros(T,nΦ)
    values[1:2] .= transpose([a1t1 a1t2])*(C_hen*q_hen.-C_egg*q_egg)
    PrismaticSphericalJoint(id,nΦ,values,e2e,axes,order)
end

function make_Φ(cst::PrismaticSphericalJoint,indexed,numbered)
    (;nconstraints,values,e2e,axes,order) = cst
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
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        a1t1 = X_hen*axes.t1
        a1t2 = X_hen*axes.t2
        ret[1:2] .= transpose([a1t1 a1t2])*(C_hen*q_hen.-C_egg*q_egg) .- values[1:2]
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a1t1 = X_hen*axes.t1
        a1t2 = X_hen*axes.t2
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)*q_hen
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)*q_egg
        ret[1:2] .= transpose([a1t1 a1t2])*(C_hen*q_hen.-C_egg*q_egg) .- values[1:2]
        ret
    end
    inner_Φ
end

function make_A(cst::PrismaticSphericalJoint,indexed,numbered)
    (;nconstraints,e2e,axes,order) = cst
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
    function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a1t1 = X_hen*axes.t1
        a1t2 = X_hen*axes.t2
        ret[1:2,free_hen] .= kron(
                [
                    0 transpose(axes.t1); 
                    0 transpose(axes.t2)
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1:2,free_hen] .+= transpose([a1t1 a1t2])*C_hen[:,uci_hen]
        ret[1:2,free_egg] .-= transpose([a1t1 a1t2])*C_egg[:,uci_egg]
        ret
    end
    function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a1t1 = X_hen*axes.t1
        a1t2 = X_hen*axes.t2
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)
        ret[1:2,free_hen] .= kron(
                [
                    0 transpose(axes.t1); 
                    0 transpose(axes.t2)
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1:2,free_hen] .+= transpose([a1t1 a1t2])*C_hen[:,uci_hen]
        ret[1:2,free_egg] .-= transpose([a1t1 a1t2])*C_egg[:,uci_egg]
        ret
    end
    inner_A
end


"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct PrismaticJoint{valueType,e2eType,axesType,orderType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::valueType
    e2e::e2eType
    axes::axesType
    order::orderType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function PrismaticJoint(id,e2e)
    hen = e2e.end1
    egg = e2e.end2
    nΦ = 5
    T = get_numbertype(hen.rbsig)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    axes = SpatialFrame(ā_hen)
    q_hen = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R)
    q_egg = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R)
    X_hen = NCF.make_X(nmcs_hen,q_hen)
    X_egg = NCF.make_X(nmcs_egg,q_egg)
    a1t1 = X_hen*axes.t1
    a1t2 = X_hen*axes.t2
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    values = zeros(T,nΦ)
    values[1:2] .= transpose([a1t1 a1t2])*(C_hen*q_hen.-C_egg*q_egg)
    nΦ1 = NCF.get_nconstraints(nmcs_hen)
    nΦ2 = NCF.get_nconstraints(nmcs_egg)
    inprods = transpose(X_egg)*X_hen
    order = find_order(nmcs_hen,nmcs_egg,q_hen,q_egg,X_hen,X_egg,nΦ1,nΦ2)
    values[3:5] .= inprods[order]
    PrismaticJoint(id,nΦ,values,e2e,axes,order)
end

function make_Φ(cst::PrismaticJoint,indexed,numbered)
    (;nconstraints,values,e2e,axes,order) = cst
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
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        a1t1 = X_hen*axes.t1
        a1t2 = X_hen*axes.t2
        ret[1:2] .= transpose([a1t1 a1t2])*(C_hen*q_hen.-C_egg*q_egg) .- values[1:2]
        inprods = transpose(X_egg)*X_hen
        ret[3:5] .= inprods[order] .- values[3:5]
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a1t1 = X_hen*axes.t1
        a1t2 = X_hen*axes.t2
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)*q_hen
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)*q_egg
        ret[1:2] .= transpose([a1t1 a1t2])*(C_hen*q_hen.-C_egg*q_egg) .- values[1:2]
        inprods = transpose(X_egg)*X_hen
        ret[3:5] .= inprods[order] .- values[3:5]
        ret
    end
    inner_Φ
end

function make_A(cst::PrismaticJoint,indexed,numbered)
    (;nconstraints,e2e,axes,order) = cst
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
    function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a1t1 = X_hen*axes.t1
        a1t2 = X_hen*axes.t2
        ret[1:2,free_hen] .= kron(
                [
                    0 transpose(axes.t1); 
                    0 transpose(axes.t2)
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1:2,free_hen] .+= transpose([a1t1 a1t2])*C_hen[:,uci_hen]
        ret[1:2,free_egg] .-= transpose([a1t1 a1t2])*C_egg[:,uci_egg]
        k = 2
        ret_rest = @view ret[k+1:end,:]
        rot_jac!(ret_rest,order,q_hen,q_egg,X_hen,X_egg,free_hen,free_egg,uci_hen,uci_egg,)
        ret
    end
    function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a1t1 = X_hen*axes.t1
        a1t2 = X_hen*axes.t2
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)
        ret[1:2,free_hen] .= kron(
                [
                    0 transpose(axes.t1); 
                    0 transpose(axes.t2)
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1:2,free_hen] .+= transpose([a1t1 a1t2])*C_hen[:,uci_hen]
        ret[1:2,free_egg] .-= transpose([a1t1 a1t2])*C_egg[:,uci_egg]
        k = 2
        ret_rest = @view ret[k+1:end,:]
        rot_jac!(ret_rest,order,q_hen,q_egg,X_hen,X_egg,free_hen,free_egg,uci_hen,uci_egg,)
        ret
    end
    inner_A
end


"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct PlanarSphericalJoint{valueType,e2eType,axesType,orderType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::valueType
    e2e::e2eType
    axes::axesType
    order::orderType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function PlanarSphericalJoint(id,e2e)
    hen = e2e.end1
    egg = e2e.end2
    nΦ = 1
    T = get_numbertype(hen.rbsig)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    axes = SpatialFrame(ā_hen)
    q_hen = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R)
    q_egg = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R)
    X_hen = NCF.make_X(nmcs_hen,q_hen)
    X_egg = NCF.make_X(nmcs_egg,q_egg)
    a_hen = X_hen*ā_hen
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    values = zeros(T,nΦ)
    values[1] = transpose(a_hen)*(C_hen*q_hen.-C_egg*q_egg)
    PlanarSphericalJoint(id,nΦ,values,e2e,axes,order)
end

function make_Φ(cst::PlanarSphericalJoint,indexed,numbered)
    (;nconstraints,values,e2e,axes,order) = cst
    (;num2sys,mem2num) = numbered
    (;mem2sysfull) = indexed
    hen = e2e.end1
    egg = e2e.end2
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    function inner_Φ(q)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        ret[1] = transpose(a_hen)*(C_hen*q_hen.-C_egg*q_egg) .- values[1]
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)*q_hen
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)*q_egg
        ret[1] = transpose(a_hen)*(C_hen*q_hen.-C_egg*q_egg) .- values[1]
        ret
    end
    inner_Φ
end

function make_A(cst::PlanarSphericalJoint,indexed,numbered)
    (;nconstraints,e2e,axes,order) = cst
    (;mem2sysfree,mem2sysfull,nfree) = indexed
    (;num2sys,mem2num) = numbered
    hen = e2e.end1
    egg = e2e.end2
    uci_hen =  hen.rbsig.state.cache.free_idx
    uci_egg =  egg.rbsig.state.cache.free_idx
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    free_hen = mem2sysfree[id_hen]
    free_egg = mem2sysfree[id_egg]
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        ret[1,free_hen] .= kron(
                [
                    0 transpose(ā_hen); 
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1,free_hen] .+= transpose(a_hen)*C_hen[:,uci_hen]
        ret[1,free_egg] .-= transpose(a_hen)*C_egg[:,uci_egg]
        ret
    end
    function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)
        ret[1,free_hen] .= kron(
                [
                    0 transpose(ā_hen); 
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1,free_hen] .+= transpose(a_hen)*C_hen[:,uci_hen]
        ret[1,free_egg] .-= transpose(a_hen)*C_egg[:,uci_egg]
        ret
    end
    inner_A
end


"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct PlanarJoint{valueType,e2eType,axesType,orderType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::valueType
    e2e::e2eType
    axes::axesType
    order::orderType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function PlanarJoint(id,e2e)
    hen = e2e.end1
    egg = e2e.end2
    nΦ = 4
    T = get_numbertype(hen.rbsig)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    axes = SpatialFrame(ā_hen)
    q_hen = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R)
    q_egg = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R)
    X_hen = NCF.make_X(nmcs_hen,q_hen)
    X_egg = NCF.make_X(nmcs_egg,q_egg)
    a_hen = X_hen*ā_hen
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    values = zeros(T,nΦ)
    values[1] = transpose(a_hen)*(C_hen*q_hen.-C_egg*q_egg)
    nΦ1 = NCF.get_nconstraints(nmcs_hen)
    nΦ2 = NCF.get_nconstraints(nmcs_egg)
    inprods = transpose(X_egg)*X_hen
    order = find_order(nmcs_hen,nmcs_egg,q_hen,q_egg,X_hen,X_egg,nΦ1,nΦ2)
    values[2:4] .= inprods[order]
    PlanarJoint(id,nΦ,values,e2e,axes,order)
end

function make_Φ(cst::PlanarJoint,indexed,numbered)
    (;nconstraints,values,e2e,axes,order) = cst
    (;num2sys,mem2num) = numbered
    (;mem2sysfull) = indexed
    hen = e2e.end1
    egg = e2e.end2
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    function inner_Φ(q)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        ret[1] = transpose(a_hen)*(C_hen*q_hen.-C_egg*q_egg) .- values[1]
        inprods = transpose(X_egg)*X_hen
        ret[2:4] .= inprods[order] .- values[2:4]
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)*q_hen
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)*q_egg
        ret[1] = transpose(a_hen)*(C_hen*q_hen.-C_egg*q_egg) .- values[1]
        inprods = transpose(X_egg)*X_hen
        ret[2:4] .= inprods[order] .- values[2:4]
        ret
    end
    inner_Φ
end

function make_A(cst::PlanarJoint,indexed,numbered)
    (;nconstraints,e2e,axes,order) = cst
    (;mem2sysfree,mem2sysfull,nfree) = indexed
    (;num2sys,mem2num) = numbered
    hen = e2e.end1
    egg = e2e.end2
    uci_hen =  hen.rbsig.state.cache.free_idx
    uci_egg =  egg.rbsig.state.cache.free_idx
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    free_hen = mem2sysfree[id_hen]
    free_egg = mem2sysfree[id_egg]
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        ret[1,free_hen] .= kron(
                [
                    0 transpose(ā_hen); 
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1,free_hen] .+= transpose(a_hen)*C_hen[:,uci_hen]
        ret[1,free_egg] .-= transpose(a_hen)*C_egg[:,uci_egg]
        k = 1
        ret_rest = @view ret[k+1:end,:]
        rot_jac!(ret_rest,order,q_hen,q_egg,X_hen,X_egg,free_hen,free_egg,uci_hen,uci_egg,)
        ret
    end
    function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)
        ret[1,free_hen] .= kron(
                [
                    0 transpose(ā_hen); 
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1,free_hen] .+= transpose(a_hen)*C_hen[:,uci_hen]
        ret[1,free_egg] .-= transpose(a_hen)*C_egg[:,uci_egg]
        k = 1
        ret_rest = @view ret[k+1:end,:]
        rot_jac!(ret_rest,order,q_hen,q_egg,X_hen,X_egg,free_hen,free_egg,uci_hen,uci_egg,)
        ret
    end
    inner_A
end


"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct PlanarRevoluteJoint{valueType,e2eType,axesType,maskType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::valueType
    e2e::e2eType
    axes::axesType
    mask::maskType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function PlanarRevoluteJoint(id,e2e)
    hen = e2e.end1
    egg = e2e.end2
    nΦ = 3
    T = get_numbertype(hen.rbsig)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    axes = SpatialFrame(ā_hen)
    q_hen = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R)
    q_egg = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R)
    X_hen = NCF.make_X(nmcs_hen,q_hen)
    X_egg = NCF.make_X(nmcs_egg,q_egg)
    a_hen = X_hen*ā_hen
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    values = zeros(T,nΦ)
    values[1] = transpose(a_hen)*(C_hen*q_hen.-C_egg*q_egg)
    inprods = transpose(X_egg)*a_hen
    mask = 1:3 .!== argmax(abs.(inprods))
    values[2:3] = inprods[mask]
    PlanarRevoluteJoint(id,nΦ,values,e2e,axes,mask)
end

function make_Φ(cst::PlanarRevoluteJoint,indexed,numbered)
    (;nconstraints,values,e2e,axes,mask) = cst
    (;num2sys,mem2num) = numbered
    (;mem2sysfull) = indexed
    hen = e2e.end1
    egg = e2e.end2
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    function inner_Φ(q)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        ret[1] = transpose(a_hen)*(C_hen*q_hen.-C_egg*q_egg) .- values[1]
        ret[2:3] = transpose(X_egg[:,mask])*a_hen - values[2:3]
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)*q_hen
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)*q_egg
        ret[1] = transpose(a_hen)*(C_hen*q_hen.-C_egg*q_egg) .- values[1]
        ret[2:3] = transpose(X_egg[:,mask])*a_hen - values[2:3]
        ret
    end
    inner_Φ
end

function make_A(cst::PlanarRevoluteJoint,indexed,numbered)
    (;nconstraints,e2e,axes,order) = cst
    (;mem2sysfree,mem2sysfull,nfree) = indexed
    (;num2sys,mem2num) = numbered
    hen = e2e.end1
    egg = e2e.end2
    uci_hen =  hen.rbsig.state.cache.free_idx
    uci_egg =  egg.rbsig.state.cache.free_idx
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    free_hen = mem2sysfree[id_hen]
    free_egg = mem2sysfree[id_egg]
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        ret[1,free_hen] .= kron(
                [
                    0 transpose(ā_hen); 
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1,free_hen] .+= transpose(a_hen)*C_hen[:,uci_hen]
        ret[1,free_egg] .-= transpose(a_hen)*C_egg[:,uci_egg]
        o3 = zero(a_hen)
        ret[2:3,free_hen] .= (transpose(
            kron(vcat(0,ā_hen),X_egg)[:,mask])*Y_hen)[:,uci_hen]
        ret[2:3,free_egg] .= (transpose(
            [
                o3 o3 o3;
                a_hen o3 o3; 
                o3 a_hen o3;  
                o3 o3 a_hen; 
            ][:,mask])*Y_egg)[:,uci_egg]
        ret
        ret
    end
    function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)
        ret[1,free_hen] .= kron(
                [
                    0 transpose(ā_hen); 
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1,free_hen] .+= transpose(a_hen)*C_hen[:,uci_hen]
        ret[1,free_egg] .-= transpose(a_hen)*C_egg[:,uci_egg]
        o3 = zero(a_hen)
        ret[2:3,free_hen] .= (transpose(
            kron(vcat(0,ā_hen),X_egg)[:,mask])*Y_hen)[:,uci_hen]
        ret[2:3,free_egg] .= (transpose(
            [
                o3 o3 o3;
                a_hen o3 o3;
                o3 a_hen o3;  
                o3 o3 a_hen;
            ][:,mask])*Y_egg)[:,uci_egg]
        ret
    end
    inner_A
end

"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct CylindricalJoint{valueType,e2eType,axesType,maskType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::valueType
    e2e::e2eType
    axes::axesType
    mask::maskType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function CylindricalJoint(id,e2e)
    hen = e2e.end1
    egg = e2e.end2
    nΦ = 4
    T = get_numbertype(hen.rbsig)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    axes = SpatialFrame(ā_hen)
    q_hen = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R)
    q_egg = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R)
    X_hen = NCF.make_X(nmcs_hen,q_hen)
    X_egg = NCF.make_X(nmcs_egg,q_egg)
    a1t1 = X_hen*axes.t1
    a1t2 = X_hen*axes.t2
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    values = zeros(T,nΦ)    
    values[1:2] .= transpose([a1t1 a1t2])*(C_hen*q_hen.-C_egg*q_egg)
    a_hen = X_hen*ā_hen
    inprods = transpose(X_egg)*a_hen
    mask = 1:3 .!== argmax(abs.(inprods))
    values[3:4] .= inprods[mask]
    CylindricalJoint(id,nΦ,values,e2e,axes,mask)
end

function make_Φ(cst::CylindricalJoint,indexed,numbered)
    (;nconstraints,values,e2e,axes,mask) = cst
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
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        a1t1 = X_hen*axes.t1
        a1t2 = X_hen*axes.t2
        ret[1:2] .= transpose([a1t1 a1t2])*(C_hen*q_hen.-C_egg*q_egg) .- values[1:2]
        ret[3:4] .= transpose(X_egg[:,mask])*a_hen - values[4:5]
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a1t1 = X_hen*axes.t1
        a1t2 = X_hen*axes.t2
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)*q_hen
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)*q_egg
        ret[1:2] .= transpose([a1t1 a1t2])*(C_hen*q_hen.-C_egg*q_egg) .- values[1:2]
        ret[3:4] .= transpose(X_egg[:,mask])*a_hen - values[4:5]
        ret
    end
    inner_Φ
end

function make_A(cst::CylindricalJoint,indexed,numbered)
    (;nconstraints,e2e,axes,mask) = cst
    (;mem2sysfree,mem2sysfull,nfree) = indexed
    (;num2sys,mem2num) = numbered
    hen = e2e.end1
    egg = e2e.end2
    uci_hen =  hen.rbsig.state.cache.free_idx
    uci_egg =  egg.rbsig.state.cache.free_idx
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    ā_hen = hen.rbsig.prop.ās[hen.aid]
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
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a1t1 = X_hen*axes.t1
        a1t2 = X_hen*axes.t2
        a_hen = X_hen*ā_hen
        ret[1:2,free_hen] .= kron(
                [
                    0 transpose(axes.t1); 
                    0 transpose(axes.t2)
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1:2,free_hen] .+= transpose([a1t1 a1t2])*C_hen[:,uci_hen]
        ret[1:2,free_egg] .-= transpose([a1t1 a1t2])*C_egg[:,uci_egg]
        o3 = zero(a_hen)
        ret[3:4,mem2sysfree[id_hen]] .= (transpose(
            kron(vcat(0,ā_hen),X_egg)[:,mask])*Y_hen)[:,uci_hen]
        ret[3:4,mem2sysfree[id_egg]] .= (transpose(
            [
                o3 o3 o3;
                a_hen o3 o3; 
                o3 a_hen o3;  
                o3 o3 a_hen; 
            ][:,mask])*Y_egg)[:,uci_egg]
        ret
    end
    function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a1t1 = X_hen*axes.t1
        a1t2 = X_hen*axes.t2
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)
        ret[1:2,free_hen] .= kron(
                [
                    0 transpose(axes.t1); 
                    0 transpose(axes.t2)
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1:2,free_hen] .+= transpose([a1t1 a1t2])*C_hen[:,uci_hen]
        ret[1:2,free_egg] .-= transpose([a1t1 a1t2])*C_egg[:,uci_egg]
        a_hen = X_hen*ā_hen
        o3 = zero(a_hen)
        ret[3:4,mem2sysfree[id_hen]] .= (transpose(
            kron(vcat(0,ā_hen),X_egg)[:,mask])*Y_hen)[:,uci_hen]
        ret[3:4,mem2sysfree[id_egg]] .= (transpose(
            [
                o3 o3 o3;
                a_hen o3 o3; 
                o3 a_hen o3;  
                o3 o3 a_hen; 
            ][:,mask])*Y_egg)[:,uci_egg]
        ret
    end
    inner_A
end

"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct UniversalJoint{valueType,axesType,e2eType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::valueType
    axes_trl_hen::axesType
    axes_rot_hen::axesType
    axes_rot_egg::axesType
    e2e::e2eType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function UniversalJoint(id,e2e)
    hen = e2e.end1
    egg = e2e.end2
    nΦ = 4
    T = get_numbertype(hen.rbsig)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    q_hen = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R)
    q_egg = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R)
    X_hen = NCF.make_X(nmcs_hen,q_hen)
    X_egg = NCF.make_X(nmcs_egg,q_egg)
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    # rotate of egg
    ā_egg = egg.rbsig.prop.ās[egg.aid]
    axes_rot_egg = Axes(ā_egg)
    rot_hen_n = X_egg*axes_rot_egg.t1
    rot_egg_t2 = X_egg*axes_rot_egg.t2
    rot_hen_n = inv(X_hen)*rot_hen_n
    axes_rot_hen = Axes(rot_hen_n)
    # constraint values
    values = zeros(T,nΦ)
    values[1:3] .= C_hen*q_hen.-C_egg*q_egg
    values[4] = rot_hen_n'*rot_egg_t2
    UniversalJoint(id,nΦ,values,axes_rot_hen,axes_rot_hen,axes_rot_egg,e2e)
end

function make_Φ(cst::UniversalJoint,indexed,numbered)
    (;nconstraints,values,e2e) = cst
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
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        # rotate of egg
        rot_hen_n = X_hen*axes_rot_hen.n
        rot_egg_t2 = X_egg*axes_rot_egg.t2
        ret[1:3] .= (C_hen*q_hen.-C_egg*q_egg) .- values[1:3]
        ret[4] = rot_hen_n'*rot_egg_t2 - values[4]
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)*q_hen
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)*q_egg
        # rotate of egg
        rot_hen_n = X_hen*axes_rot_hen.n
        rot_egg_t2 = X_egg*axes_rot_egg.t2
        ret[1:3] .= C_hen*q_hen.-C_egg*q_egg .- values[1:3]
        ret[4] = rot_hen_n'*rot_egg_t2 - values[4]
        ret
    end
    inner_Φ
end

function make_A(cst::UniversalJoint,indexed,numbered)
    (;nconstraints,axes_rot_hen,axes_rot_egg,e2e) = cst
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
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        # rotate of egg
        rot_hen_n = X_hen*axes_rot_hen.n
        rot_egg_t2 = X_egg*axes_rot_egg.t2
        ret[1:3,free_hen] .=  C_hen[:,uci_hen]
        ret[1:3,free_egg] .= -C_egg[:,uci_egg]
        ret[4,free_hen] = (
                transpose(
                    kron(vcat(0,axes_rot_hen.n),rot_egg_t2)
                )*Y_hen
            )[:,uci_hen]
        ret[4,free_egg] = (
                transpose(
                    kron(vcat(0,axes_rot_egg.t2),rot_hen_n)
                )*Y_egg
            )[:,uci_egg]
        ret
    end
    function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        c_hen = c[num2sys[mem2num[id_hen][hen.pid]]]
        c_egg = c[num2sys[mem2num[id_egg][egg.pid]]]
        # rotate of egg
        rot_hen_n = X_hen*axes_rot_hen.n
        rot_egg_t2 = X_egg*axes_rot_egg.t2
        ret[1:3,free_hen] .=  hen.rbsig.state.cache.funcs.C(c_hen)[:,uci_hen]
        ret[1:3,free_egg] .= -egg.rbsig.state.cache.funcs.C(c_egg)[:,uci_egg]
        ret[4,free_hen] = (
                transpose(
                    kron(vcat(0,axes_rot_hen.n),rot_egg_t2)
                )*Y_hen
            )[:,uci_hen]
        ret[4,free_egg] = (
                transpose(
                    kron(vcat(0,axes_rot_egg.t2),rot_hen_n)
                )*Y_egg
            )[:,uci_egg]
        ret
    end
    inner_A
end


"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct PrismaticUniversalJoint{valueType,e2eType,axesType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::valueType
    axes_trl_hen::axesType
    axes_rot_hen::axesType
    axes_rot_egg::axesType
    e2e::e2eType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function PrismaticUniversalJoint(id,e2e)
    hen = e2e.end1
    egg = e2e.end2
    nΦ = 3
    T = get_numbertype(hen.rbsig)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    q_hen = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R)
    q_egg = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R)
    X_hen = NCF.make_X(nmcs_hen,q_hen)
    X_egg = NCF.make_X(nmcs_egg,q_egg)
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    # translate on hen
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    axes_trl_hen = Axes(ā_hen)
    trl_hen_t1 = X_hen*axes_trl_hen.t1
    trl_hen_t2 = X_hen*axes_trl_hen.t2
    # rotate of egg
    ā_egg = egg.rbsig.prop.ās[egg.aid]
    axes_rot_egg = Axes(ā_egg)
    rot_hen_n = X_egg*axes_rot_egg.t1
    rot_egg_t2 = X_egg*axes_rot_egg.t2
    rot_hen_n = inv(X_hen)*rot_hen_n
    axes_rot_hen = Axes(rot_hen_n)
    # constraint values
    values = zeros(T,nΦ)
    values[1:2] .= transpose([trl_hen_t1 trl_hen_t2])*(C_hen*q_hen.-C_egg*q_egg)
    values[3] = rot_hen_n'*rot_egg_t2
    PrismaticUniversalJoint(id,nΦ,values,axes_trl_hen,axes_rot_hen,axes_rot_egg,e2e)
end

function make_Φ(cst::PrismaticUniversalJoint,indexed,numbered)
    (;nconstraints,values,axes_trl_hen,axes_rot_hen,axes_rot_egg,e2e) = cst
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
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        # translate on hen
        trl_hen_t1 = X_hen*axes_trl_hen.t1
        trl_hen_t2 = X_hen*axes_trl_hen.t2
        # rotate of egg
        rot_hen_n = X_hen*axes_rot_hen.n
        rot_egg_t2 = X_egg*axes_rot_egg.t2
        # constraint values
        ret[1:2] .= transpose([trl_hen_t1 trl_hen_t2])*(C_hen*q_hen.-C_egg*q_egg) .- values[1:2]
        ret[3] = rot_hen_n'*rot_egg_t2 - values[3]
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)*q_hen
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)*q_egg
        # translate on hen
        trl_hen_t1 = X_hen*axes_trl_hen.t1
        trl_hen_t2 = X_hen*axes_trl_hen.t2
        # rotate of egg
        rot_hen_n = X_hen*axes_rot_hen.n
        rot_egg_t2 = X_egg*axes_rot_egg.t2
        # constraint values
        ret[1:2] .= transpose([trl_hen_t1 trl_hen_t2])*(C_hen*q_hen.-C_egg*q_egg) .- values[1:2]
        ret[3] = rot_hen_n'*rot_egg_t2 - values[3]
        ret
    end
    inner_Φ
end

function make_A(cst::PrismaticUniversalJoint,indexed,numbered)
    (;nconstraints,axes_trl_hen,axes_rot_hen,axes_rot_egg,e2e) = cst
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
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        # translate on hen
        trl_hen_t1 = X_hen*axes_trl_hen.t1
        trl_hen_t2 = X_hen*axes_trl_hen.t2
        # rotate of egg
        rot_hen_n = X_hen*axes_rot_hen.n
        rot_egg_t2 = X_egg*axes_rot_egg.t2
        # jacobian
        ret[1:2,free_hen] .= kron(
                [
                    0 transpose(trl_hen_t1); 
                    0 transpose(trl_hen_t2)
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1:2,free_hen] .+= transpose([trl_hen_t1 trl_hen_t2])*C_hen[:,uci_hen]
        ret[1:2,free_egg] .-= transpose([trl_hen_t1 trl_hen_t2])*C_egg[:,uci_egg]
        ret[3,free_hen] = (
                transpose(
                    kron(vcat(0,axes_rot_hen.n),rot_egg_t2)
                )*Y_hen
            )[:,uci_hen]
        ret[3,free_egg] = (
                transpose(
                    kron(vcat(0,axes_rot_egg.t2),rot_hen_n)
                )*Y_egg
            )[:,uci_egg]
        ret
    end
    function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        c_hen = c[num2sys[mem2num[id_hen][hen.pid]]]
        c_egg = c[num2sys[mem2num[id_egg][egg.pid]]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)
        # translate on hen
        trl_hen_t1 = X_hen*axes_trl_hen.t1
        trl_hen_t2 = X_hen*axes_trl_hen.t2
        # rotate of egg
        rot_hen_n = X_hen*axes_rot_hen.n
        rot_egg_t2 = X_egg*axes_rot_egg.t2
        # jacobian
        ret[1:2,free_hen] .= kron(
                [
                    0 transpose(trl_hen_t1); 
                    0 transpose(trl_hen_t2)
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1:2,free_hen] .+= transpose([trl_hen_t1 trl_hen_t2])*C_hen[:,uci_hen]
        ret[1:2,free_egg] .-= transpose([trl_hen_t1 trl_hen_t2])*C_egg[:,uci_egg]
        ret[3,free_hen] = (
                transpose(
                    kron(vcat(0,axes_rot_hen.n),rot_egg_t2)
                )*Y_hen
            )[:,uci_hen]
        ret[3,free_egg] = (
                transpose(
                    kron(vcat(0,axes_rot_egg.t2),rot_hen_n)
                )*Y_egg
            )[:,uci_egg]
        ret
    end
    inner_A
end

"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct PlanarUniversalJoint{valueType,e2eType,axesType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::valueType
    axes::axesType
    e2e::e2eType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function PlanarUniversalJoint(id,e2e)
    hen = e2e.end1
    egg = e2e.end2
    nΦ = 2
    T = get_numbertype(hen.rbsig)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    ā_egg = egg.rbsig.prop.ās[egg.aid]
    q_hen = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R)
    q_egg = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R)
    X_hen = NCF.make_X(nmcs_hen,q_hen)
    X_egg = NCF.make_X(nmcs_egg,q_egg)
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    values = zeros(T,nΦ)
    values[1] = transpose(a_hen)*(C_hen*q_hen.-C_egg*q_egg)
    a_hen = X_hen*ā_hen
    a_egg = X_egg*ā_egg
    values[2] = a_egg'*a_hen
    PlanarUniversalJoint(id,nΦ,values,e2e)
end

function make_Φ(cst::PlanarUniversalJoint,indexed,numbered)
    (;nconstraints,values,e2e) = cst
    (;num2sys,mem2num) = numbered
    (;mem2sysfull) = indexed
    hen = e2e.end1
    egg = e2e.end2
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    ā_egg = egg.rbsig.prop.ās[egg.aid]
    function inner_Φ(q)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        a_egg = X_egg*ā_egg
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        ret[1] = transpose(a_hen)*(C_hen*q_hen.-C_egg*q_egg) .- values[1]
        ret[2] = a_egg'*a_hen - values[2]
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        a_egg = X_egg*ā_egg
        c_hen = c[num2sys[mem2num[id_hen]][hen.pid]]
        c_egg = c[num2sys[mem2num[id_egg]][egg.pid]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)*q_hen
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)*q_egg
        ret[1] = transpose(a_hen)*(C_hen*q_hen.-C_egg*q_egg) .- values[1]
        ret[2] = a_egg'*a_hen - values[2]
        ret
    end
    inner_Φ
end

function make_A(cst::PlanarUniversalJoint,indexed,numbered)
    (;nconstraints,e2e) = cst
    (;mem2sysfree,mem2sysfull,nfree) = indexed
    (;num2sys,mem2num) = numbered
    hen = e2e.end1
    egg = e2e.end2
    uci_hen =  hen.rbsig.state.cache.free_idx
    uci_egg =  egg.rbsig.state.cache.free_idx
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    ā_egg = egg.rbsig.prop.ās[egg.aid]
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
        C_hen = hen.rbsig.state.cache.Cps[hen.pid]
        C_egg = egg.rbsig.state.cache.Cps[egg.pid]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        a_egg = X_egg*ā_egg
        ret[1,free_hen] .= kron(
                [
                    0 transpose(ā_hen); 
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1,free_hen] .+= transpose(a_hen)*C_hen[:,uci_hen]
        ret[1,free_egg] .-= transpose(a_hen)*C_egg[:,uci_egg]
        ret[2,free_hen] = (
                transpose(
                    kron(vcat(0,ā_hen),a_egg)
                )*Y_hen
            )[:,uci_hen]
        ret[2,free_egg] = (
                transpose(
                    kron(vcat(0,ā_egg),a_hen)
                )*Y_egg
            )[:,uci_egg]
        ret
    end
    function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a1t1 = X_hen*axes.t1
        a1t2 = X_hen*axes.t2
        c_hen = c[num2sys[mem2num[id_hen][hen.pid]]]
        c_egg = c[num2sys[mem2num[id_egg][egg.pid]]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)
        ret[1,free_hen] .= kron(
                [
                    0 transpose(ā_hen); 
                ],
                transpose(C_hen*q_hen.-C_egg*q_egg)
            )[:,uci_hen]
        ret[1,free_hen] .+= transpose(a_hen)*C_hen[:,uci_hen]
        ret[1,free_egg] .-= transpose(a_hen)*C_egg[:,uci_egg]
        a_hen = X_hen*ā_hen
        a_egg = X_egg*ā_egg
        ret[2,mem2sysfree[id_hen]] = (
                transpose(
                    kron(vcat(0,ā_hen),a_egg)
                )*Y_hen
            )[:,uci_hen]
        ret[2,mem2sysfree[id_egg]] = (
                transpose(
                    kron(vcat(0,ā_egg),a_hen)
                )*Y_egg
            )[:,uci_egg]
        ret
    end
    inner_A
end


"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct FloatingUniversalJoint{valueType,e2eType,axesType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::valueType
    axes::axesType
    e2e::e2eType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function FloatingUniversalJoint(id,e2e)
    hen = e2e.end1
    egg = e2e.end2
    nΦ = 1
    T = get_numbertype(hen.rbsig)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    ā_egg = egg.rbsig.prop.ās[egg.aid]
    q_hen = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R)
    q_egg = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R)
    X_hen = NCF.make_X(nmcs_hen,q_hen)
    X_egg = NCF.make_X(nmcs_egg,q_egg)
    values = zeros(T,nΦ)
    a_hen = X_hen*ā_hen
    a_egg = X_egg*ā_egg
    values[1] = a_egg'*a_hen
    FloatingUniversalJoint(id,nΦ,values,e2e)
end

function make_Φ(cst::FloatingUniversalJoint,indexed,numbered)
    (;nconstraints,values,e2e) = cst
    (;num2sys,mem2num) = numbered
    (;mem2sysfull) = indexed
    hen = e2e.end1
    egg = e2e.end2
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    ā_egg = egg.rbsig.prop.ās[egg.aid]
    function inner_Φ(q)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        a_egg = X_egg*ā_egg
        ret[1] = a_egg'*a_hen - values[1]
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        a_egg = X_egg*ā_egg
        ret[1] = a_egg'*a_hen - values[1]
        ret
    end
    inner_Φ
end

function make_A(cst::FloatingUniversalJoint,indexed,numbered)
    (;nconstraints,e2e) = cst
    (;mem2sysfree,mem2sysfull,nfree) = indexed
    (;num2sys,mem2num) = numbered
    hen = e2e.end1
    egg = e2e.end2
    uci_hen =  hen.rbsig.state.cache.free_idx
    uci_egg =  egg.rbsig.state.cache.free_idx
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    ā_hen = hen.rbsig.prop.ās[hen.aid]
    ā_egg = egg.rbsig.prop.ās[egg.aid]
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
        a_hen = X_hen*ā_hen
        a_egg = X_egg*ā_egg
        ret[1,free_hen] = (
                transpose(
                    kron(vcat(0,ā_hen),a_egg)
                )*Y_hen
            )[:,uci_hen]
        ret[1,free_egg] = (
                transpose(
                    kron(vcat(0,ā_egg),a_hen)
                )*Y_egg
            )[:,uci_egg]
        ret
    end
    function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        a_hen = X_hen*ā_hen
        a_egg = X_egg*ā_egg
        ret[1,mem2sysfree[id_hen]] = (
                transpose(
                    kron(vcat(0,ā_hen),a_egg)
                )*Y_hen
            )[:,uci_hen]
        ret[1,mem2sysfree[id_egg]] = (
                transpose(
                    kron(vcat(0,ā_egg),a_hen)
                )*Y_egg
            )[:,uci_egg]
        ret
    end
    inner_A
end

"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct FloatingJoint{valueType,e2eType,orderType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::valueType
    e2e::e2eType
    order::orderType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function FloatingJoint(id,e2e)
    hen = e2e.end1
    egg = e2e.end2
    rb1 = hen.rbsig
    nΦ = 3
    T = get_numbertype(rb1)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    q_hen = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R)
    q_egg = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R)
    X_hen = NCF.make_X(nmcs_hen,q_hen)
    X_egg = NCF.make_X(nmcs_egg,q_egg)
    values = zeros(T,nΦ)
    nΦ1 = NCF.get_nconstraints(nmcs_hen)
    nΦ2 = NCF.get_nconstraints(nmcs_egg)
    order = find_order(nmcs_hen,nmcs_egg,q_hen,q_egg,X_hen,X_egg,nΦ1,nΦ2)
    inprods = transpose(X_egg)*X_hen
    values[1:3] .= inprods[order]
    FloatingJoint(id,nΦ,values,e2e,order)
end

function make_Φ(cst::FloatingJoint,indexed,numbered)
    (;nconstraints,values,e2e,order) = cst
    (;mem2num,num2sys) = numbered
    (;mem2sysfull) = indexed
    hen = e2e.end1
    egg = e2e.end2
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    function inner_Φ(q)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[hen.rbsig.prop.id]]
        q_egg = @view q[mem2sysfull[egg.rbsig.prop.id]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        inprods = transpose(X_egg)*X_hen
        ret[1:3] .= inprods[order] .- values[1:3]
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        id_hen = hen.rbsig.prop.id
        id_egg = egg.rbsig.prop.id
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        inprods = transpose(X_egg)*X_hen
        ret[1:3] .= inprods[order] .- values[1:3]
        ret
    end
    inner_Φ
end

function make_A(cst::FloatingJoint,indexed,numbered)
    (;nconstraints,e2e,order) = cst
    (;mem2sysfree,mem2sysfull,nfree) = indexed
    (;mem2num,num2sys) = numbered
    hen = e2e.end1
    egg = e2e.end2
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    uci_hen =  hen.rbsig.state.cache.free_idx
    uci_egg =  egg.rbsig.state.cache.free_idx
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    free_hen = mem2sysfree[id_hen]
    free_egg = mem2sysfree[id_egg]
    function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        k = 0
        ret_rest = @view ret[k+1:end,:]
        rot_jac!(ret_rest,order,q_hen,q_egg,X_hen,X_egg,free_hen,free_egg,uci_hen,uci_egg,)
        ret
    end
    function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        k = 0
        ret_rest = @view ret[k+1:end,:]
        rot_jac!(ret_rest,order,q_hen,q_egg,X_hen,X_egg,free_hen,free_egg,uci_hen,uci_egg,)
        ret
    end
    inner_A
end

"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct FixedJoint{valueType,e2eType,orderType} <: ExternalConstraints{valueType}
    id::Int
    nconstraints::Int
    values::valueType
    e2e::e2eType
    order::orderType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function FixedJoint(id,e2e)
    hen = e2e.end1
    egg = e2e.end2
    rb1 = hen.rbsig
    nΦ = 6
    T = get_numbertype(rb1)
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    state_hen = hen.rbsig.state
    state_egg = egg.rbsig.state
    q_hen = NCF.rigidstate2naturalcoords(nmcs_hen,state_hen.ro,state_hen.R)
    q_egg = NCF.rigidstate2naturalcoords(nmcs_egg,state_egg.ro,state_egg.R)
    X_hen = NCF.make_X(nmcs_hen,q_hen)
    X_egg = NCF.make_X(nmcs_egg,q_egg)
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    values = zeros(T,nΦ)
    values[1:3] .= C_hen*q_hen.-C_egg*q_egg
    nΦ1 = NCF.get_nconstraints(nmcs_hen)
    nΦ2 = NCF.get_nconstraints(nmcs_egg)
    order = find_order(nmcs_hen,nmcs_egg,q_hen,q_egg,X_hen,X_egg,nΦ1,nΦ2)
    inprods = transpose(X_egg)*X_hen
    values[4:6] .= inprods[order]
    FixedJoint(id,nΦ,values,e2e,order)
end

function make_Φ(cst::FixedJoint,indexed,numbered)
    (;nconstraints,values,e2e,order) = cst
    (;mem2num,num2sys) = numbered
    (;mem2sysfull) = indexed
    hen = e2e.end1
    egg = e2e.end2
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    function inner_Φ(q)
        ret = zeros(eltype(q),nconstraints)
        q_hen = @view q[mem2sysfull[hen.rbsig.prop.id]]
        q_egg = @view q[mem2sysfull[egg.rbsig.prop.id]]
        ret[1:3] .= C_hen*q_hen.-C_egg*q_egg .- values[1:3]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        inprods = transpose(X_egg)*X_hen
        ret[4:6] .= inprods[order] .- values[4:6]
        ret
    end
    function inner_Φ(q,d,c)
        ret = zeros(eltype(q),nconstraints)
        id_hen = hen.rbsig.prop.id
        id_egg = egg.rbsig.prop.id
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        c_hen = c[num2sys[mem2num[id_hen][hen.pid]]]
        c_egg = c[num2sys[mem2num[id_egg][egg.pid]]]
        C_hen = hen.rbsig.state.cache.funcs.C(c_hen)*q_hen
        C_egg = egg.rbsig.state.cache.funcs.C(c_egg)*q_egg
        ret[1:3] .= C_hen*q_hen.-C_egg*q_egg .- values[1:3]
        inprods = transpose(X_egg)*X_hen
        ret[4:6] .= inprods[order] .- values[4:6]
        ret
    end
    inner_Φ
end

function make_A(cst::FixedJoint,indexed,numbered)
    (;nconstraints,e2e,order) = cst
    (;mem2sysfree,mem2sysfull,nfree) = indexed
    (;mem2num,num2sys) = numbered
    hen = e2e.end1
    egg = e2e.end2
    nmcs_hen = hen.rbsig.state.cache.funcs.nmcs
    nmcs_egg = egg.rbsig.state.cache.funcs.nmcs
    C_hen = hen.rbsig.state.cache.Cps[hen.pid]
    C_egg = egg.rbsig.state.cache.Cps[egg.pid]
    uci_hen =  hen.rbsig.state.cache.free_idx
    uci_egg =  egg.rbsig.state.cache.free_idx
    id_hen = hen.rbsig.prop.id
    id_egg = egg.rbsig.prop.id
    free_hen = mem2sysfree[id_hen]
    free_egg = mem2sysfree[id_egg]
    function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        ret[1:3,free_hen] .=  C_hen[:,uci_hen]
        ret[1:3,free_egg] .= -C_egg[:,uci_egg]
        k = 3
        ret_rest = @view ret[k+1:end,:]
        rot_jac!(ret_rest,order,q_hen,q_egg,X_hen,X_egg,free_hen,free_egg,uci_hen,uci_egg,)
        ret
    end
    function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
        q_hen = @view q[mem2sysfull[id_hen]]
        q_egg = @view q[mem2sysfull[id_egg]]
        X_hen = NCF.make_X(nmcs_hen,q_hen)
        X_egg = NCF.make_X(nmcs_egg,q_egg)
        c_hen = c[num2sys[mem2num[id_hen][hen.pid]]]
        c_egg = c[num2sys[mem2num[id_egg][egg.pid]]]
        ret[1:3,free_hen] .=  hen.rbsig.state.cache.funcs.C(c_hen)[:,uci_hen]
        ret[1:3,free_egg] .= -egg.rbsig.state.cache.funcs.C(c_egg)[:,uci_egg]
        k = 3
        ret_rest = @view ret[k+1:end,:]
        rot_jac!(ret_rest,order,q_hen,q_egg,X_hen,X_egg,free_hen,free_egg,uci_hen,uci_egg,)
        ret
    end
    inner_A
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

function get_jointed_free_idx(cst::LinearJoint)
    uci = collect(1:size(cst.A,2))
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

function get_jointed_free(cst::LinearJoint,indexed)
    cst2sysfree = collect(1:indexed.nfree)
end

function make_Φqᵀq(cst,numbered,c)
    (;axes_trl_hen,axes_rot_hen,axes_rot_egg,e2e,) = cst
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
    c_hen = c[num2sys[mem2num[id_hen][hen.pid]]]
    c_egg = c[num2sys[mem2num[id_egg][egg.pid]]]
    ā_hen = axes_rot_hen.n
    ā_egg = axes_rot_egg.t2
    D̄ = axes_trl_hen.X

    # Φqᵀq_first = 0

    Φqᵀq_second = 
        begin
            ret_raw = zeros(
                T,
                (1+nld_hen)+(1+nld_egg),
                (1+nld_hen)+(1+nld_egg),
            )
            ret_raw[2:1+nld_hen,(1+nld_hen)+1+id] = ā_hen*ā_egg'
            ret_raw[(1+nld_hen)+1+id,2:1+nld_hen] = ā_egg*ā_hen'
            transpose(cv)*kron(ret_raw,I_Int)*cv
        end
    
    Φqᵀq_third = 
        begin
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

    Φqᵀq_rot_fix = [
        begin
            ret_raw = zeros(
                T,
                (1+nld_hen)+(1+nld_egg),
                (1+nld_hen)+(1+nld_egg),
            )
            ret_raw[       1+i,(1+nld_hen)+j] = 1
            ret_raw[(1+nld_hen)+j,       1+j] = 1
            transpose(cv)*kron(ret_raw,I_Int)*cv
        end
        for i = 1:3 for j = 1:3
    ]

    if cst isa FixedJoint
        Φqᵀq = vcat(zero.(Φqᵀq_trans),Φqᵀq_rot_fix[cst.order])
    elseif cst isa RevoluteJoint
        Φqᵀq = vcat(zero.(Φqᵀq_trans),Φqᵀq_rot_ā_hen[cst.mask])
    elseif cst isa PrismaticJoint
        Φqᵀq = vcat(Φqᵀq_trans[[2,3]],Φqᵀq_rot_ā_hen)
    elseif cst isa CylindricalJoint
        Φqᵀq = vcat(Φqᵀq_trans[[2,3]],Φqᵀq_rot_ā_hen[cst.mask])
    elseif cst isa UniversalJoint
        Φqᵀq = vcat(zero.(Φqᵀq_trans),Φqᵀq_rot_ā_hen*ā_egg.+Φqᵀq_rot_ā_egg*ā_hen)
    elseif cst isa PrismaticUniversalJoint
        Φqᵀq_rot = sum(
                Φqᵀq_rot_ā_hen[i]*ā_egg[i].+Φqᵀq_rot_ā_egg[i]*ā_hen[i]
                for i = 1:3
        )
        Φqᵀq = vcat(
            Φqᵀq_trans[[2,3]],
            [Φqᵀq_rot]
        )
    elseif cst isa PinJoint
        Φqᵀq = vcat(zero.(Φqᵀq_trans),zero.(Φqᵀq_rot_ā1))
    else
        @error "Unknown Joint"
    end
    Φqᵀq
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

function make_∂Aᵀλ∂q(cst::PinJoint,numbered,c)
    (;nconstraints) = cst
    uci = get_jointed_free_idx(cst)
    function ∂Aᵀλ∂q(λ)
        zeros(eltype(λ),length(uci),length(uci))
    end
end

function make_∂Aᵀλ∂q(cst::LinearJoint,numbered,c)
    (;nconstraints) = cst
    uci = get_jointed_free_idx(cst)
    function ∂Aᵀλ∂q(λ)
        zeros(eltype(λ),length(uci),length(uci))
    end
end

function rot_jac!(ret,order,q_hen,q_egg,X_hen,X_egg,free_hen,free_egg,uci_hen,uci_egg,)
    k = 0
    for i = 1:3
        for j = 1:3
            if 3(i-1)+j in order
                k += 1
                tpl1 = zero(q_hen)
                tpl2 = zero(q_egg)         
                tpl1[3+3(i-1)+1:3+3i] .= X_egg[:,j]
                tpl2[3+3(j-1)+1:3+3j] .= X_hen[:,i]
                ret[k, free_hen] .= tpl1[uci_hen] #todo cv
                ret[k, free_egg] .= tpl2[uci_egg] #todo cv
            end
        end
    end
end

#todo 
#done (joint_type == :FloatingSpherical) && (return 3, 3) #dof=6, ncsts=0, 
#todo (joint_type == :OrbitalSpherical)  && (return 2, 3) #dof=5, ncsts=1, 
#done (joint_type == :PlanarSpherical)   && (return 2, 3) #dof=5, ncsts=1, 
#done (joint_type == :PrismaticSpherical)&& (return 1, 3) #dof=4, ncsts=2, 
#done (joint_type == :Spherical)         && (return 0, 3) #dof=3, ncsts=3,  linear trans
#done (joint_type == :FloatingUniversal) && (return 3, 2) #dof=5, ncsts=1, 
#todo (joint_type == :OrbitalUniversal)  && (return 2, 2) #dof=4, ncsts=2,
#done (joint_type == :PlanarUniversal)   && (return 2, 2) #dof=4, ncsts=2, 
#done (joint_type == :PrismaticUniversal)&& (return 1, 2) #dof=3, ncsts=3, 
#done (joint_type == :Universal)         && (return 0, 2) #dof=2, ncsts=4,  linear trans
#done (joint_type == :FloatingRevolute)  && (return 3, 1) #dof=4, ncsts=2, 
#todo (joint_type == :OrbitalRevolute)   && (return 2, 1) #dof=3, ncsts=3, 
#done (joint_type == :PlanarRevolute)    && (return 2, 1) #dof=3, ncsts=3, 
#done (joint_type == :Cylindrical)       && (return 1, 1) #dof=2, ncsts=4, 
#done (joint_type == :Revolute)          && (return 0, 1) #dof=1, ncsts=5,  linear trans
#done (joint_type == :Floating)          && (return 3, 0) #dof=3, ncsts=3,  rot order
#todo (joint_type == :Orbital)           && (return 2, 0) #dof=2, ncsts=4,  rot order
#done (joint_type == :Planar)            && (return 2, 0) #dof=2, ncsts=4,  rot order
#done (joint_type == :Prismatic)         && (return 1, 0) #dof=1, ncsts=5,  rot order
#done (joint_type == :Fixed)             && (return 0, 0) #dof=0, ncsts=6,  rot order #linear trans

