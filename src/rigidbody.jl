abstract type AbstractRigidBody{N,T} end
abstract type AbstractRigidBodyProperty{N,T} end
abstract type AbstractRigidBodyState{N,T} end
abstract type ExternalConstraints{T} end

"""
刚体属性类
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct RigidBodyProperty{N,T,L} <: AbstractRigidBodyProperty{N,T}
	"是否可动（大概的确已经没用了）"
    movable::Bool
	"是否受约束"
    constrained::Bool
	"编号。必须在系统中唯一。"
    id::Int
	"类型。暂时没有。"
    type::Symbol
	"质量"
    mass::T
	"转动惯量"
    inertia::Union{T,SArray{Tuple{N,N},T,2,L}}
	"质心（在固连系中的）坐标"
    r̄g::SArray{Tuple{N},T,1,N}
	"连接点个数"
    nr̄ps::Int
	"连接点（在固连系中的）坐标"
    r̄ps::Vector{SArray{Tuple{N},T,1,N}}
end

"""
刚体属性构造子
$(TYPEDSIGNATURES)
"""
function RigidBodyProperty(id::Integer,movable::Bool,
                            mass::Real,inertia_input,
                            r̄g::AbstractVector,r̄ps;constrained=false)
    type = :generic
    nr̄ps = length(r̄ps)
    if !movable
        constrained = true
    end
	mtype = similar_type(inertia_input)
	vtype = SVector{size(inertia_input)[1]}
	return RigidBodyProperty(movable,constrained,id,type,mass,
							mtype(inertia_input),vtype(r̄g),nr̄ps,vtype.(r̄ps))
end

struct NaturalCoordinatesCache{ArrayT,MT,fT}
    constrained_index::Vector{Int}
    unconstrained_index::Vector{Int}
	Φi::Vector{Int}
    nΦ::Int
    funcs::fT
    M::MT
    Co::ArrayT
    Cg::ArrayT
    Cps::Vector{ArrayT}
end

function NaturalCoordinatesCache(prop::RigidBodyProperty{N,T,iT},
                                 lncs::NaturalCoordinates.LNC,
								 ci=Vector{Int}(),
								 Φi=get_all_Φi(lncs)) where {N,T,iT}
    (;mass,inertia,r̄g,nr̄ps,r̄ps) = prop
	uci = NaturalCoordinates.get_unconstrained_indices(lncs,ci)
	nΦ = length(Φi)
	cf = NaturalCoordinates.CoordinateFunctions(lncs,uci,Φi)
    M = NaturalCoordinates.make_M(cf,mass,inertia,r̄g)
    (;C,c) = cf
    Co = C(c(zeros(T,N)))
    Cg = C(c(r̄g))
    Cps = [typeof(Cg)(C(c(r̄ps[i]))) for i in 1:nr̄ps]
    if prop.movable
        if prop.constrained && ci == Vector{Int}()
            @error "Rigid body constrained, but no index specified."
        elseif !prop.constrained && !(ci == Vector{Int}())
            @error "Rigid body not constrained. No index should be specified."
        end
    end
    NaturalCoordinatesCache(ci,uci,Φi,nΦ,cf,M,Co,Cg,Cps)
end

"""
刚体状态mutable类。所有坐标在同一个惯性系中表达。
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
mutable struct RigidBodyState{N,T,L,M,cacheType} <: AbstractRigidBodyState{N,T}
	"固连系原点"
    ro::MArray{Tuple{N},T,1,N}
	"固连系旋转矩阵"
    R::MArray{Tuple{N,N},T,2,L}
	"固连系原点平移速度"
    ṙo::MArray{Tuple{N},T,1,N}
	"固连系角速度"
    ω::MArray{Tuple{M},T,1,M}
	"质心坐标"
    rg::MArray{Tuple{N},T,1,N}
	"质心速度"
    ṙg::MArray{Tuple{N},T,1,N}
	"各连接点坐标"
    rps::Vector{MArray{Tuple{N},T,1,N}} # Anchor Points in global frame
	"各连接点速度"
    ṙps::Vector{MArray{Tuple{N},T,1,N}} # Anchor Points in global frame
	"合力"
    f::MArray{Tuple{N},T,1,N}
	"合力矩"
    τ::MArray{Tuple{M},T,1,M}
	"各连接点所受作用力"
    fps::Vector{MArray{Tuple{N},T,1,N}}
	"各连接点所受力矩"
    τps::Vector{MArray{Tuple{M},T,1,M}}
	"其他重要信息"
    cache::cacheType
end

"""
刚体状态构造子
$(TYPEDSIGNATURES)
---
`ci`为约束坐标的索引。
`Φi`为约束方程的索引。
"""
function RigidBodyState(prop::RigidBodyProperty{N,T},
                        lncs::NaturalCoordinates.LNC,
                        r_input,rotation_input,ṙ_input,ω_input,
                        ci=Vector{Int}(),
						Φi=get_all_Φi(lncs)) where {N,T}
    (;r̄g,nr̄ps,r̄ps) = prop
    ro = MVector{N}(r_input)
    ṙo = MVector{N}(ṙ_input)
	if rotation_input isa Number
		rotation = rotation_matrix(rotation_input)
	else
		rotation = Matrix(rotation_input)
	end
	R = MMatrix{N,N,T,N*N}(rotation)
	if N==3
	    ω = MVector{N}(ω_input)
	else
		ω =MVector{1}([ω_input])
	end
    rg = MVector{N}(ro+R*r̄g)
    ṙg = MVector{N}(ṙo+ω×(rg-ro))
    f = @MVector zeros(T,N)
    τ = zero(ω) ###
    nr̄ps = prop.nr̄ps
    rps = [MVector{N}(ro+R*r̄ps[i]) for i in 1:nr̄ps]
    ṙps = [MVector{N}(ṙo+ω×(rps[i]-ro)) for i in 1:nr̄ps]
    fps = [@MVector zeros(T,N) for i in 1:nr̄ps]
    τps = [zero(τ) for i in 1:nr̄ps]

	cache = NaturalCoordinatesCache(prop,lncs,ci,Φi)

    RigidBodyState(ro,R,ṙo,ω,rg,ṙg,rps,ṙps,f,τ,fps,τps,cache)
end

"""
通用刚体类
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct RigidBody{N,T,L,M,cacheType,meshType} <: AbstractRigidBody{N,T}
	"刚体属性"
    prop::RigidBodyProperty{N,T,L}
	"刚体状态"
    state::RigidBodyState{N,T,L,M,cacheType}
	"可视化网格"
	mesh::meshType
end

RigidBody(prop,state) = RigidBody(prop,state,nothing)

# kinematic joint constraints

function get_rbids(rbs)
    ids = mapreduce((rb)->rb.prop.id,vcat,rbs;init=Vector{Int}())
    nb = length(ids)
	ids,nb
end

"""
刚体上任一点类
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct ID{RBType,APType}
	"指向点所在的刚体"
    rbsig::RBType
	"点在刚体上的编号"
    pid::APType
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

get_numbertype(cst::ExternalConstraints{T}) where T = T

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

function make_Φ(cst::FixedIndicesConstraint)
	@unpack indices, values = cst
	@inline @inbounds inner_Φ(q)   = q[indices]-values
	@inline @inbounds inner_Φ(q,d) = q[indices]-d
	inner_Φ
end

function make_A(cst::FixedIndicesConstraint)
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
function FixedBodyConstraint(rbs,body2q,rbid)
	rb = rbs[rbid]
	lncs = rb.state.cache.funcs.lncs
	q_rb = rb.state.coords.q
	constrained_index = find_full_constrained_index(lncs,q_rb)
	indices = body2q[rbid][constrained_index]
	values = q_rb[constrained_index]
	FixedBodyConstraint(length(indices),indices,values)
end

function make_Φ(cst::FixedBodyConstraint)
	@unpack indices, values = cst
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
刚体铰接约束类。
$(TYPEDEF)
"""
struct PinJoint{valueType,p2pType} <: ExternalConstraints{valueType}
	nconstraints::Int
	values::valueType
	p2p::p2pType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function PinJoint(p2p)
	rb1 = p2p.end1.rbsig
    nΦ = get_ndim(rb1)
	T = get_numbertype(rb1)
	values = zeros(T,nΦ)
	PinJoint(nΦ,values,p2p)
end

function make_Φ(cst::PinJoint,mem2sysfull)
	(;nconstraints,values,p2p) = cst
	(;end1,end2) = p2p
	pid1 = end1.pid
	pid2 = end2.pid
	C1 = end1.rbsig.state.cache.Cps[pid1]
	C2 = end2.rbsig.state.cache.Cps[pid2]
	function _inner_Φ(q,d)
		ret = zeros(eltype(q),nconstraints)
		q1 = @view q[mem2sysfull[end1.rbsig.prop.id]]
		q2 = @view q[mem2sysfull[end2.rbsig.prop.id]]
		ret .= C1*q1.-C2*q2
		ret
	end
	inner_Φ(q)   = _inner_Φ(q,values)
	inner_Φ(q,d) = _inner_Φ(q,d)
	inner_Φ
end

function make_A(cst::PinJoint,mem2sysfree,nq)
	(;nconstraints,values,p2p) = cst
	(;end1,end2) = p2p
	pid1 = end1.pid
	pid2 = end2.pid
	C1 = end1.rbsig.state.cache.Cps[pid1]
	C2 = end2.rbsig.state.cache.Cps[pid2]
	uci1 =  end1.rbsig.state.cache.unconstrained_index
	uci2 =  end2.rbsig.state.cache.unconstrained_index
	function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nq)
        ret[:,mem2sysfree[end1.rbsig.prop.id]] =  C1[:,uci1]
        ret[:,mem2sysfree[end2.rbsig.prop.id]] = -C2[:,uci2]
        ret
    end
end

"""
返回约束方程编号。
$(TYPEDSIGNATURES)
"""
function get_all_Φi(lncs::NaturalCoordinates.LNC)
	nΦ = NaturalCoordinates.get_nconstraints(lncs)
	collect(1:nΦ)
end
##

# operations on rigid body
"""
返回刚体平移动能。
$(TYPEDSIGNATURES)
"""
function kinetic_energy_translation(rb::AbstractRigidBody)
    (;mass) = rb.prop
	(;ṙg) = rb.state
    T = 1/2*transpose(ṙg)*mass*ṙg
end

"""
返回刚体旋转动能。
$(TYPEDSIGNATURES)
"""
function kinetic_energy_rotation(rb::AbstractRigidBody)
    (;inertia) = rb.prop
	(;ω) = rb.state
	T = 1/2*transpose(ω)*inertia*ω
end

"""
返回刚体重力势能。
$(TYPEDSIGNATURES)
"""
function potential_energy_gravity(rb::AbstractRigidBody)
	(;mass) = rb.prop
    (;rg) = rb.state
    gravity_acceleration = get_gravity(rb)
    -transpose(rg)*gravity_acceleration*mass
end
