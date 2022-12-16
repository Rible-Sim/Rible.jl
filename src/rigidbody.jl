abstract type AbstractBody{N,T} end
abstract type AbstractBodyProperty{N,T} end
abstract type AbstractBodyState{N,T} end

abstract type AbstractRigidBody{N,T} <: AbstractBody{N,T}end
abstract type AbstractRigidBodyProperty{N,T} <: AbstractBodyProperty{N,T} end
abstract type AbstractRigidBodyState{N,T} <: AbstractBodyState{N,T} end
abstract type ExternalConstraints{T} end

"""
刚体属性类
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct RigidBodyProperty{N,T} <: AbstractRigidBodyProperty{N,T}
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
    inertia::SMatrix{N,N,T}
	"质心（在固连系中的）坐标"
    r̄g::SVector{N,T}
	"连接点个数"
    nr̄ps::Int
	"连接点（在固连系中的）坐标"
    r̄ps::Vector{SVector{N,T}}
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

struct NonminimalCoordinatesCache{fT,MT,JT,VT,ArrayT}
    pres_idx::Vector{Int}
    free_idx::Vector{Int}
	Φ_mask::Vector{Int}
    nΦ::Int
    funcs::fT
    M::MT
	M⁻¹::MT
	∂Mq̇∂q::JT
	∂M⁻¹p∂q::JT
	Ṁq̇::VT
    ∂T∂qᵀ::VT
    Co::ArrayT
    Cg::ArrayT
    Cps::Vector{ArrayT}
end

function get_CoordinatesCache(prop::RigidBodyProperty{N,T},
                                 lncs::NaturalCoordinates.LNC,
								 pres_idx=Int[],
								 Φ_mask=get_Φ_mask(lncs)) where {N,T}
    (;mass,inertia,r̄g,nr̄ps,r̄ps) = prop
	free_idx = NaturalCoordinates.get_free_idx(lncs,pres_idx)
	nΦ = length(Φ_mask)
	cf = NaturalCoordinates.CoordinateFunctions(lncs,free_idx,Φ_mask)
    M = NaturalCoordinates.make_M(cf,mass,inertia,r̄g)
	M⁻¹ = inv(M)
	∂Mq̇∂q = zero(M)
	∂M⁻¹p∂q = zero(M)
	Ṁq̇ = @MVector zeros(T,size(M,2))
	∂T∂qᵀ = @MVector zeros(T,size(M,2))
    (;C,c) = cf
    Co = C(c(zeros(T,N)))
    Cg = C(c(r̄g))
    Cps = [typeof(Cg)(C(c(r̄ps[i]))) for i in 1:nr̄ps]
    if prop.movable
        if prop.constrained && pres_idx == Int[]
            @error "Rigid body constrained, but no index specified."
        elseif !prop.constrained && !(pres_idx == Int[])
            @error "Rigid body not constrained. No index should be specified."
        end
    end
    NonminimalCoordinatesCache(pres_idx,free_idx,Φ_mask,nΦ,cf,M,M⁻¹,∂Mq̇∂q,∂M⁻¹p∂q,Ṁq̇,∂T∂qᵀ,Co,Cg,Cps)
end

function get_CoordinatesCache(prop::RigidBodyProperty{N,T},
								qcs::QuaternionCoordinates.QC,
								pres_idx=Int[],
								Φ_mask=get_Φ_mask(qcs)) where {N,T}
	(;mass,inertia,r̄g,nr̄ps,r̄ps) = prop
	free_idx = deleteat!(collect(1:7),pres_idx)
	nΦ = length(Φ_mask)
	cf = QuaternionCoordinates.CoordinateFunctions(qcs)
	M = MMatrix{7,7}(Matrix(one(T)*I,7,7))
	M⁻¹ = MMatrix{7,7}(Matrix(one(T)*I,7,7))
	∂Mq̇∂q = @MMatrix zeros(T,7,7)
	∂M⁻¹p∂q = @MMatrix zeros(T,7,7)
	Ṁq̇ = @MVector zeros(T,7)
	∂T∂qᵀ = @MVector zeros(T,7)
    Co = MMatrix{3,7}(zeros(T,3,7))
	for i = 1:3
		Co[i,i] = 1
	end
    Cg = deepcopy(Co)
    Cps = [deepcopy(Co) for i in 1:nr̄ps]
    if prop.movable
        if prop.constrained && pres_idx == Int[]
            @error "Rigid body constrained, but no index specified."
        elseif !prop.constrained && !(pres_idx == Int[])
            @error "Rigid body not constrained. No index should be specified."
        end
    end
	NonminimalCoordinatesCache(pres_idx,free_idx,Φ_mask,nΦ,cf,M,M⁻¹,∂Mq̇∂q,∂M⁻¹p∂q,Ṁq̇,∂T∂qᵀ,Co,Cg,Cps)
end

"""
刚体状态mutable类。所有坐标在同一个惯性系中表达。
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
mutable struct RigidBodyState{N,T,M,cacheType} <: AbstractRigidBodyState{N,T}
	"固连系原点"
    ro::MVector{N,T}
	"固连系旋转矩阵"
    R::MMatrix{N,N,T}
	"固连系原点平移速度"
    ṙo::MVector{N,T}
	"固连系角速度"
    ω::MVector{M,T}
	"质心坐标"
    rg::MVector{N,T}
	"质心速度"
    ṙg::MVector{N,T}
	"各连接点坐标"
    rps::Vector{MVector{N,T}} # Anchor Points in global frame
	"各连接点速度"
    ṙps::Vector{MVector{N,T}} # Anchor Points in global frame
	"合力"
    f::MVector{N,T}
	"合力矩"
    τ::MVector{M,T}
	"各连接点所受作用力"
    fps::Vector{MVector{N,T}}
	"各连接点所受力矩"
    τps::Vector{MVector{M,T}}
	"其他重要信息"
    cache::cacheType
end

"""
刚体状态构造子
$(TYPEDSIGNATURES)
---
`pres_idx`为约束坐标的索引。
`Φ_mask`为约束方程的索引。
"""
function RigidBodyState(prop::RigidBodyProperty{N,T},
                        lncs,
                        r_input,rotation_input,ṙ_input,ω_input,
                        pres_idx=Int[],
						Φ_mask=get_Φ_mask(lncs)) where {N,T}
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
		ω = MVector{1}([ω_input])
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

	cache = get_CoordinatesCache(prop,lncs,pres_idx,Φ_mask)

    RigidBodyState(ro,R,ṙo,ω,rg,ṙg,rps,ṙps,f,τ,fps,τps,cache)
end

"""
通用刚体类
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct RigidBody{N,T,M,cacheType,meshType} <: AbstractRigidBody{N,T}
	"刚体属性"
    prop::RigidBodyProperty{N,T}
	"刚体状态"
    state::RigidBodyState{N,T,M,cacheType}
	"可视化网格"
	mesh::meshType
end

RigidBody(prop,state) = RigidBody(prop,state,nothing)


update_rigid!(rb::AbstractBody,q,q̇) = update_rigid!(rb.state,rb.state.cache,rb.prop,q,q̇)
move_rigid!(rb::AbstractBody,q,q̇)	= move_rigid!(rb.state,rb.state.cache,rb.prop,q,q̇)
stretch_rigid!(rb::AbstractBody,c) = stretch_rigid!(rb.state.cache,rb.prop,c)
update_transformations!(rb::AbstractBody,q) = update_transformations!(rb.state.cache,rb.state,rb.prop,q)

function update_rigid!(state::RigidBodyState,
			cache::NonminimalCoordinatesCache{<:NaturalCoordinates.CoordinateFunctions},
			prop::RigidBodyProperty,q,q̇)
	(;cache,ro,R,ṙo,ω,rg,ṙg) = state
	(;funcs,Co,Cg) = cache
	(;nmcs) = funcs
	mul!(ro, Co, q)
	mul!(ṙo, Co, q̇)
	mul!(rg, Cg, q)
	mul!(ṙg, Cg, q̇)
	R .= NaturalCoordinates.find_R(nmcs,q)
	ω .= NaturalCoordinates.find_ω(nmcs,q,q̇)
end

function update_rigid!(state::RigidBodyState,
			cache::NonminimalCoordinatesCache{<:QuaternionCoordinates.CoordinateFunctions},
			prop::RigidBodyProperty,x,ẋ)
	(;cache,ro,R,ṙo,ω,rg,ṙg) = state
	(;M,M⁻¹,∂Mq̇∂q,∂M⁻¹p∂q,∂T∂qᵀ,funcs) = cache
	# update inertia
	# @show x[4:7],ẋ[4:7]
	M .= funcs.build_M(x)
	M⁻¹ .= funcs.build_M⁻¹(x)	
	∂Mq̇∂q .= funcs.build_∂Mẋ∂x(x,ẋ)
	∂M⁻¹p∂q .= funcs.build_∂M⁻¹y∂x(x,M*ẋ)
	∂T∂qᵀ .= funcs.build_∂T∂xᵀ(x,ẋ)
	ro .= rg .= x[1:3]
	ṙo .= ṙg .= ẋ[1:3]
	R .= QuaternionCoordinates.find_R(x)
	ω .= R*QuaternionCoordinates.find_Ω(x,ẋ)
end

function stretch_rigid!(cache::NonminimalCoordinatesCache{<:NaturalCoordinates.CoordinateFunctions},
					prop::RigidBodyProperty,c)
	(;nr̄ps) = prop
	(;Cps,funcs) = cache
	nlocaldim = get_nlocaldim(cache)
	for pid in 1:nr̄ps
		Cps[pid] = funcs.C(c[nlocaldim*(pid-1)+1:nlocaldim*pid])
	end
end

function stretch_rigid!(cache::NonminimalCoordinatesCache{<:QuaternionCoordinates.CoordinateFunctions},
				prop::RigidBodyProperty,c)
end

function move_rigid!(state::RigidBodyState,
					cache::NonminimalCoordinatesCache{<:NaturalCoordinates.CoordinateFunctions},
					prop::RigidBodyProperty,q,q̇)
	(;rps,ṙps) = state
	(;Cps) = cache
	for (i,(rp,ṙp)) in enumerate(zip(rps,ṙps))
		mul!(rp, Cps[i], q)
		mul!(ṙp, Cps[i], q̇)
	end
end

function move_rigid!(state::RigidBodyState,
					cache::NonminimalCoordinatesCache{<:QuaternionCoordinates.CoordinateFunctions},
					prop::RigidBodyProperty,q,q̇)
	update_rigid!(state,cache,prop,q,q̇)
	(;r̄ps) = prop
	(;ro,R,ṙo,ω,rps,ṙps) = state
	for (i,(rp,ṙp)) in enumerate(zip(rps,ṙps))
		vp = R * r̄ps[i]
		rp .= ro .+ vp
		ṙp .= ṙo .+ ω × vp
	end
end

function update_transformations!(
		cache::NonminimalCoordinatesCache{<:NaturalCoordinates.CoordinateFunctions},
		state::RigidBodyState,
		prop::RigidBodyProperty,q)
end

function update_transformations!(
		cache::NonminimalCoordinatesCache{<:QuaternionCoordinates.CoordinateFunctions},
		state::RigidBodyState,
		prop::RigidBodyProperty,q)
	(;r̄g,r̄ps) = prop
	(;R,) = state
	(;Cg,Cps) = cache
	L = QuaternionCoordinates.Lmat(q[4:7])
	for (Cp,r̄p) in zip(Cps,r̄ps)
		# for i = 1:3 
		# 	Cp[i,i] = 1
		# end
		Cp[1:3,4:7] .= -2R*NaturalCoordinates.skew(r̄p)*L
	end 
	# for i = 1:3 
	# 	Cg[i,i] = 1
	# end
	Cg[1:3,4:7] .= -2R*NaturalCoordinates.skew(r̄g)*L
end

function generalize_force!(F,state::AbstractRigidBodyState)
	(;cache,f,fps) = state
	(;Cps,Cg,∂T∂qᵀ) = cache
	for (pid,fp) in enumerate(fps)
		# F .+= transpose(Cps[pid])*fp
		mul!(F,transpose(Cps[pid]),fp,1,1)
	end
	# F .+= transpose(Cg)*f
	mul!(F,transpose(Cg),f,1,1)
	F .+= ∂T∂qᵀ
	F
end

# kinematic joint constraints

function get_rbids(rbs)
    ids = mapreduce((rb)->rb.prop.id,vcat,rbs;init=Int[])
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
刚体铰接约束类。
$(TYPEDEF)
"""
struct PinJoint{valueType,e2eType} <: ExternalConstraints{valueType}
	nconstraints::Int
	values::valueType
	e2e::e2eType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function PinJoint(e2e)
	rb1 = e2e.end1.rbsig
    nΦ = get_ndim(rb1)
	T = get_numbertype(rb1)
	values = zeros(T,nΦ)
	PinJoint(nΦ,values,e2e)
end

function make_Φ(cst::PinJoint,indexed,numbered)
	(;nconstraints,values,e2e) = cst
	(;mem2num,num2sys) = numbered
	(;mem2sysfull) = indexed
	(;end1,end2) = e2e
	pid1 = end1.pid
	pid2 = end2.pid
	C1 = end1.rbsig.state.cache.Cps[pid1]
	C2 = end2.rbsig.state.cache.Cps[pid2]
	function inner_Φ(q)
		ret = zeros(eltype(q),nconstraints)
		q1 = @view q[mem2sysfull[end1.rbsig.prop.id]]
		q2 = @view q[mem2sysfull[end2.rbsig.prop.id]]
		ret .= C1*q1.-C2*q2
		ret
	end
	function inner_Φ(q,d,c)
		ret = zeros(eltype(q),nconstraints)
		rbid1 = end1.rbsig.prop.id
		rbid2 = end2.rbsig.prop.id
		q1 = @view q[mem2sysfull[rbid1]]
		q2 = @view q[mem2sysfull[rbid2]]
		c1 = c[num2sys[mem2num[rbid1][pid1]]]
		c2 = c[num2sys[mem2num[rbid2][pid2]]]
		ret .= end1.rbsig.state.cache.funcs.C(c1)*q1 .-
		       end2.rbsig.state.cache.funcs.C(c2)*q2
		ret
	end
	inner_Φ
end

function make_A(cst::PinJoint,indexed,numbered)
	(;nconstraints,values,e2e) = cst
	(;mem2sysfree,nfree) = indexed
	(;mem2num,num2sys) = numbered
	(;end1,end2) = e2e
	pid1 = end1.pid
	pid2 = end2.pid
	C1 = end1.rbsig.state.cache.Cps[pid1]
	C2 = end2.rbsig.state.cache.Cps[pid2]
	uci1 =  end1.rbsig.state.cache.free_idx
	uci2 =  end2.rbsig.state.cache.free_idx
	function inner_A(q)
		ret = zeros(eltype(q),nconstraints,nfree)
		ret[:,mem2sysfree[end1.rbsig.prop.id]] =  C1[:,uci1]
		ret[:,mem2sysfree[end2.rbsig.prop.id]] = -C2[:,uci2]
		ret
	end
	function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
		rbid1 = end1.rbsig.prop.id
		rbid2 = end2.rbsig.prop.id
		c1 = c[num2sys[mem2num[rbid1][pid1]]]
		c2 = c[num2sys[mem2num[rbid2][pid2]]]
        ret[:,mem2sysfree[end1.rbsig.prop.id]] =  end1.rbsig.state.cache.funcs.C(c1)[:,uci1]
        ret[:,mem2sysfree[end2.rbsig.prop.id]] = -end2.rbsig.state.cache.funcs.C(c2)[:,uci2]
        ret
    end
	inner_A
end


"""
刚体通用线性约束类。
$(TYPEDEF)
"""
struct LinearJoint{valueType} <: ExternalConstraints{valueType}
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
	LinearJoint(nΦ,values,A)
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
	(;mem2sysfree,nfull) = indexed
	(;nconstraints,values,A) = cst
	function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nfull)
        ret .= A
        ret
    end
end


"""
返回约束方程编号。
$(TYPEDSIGNATURES)
"""
function get_Φ_mask(lncs::NaturalCoordinates.LNC)
	nΦ = NaturalCoordinates.get_nconstraints(lncs)
	collect(1:nΦ)
end

get_Φ_mask(::QuaternionCoordinates.QC) = collect(1:1)
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
	(;R,ω) = rb.state
	Ω = inv(R)*ω
	T = 1/2*transpose(Ω)*inertia*Ω
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

"""
返回刚体应变势能。
$(TYPEDSIGNATURES)
"""
function potential_energy_strain(rb::AbstractRigidBody)
	zero(get_numbertype(rb))
end

