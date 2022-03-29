abstract type AbstractRigidBody{N,T} end
abstract type AbstractRigidBodyProperty{N,T} end
abstract type AbstractRigidBodyState{N,T} end
abstract type ExternalConstraints{T} end

struct RigidBodyProperty{N,T,L} <: AbstractRigidBodyProperty{N,T}
    movable::Bool
    constrained::Bool
    id::Int
    type::Symbol
    mass::T
    inertia::Union{T,SArray{Tuple{N,N},T,2,L}}
    r̄g::SArray{Tuple{N},T,1,N}
    nr̄ps::Int
    r̄ps::Vector{SArray{Tuple{N},T,1,N}}
end

mutable struct RigidBodyState{N,T,L,M,cacheType} <: AbstractRigidBodyState{N,T}
    ro::MArray{Tuple{N},T,1,N}
    R::MArray{Tuple{N,N},T,2,L}
    ṙo::MArray{Tuple{N},T,1,N}
    ω::MArray{Tuple{M},T,1,M}
    rg::MArray{Tuple{N},T,1,N}
    ṙg::MArray{Tuple{N},T,1,N}
    rps::Vector{MArray{Tuple{N},T,1,N}} # Anchor Points in global frame
    ṙps::Vector{MArray{Tuple{N},T,1,N}} # Anchor Points in global frame
    f::MArray{Tuple{N},T,1,N}
    τ::MArray{Tuple{M},T,1,M}
    fps::Vector{MArray{Tuple{N},T,1,N}}
    τps::Vector{MArray{Tuple{M},T,1,M}}
    cache::cacheType
end

struct RigidBody{N,T,L,M,cacheType} <: AbstractRigidBody{N,T}
    prop::RigidBodyProperty{N,T,L}
    state::RigidBodyState{N,T,L,M,cacheType}
end

struct ConstrainedFunctions{ΦT,ΦqT,∂Aᵀλ∂qT,∂Aq̇∂qT}
    Φ::ΦT
    Φq::ΦqT
    ∂Aᵀλ∂q::∂Aᵀλ∂qT
    ∂Aq̇∂q::∂Aq̇∂qT
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

# 2D and 3D
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

function NaturalCoordinatesCache(prop::RigidBodyProperty{N,T,iT},
                                 lncs::NaturalCoordinates.LocalNaturalCoordinates,
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

function RigidBodyState(prop::RigidBodyProperty{N,T},
                        lncs::NaturalCoordinates.LocalNaturalCoordinates,
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
    ṙps = [MVector{N}(ṙo+ω×(r̄ps[i]-ro)) for i in 1:nr̄ps]
    fps = [@MVector zeros(T,N) for i in 1:nr̄ps]
    τps = [zero(τ) for i in 1:nr̄ps]

	cache = NaturalCoordinatesCache(prop,lncs,ci,Φi)

    RigidBodyState(ro,R,ṙo,ω,rg,ṙg,rps,ṙps,f,τ,fps,τps,cache)
end

# kinematic joint constraints

function get_rbids(rbs)
    ids = mapreduce((rb)->rb.prop.id,vcat,rbs;init=Vector{Int}())
    nb = length(ids)
	ids,nb
end

struct ID{RBType,APType}
    rbsig::RBType
    pid::APType
end

struct Point2Point{end1Type<:ID,end2Type<:ID}
	id::Int
	end1::end1Type
	end2::end2Type
end

struct EmptyConstraint{T} <: ExternalConstraints{T}
	nconstraints::Int64
    indices::Vector{Int64}
	values::T
end

get_numbertype(cst::ExternalConstraints{T}) where T = T

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

struct FixedIndicesConstraint{T} <: ExternalConstraints{T}
	nconstraints::Int64
    indices::Vector{Int64}
	values::T
end

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

struct FixedBodyConstraint{T} <: ExternalConstraints{T}
	nconstraints::Int64
    indices::Vector{Int64}
	values::T
end

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

struct PinJoint{valueType,p2pType} <: ExternalConstraints{valueType}
	nconstraints::Int
	values::valueType
	p2p::p2pType
end

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

function get_all_Φi(lncs::NaturalCoordinates.LocalNaturalCoordinates)
	nΦ = NaturalCoordinates.get_nconstraints(lncs)
	collect(1:nΦ)
end
##

# operations on rigid body
function kinetic_energy_translation(rb::AbstractRigidBody)
    (;mass) = rb.prop
	(;ṙg) = rb.state
    T = 1/2*transpose(ṙg)*mass*ṙg
end

function kinetic_energy_rotation(rb::AbstractRigidBody)
    (;inertia) = rb.prop
	(;ω) = rb.state
	T = 1/2*transpose(ω)*inertia*ω
end

function potential_energy_gravity(rb::AbstractRigidBody)
	(;mass) = rb.prop
    (;rg) = rb.state
    gravity_acceleration = get_gravity(rb)
    -transpose(rg)*gravity_acceleration*mass
end
