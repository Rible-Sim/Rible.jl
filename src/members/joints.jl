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
	ID(rbsig,pid,0)
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

function find_order(nmcs1,nmcs2,q1,q2,X1,X2,nΦ1,nΦ2)
	T = eltype(q1)
	nq1 = length(q1)
	nq2 = length(q2)
	A = zeros(T,nΦ1+nΦ2+9,nq1+nq2)
	A[    1:nΦ1,          1:nq1    ] = NCF.make_Φq(nmcs1,collect(1:nq1),collect(1:nΦ1))(q1)
	A[nΦ1+1:nΦ1+nΦ2,  nq1+1:nq1+nq2] = NCF.make_Φq(nmcs2,collect(1:nq2),collect(1:nΦ2))(q2)
	k = nΦ1+nΦ2	
	A_rest = @view A[k+1:end,:]
	rot_jac!(A_rest,collect(1:9),
		q1,q2,
		X1,X2,
		collect(1:nq1),
		collect(1:nq2),
		collect(1:nq1),
		collect(1:nq2),
	)
	_,pidx = rref_with_pivots!(Matrix(transpose(A)),min(size(A)...)*eps(real(float(one(eltype(A))))))
	@assert length(pidx) >= 15
	order = pidx[end-2:end] .- 12
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
	rb1 = e2e.end1.rbsig
    nΦ = get_ndim(rb1)
	T = get_numbertype(rb1)
	values = zeros(T,nΦ)
	PinJoint(id,nΦ,values,e2e)
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
	(;end1,end2) = e2e
	pid1 = end1.pid
	pid2 = end2.pid
	rb1 = end1.rbsig
    nΦ = 5
	T = get_numbertype(rb1)
	nmcs1 = end1.rbsig.state.cache.funcs.nmcs
	nmcs2 = end2.rbsig.state.cache.funcs.nmcs
	state1 = end1.rbsig.state
	state2 = end2.rbsig.state
	aid1 = end1.aid
	ā1 = end1.rbsig.prop.ās[aid1]
	q1,_ = NCF.rigidstate2naturalcoords(nmcs1,state1.ro,state1.R,state1.ṙo,state1.ω)
	q2,_ = NCF.rigidstate2naturalcoords(nmcs2,state2.ro,state2.R,state2.ṙo,state2.ω)
	C1 = end1.rbsig.state.cache.Cps[pid1]
	C2 = end2.rbsig.state.cache.Cps[pid2]
	values = zeros(T,nΦ)
	values[1:3] .= C1*q1 .- C2*q2
	X1 = NCF.make_X(nmcs1,q1)
	X2 = NCF.make_X(nmcs2,q2)
	a1 = X1*ā1
	inprods = transpose(X2)*a1
	mask = 1:3 .!== argmax(abs.(inprods))
	values[4:5] = inprods[mask]
	RevoluteJoint(id,nΦ,values,e2e,mask)
end

function make_Φ(cst::RevoluteJoint,indexed,numbered)
	(;nconstraints,values,e2e,mask) = cst
	(;mem2num,num2sys) = numbered
	(;mem2sysfull) = indexed
	(;end1,end2) = e2e
	nmcs1 = end1.rbsig.state.cache.funcs.nmcs
	nmcs2 = end2.rbsig.state.cache.funcs.nmcs
	pid1 = end1.pid
	pid2 = end2.pid
	aid1 = end1.aid
	C1 = end1.rbsig.state.cache.Cps[pid1]
	C2 = end2.rbsig.state.cache.Cps[pid2]
	ā1 = end1.rbsig.prop.ās[aid1]
	function inner_Φ(q)
		ret = zeros(eltype(q),nconstraints)
		q1 = @view q[mem2sysfull[end1.rbsig.prop.id]]
		q2 = @view q[mem2sysfull[end2.rbsig.prop.id]]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		a1 = X1*ā1
		ret[1:3] .= C1*q1.-C2*q2 - values[1:3]
		ret[4:5] = transpose(X2[:,mask])*a1 - values[4:5]
		ret
	end
	function inner_Φ(q,d,c)
		ret = zeros(eltype(q),nconstraints)
		rbid1 = end1.rbsig.prop.id
		rbid2 = end2.rbsig.prop.id
		q1 = @view q[mem2sysfull[rbid1]]
		q2 = @view q[mem2sysfull[rbid2]]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		a1 = X1*ā1
		c1 = c[num2sys[mem2num[rbid1][pid1]]]
		c2 = c[num2sys[mem2num[rbid2][pid2]]]
		ret[1:3] .= end1.rbsig.state.cache.funcs.C(c1)*q1 .-
		            end2.rbsig.state.cache.funcs.C(c2)*q2
		ret[1:3] .= C1*q1.-C2*q2 .- values[1:3]
		ret[4:5] = transpose(X2[:,mask])*a1 - values[4:5]
		ret
	end
	inner_Φ
end

function make_A(cst::RevoluteJoint,indexed,numbered)
	(;nconstraints,e2e,mask) = cst
	(;mem2sysfree,mem2sysfull,nfree) = indexed
	(;mem2num,num2sys) = numbered
	(;end1,end2) = e2e
	nmcs1 = end1.rbsig.state.cache.funcs.nmcs
	nmcs2 = end2.rbsig.state.cache.funcs.nmcs
	pid1 = end1.pid
	aid1 = end1.aid
	pid2 = end2.pid
	C1 = end1.rbsig.state.cache.Cps[pid1]
	C2 = end2.rbsig.state.cache.Cps[pid2]
	uci1 =  end1.rbsig.state.cache.free_idx
	uci2 =  end2.rbsig.state.cache.free_idx
	ā1 = end1.rbsig.prop.ās[aid1]
	Y1 = NCF.get_conversion(end1.rbsig.state.cache.funcs.nmcs)
	Y2 = NCF.get_conversion(end2.rbsig.state.cache.funcs.nmcs)
	function inner_A(q)
		ret = zeros(eltype(q),nconstraints,nfree)
		rbid1 = end1.rbsig.prop.id
		rbid2 = end2.rbsig.prop.id
		q1 = @view q[mem2sysfull[rbid1]]
		q2 = @view q[mem2sysfull[rbid2]]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		a1 = X1*ā1
		ret[1:3,mem2sysfree[rbid1]] .=  C1[:,uci1]
		ret[1:3,mem2sysfree[rbid2]] .= -C2[:,uci2]
		o3 = zero(a1)
		ret[4:5,mem2sysfree[rbid1]] .= (transpose(
			kron(vcat(0,ā1),X2)[:,mask])*Y1)[:,uci1]
		ret[4:5,mem2sysfree[rbid2]] .= (transpose(
			[
				o3 o3 o3;
				a1 o3 o3; 
				o3 a1 o3;  
				o3 o3 a1; 
			][:,mask])*Y2)[:,uci2]
		ret
	end
	function inner_A(q,c)
		T = eltype(q)
        ret = zeros(eltype(q),nconstraints,nfree)
		rbid1 = end1.rbsig.prop.id
		rbid2 = end2.rbsig.prop.id
		q1 = @view q[mem2sysfull[rbid1]]
		q2 = @view q[mem2sysfull[rbid2]]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		a1 = X1*ā1
		c1 = c[num2sys[mem2num[rbid1][pid1]]]
		c2 = c[num2sys[mem2num[rbid2][pid2]]]
        ret[1:3,mem2sysfree[rbid1]] .=  end1.rbsig.state.cache.funcs.C(c1)[:,uci1]
        ret[1:3,mem2sysfree[rbid2]] .= -end2.rbsig.state.cache.funcs.C(c2)[:,uci2]
		o3 = zero(a1)
		ret[4:5,mem2sysfree[rbid1]] .= (transpose(
			kron(vcat(0,ā1),X2)[:,mask])*Y1)[:,uci1]
		ret[4:5,mem2sysfree[rbid2]] .= (transpose(
			[
				o3 o3 o3;
				a1 o3 o3;
				o3 a1 o3;  
				o3 o3 a1;
			][:,mask])*Y2)[:,uci2]
        ret
    end
	inner_A
end

"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct PrismaticJoint{valueType,e2eType,frameType,orderType} <: ExternalConstraints{valueType}
	id::Int
	nconstraints::Int
	values::valueType
	e2e::e2eType
	frame::frameType
	order::orderType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function PrismaticJoint(id,e2e)
	(;end1,end2) = e2e
    nΦ = 5
	T = get_numbertype(end1.rbsig)
	pid1 = end1.pid
	pid2 = end2.pid
	nmcs1 = end1.rbsig.state.cache.funcs.nmcs
	nmcs2 = end2.rbsig.state.cache.funcs.nmcs
	state1 = end1.rbsig.state
	state2 = end2.rbsig.state
	aid1 = end1.aid
	ā1 = end1.rbsig.prop.ās[aid1]
	frame = SpatialFrame(ā1)
	q1 = NCF.rigidstate2naturalcoords(nmcs1,state1.ro,state1.R)
	q2 = NCF.rigidstate2naturalcoords(nmcs2,state2.ro,state2.R)
	X1 = NCF.make_X(nmcs1,q1)
	X2 = NCF.make_X(nmcs2,q2)
	a1t1 = X1*frame.t1
	a1t2 = X1*frame.t2
	C1 = end1.rbsig.state.cache.Cps[pid1]
	C2 = end2.rbsig.state.cache.Cps[pid2]
	values = zeros(T,nΦ)
	values[1:2] .= transpose([a1t1 a1t2])*(C1*q1.-C2*q2)
	nΦ1 = NCF.get_nconstraints(nmcs1)
	nΦ2 = NCF.get_nconstraints(nmcs2)
	inprods = transpose(X2)*X1
	order = find_order(nmcs1,nmcs2,q1,q2,X1,X2,nΦ1,nΦ2)
	values[3:5] .= inprods[order]
	PrismaticJoint(id,nΦ,values,e2e,frame,order)
end

function make_Φ(cst::PrismaticJoint,indexed,numbered)
	(;nconstraints,values,e2e,frame,order) = cst
	(;num2sys,mem2num) = numbered
	(;mem2sysfull) = indexed
	(;end1,end2) = e2e
	pid1 = end1.pid
	pid2 = end2.pid
	nmcs1 = end1.rbsig.state.cache.funcs.nmcs
	nmcs2 = end2.rbsig.state.cache.funcs.nmcs
	rbid1 = end1.rbsig.prop.id
	rbid2 = end2.rbsig.prop.id
	function inner_Φ(q)
		ret = zeros(eltype(q),nconstraints)
		q1 = @view q[mem2sysfull[rbid1]]
		q2 = @view q[mem2sysfull[rbid2]]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		C1 = end1.rbsig.state.cache.Cps[pid1]
		C2 = end2.rbsig.state.cache.Cps[pid2]
		a1t1 = X1*frame.t1
		a1t2 = X1*frame.t2
		ret[1:2] .= transpose([a1t1 a1t2])*(C1*q1.-C2*q2) .- values[1:2]
		inprods = transpose(X2)*X1
		ret[3:5] .= inprods[order] .- values[3:5]
		ret
	end
	function inner_Φ(q,d,c)
		ret = zeros(eltype(q),nconstraints)
		q1 = @view q[mem2sysfull[rbid1]]
		q2 = @view q[mem2sysfull[rbid2]]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		a1t1 = X1*frame.t1
		a1t2 = X1*frame.t2
		c1 = c[num2sys[mem2num[rbid1]][pid1]]
		c2 = c[num2sys[mem2num[rbid2]][pid2]]
		C1 = end1.rbsig.state.cache.funcs.C(c1)*q1
		C2 = end2.rbsig.state.cache.funcs.C(c2)*q2
		ret[1:2] .= transpose([a1t1 a1t2])*(C1*q1.-C2*q2) .- values[1:2]
		inprods = transpose(X2)*X1
		ret[3:5] .= inprods[order] .- values[3:5]
		ret
	end
	inner_Φ
end

function make_A(cst::PrismaticJoint,indexed,numbered)
	(;nconstraints,e2e,frame,order) = cst
	(;mem2sysfree,mem2sysfull,nfree) = indexed
	(;num2sys,mem2num) = numbered
	(;end1,end2) = e2e
	pid1 = end1.pid
	pid2 = end2.pid
	uci1 =  end1.rbsig.state.cache.free_idx
	uci2 =  end2.rbsig.state.cache.free_idx
	rbid1 = end1.rbsig.prop.id
	rbid2 = end2.rbsig.prop.id
	memfree1 = mem2sysfree[rbid1]
	memfree2 = mem2sysfree[rbid2]
	nmcs1 = end1.rbsig.state.cache.funcs.nmcs
	nmcs2 = end2.rbsig.state.cache.funcs.nmcs
	function inner_A(q)
		ret = zeros(eltype(q),nconstraints,nfree)
		q1 = @view q[mem2sysfull[rbid1]]
		q2 = @view q[mem2sysfull[rbid2]]
		C1 = end1.rbsig.state.cache.Cps[pid1]
		C2 = end2.rbsig.state.cache.Cps[pid2]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		a1t1 = X1*frame.t1
		a1t2 = X1*frame.t2
		ret[1:2,memfree1] .= kron(
				[
					0 transpose(frame.t1); 
					0 transpose(frame.t2)
				],
				transpose(C1*q1.-C2*q2)
			)[:,uci1]
		ret[1:2,memfree1] .+= transpose([a1t1 a1t2])*C1[:,uci1]
		ret[1:2,memfree2] .-= transpose([a1t1 a1t2])*C2[:,uci2]
		k = 2
		ret_rest = @view ret[k+1:end,:]
		rot_jac!(ret_rest,order,q1,q2,X1,X2,memfree1,memfree2,uci1,uci2,)
		ret
	end
	function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
		q1 = @view q[mem2sysfull[rbid1]]
		q2 = @view q[mem2sysfull[rbid2]]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		a1t1 = X1*frame.t1
		a1t2 = X1*frame.t2
		c1 = c[num2sys[mem2num[rbid1]][pid1]]
		c2 = c[num2sys[mem2num[rbid2]][pid2]]
		C1 = end1.rbsig.state.cache.funcs.C(c1)
		C2 = end2.rbsig.state.cache.funcs.C(c2)
		ret[1:2,memfree1] .= kron(
				[
					0 transpose(frame.t1); 
					0 transpose(frame.t2)
				],
				transpose(C1*q1.-C2*q2)
			)[:,uci1]
		ret[1:2,memfree1] .+= transpose([a1t1 a1t2])*C1[:,uci1]
		ret[1:2,memfree2] .-= transpose([a1t1 a1t2])*C2[:,uci2]
		k = 2
		ret_rest = @view ret[k+1:end,:]
		rot_jac!(ret_rest,order,q1,q2,X1,X2,memfree1,memfree2,uci1,uci2,)
        ret
    end
	inner_A
end

"""
刚体铰接约束类。
$(TYPEDEF)
"""
struct CylindricalJoint{valueType,e2eType,frameType} <: ExternalConstraints{valueType}
	id::Int
	nconstraints::Int
	values::valueType
	e2e::e2eType
	frame::frameType
end

"""
铰接约束构造子。
$(TYPEDSIGNATURES)
"""
function CylindricalJoint(id,e2e)
	(;end1,end2) = e2e
    nΦ = 4
	T = get_numbertype(end1.rbsig)
	pid1 = end1.pid
	pid2 = end2.pid
	nmcs1 = end1.rbsig.state.cache.funcs.nmcs
	nmcs2 = end2.rbsig.state.cache.funcs.nmcs
	state1 = end1.rbsig.state
	state2 = end2.rbsig.state
	aid1 = end1.aid
	ā1 = end1.rbsig.prop.ās[aid1]
	frame = SpatialFrame(ā1)
	q1 = NCF.rigidstate2naturalcoords(nmcs1,state1.ro,state1.R)
	q2 = NCF.rigidstate2naturalcoords(nmcs2,state2.ro,state2.R)
	X1 = NCF.make_X(nmcs1,q1)
	X2 = NCF.make_X(nmcs2,q2)
	a1t1 = X1*frame.t1
	a1t2 = X1*frame.t2
	C1 = end1.rbsig.state.cache.Cps[pid1]
	C2 = end2.rbsig.state.cache.Cps[pid2]
	values = zeros(T,nΦ)
	values[1:2] .= transpose([a1t1 a1t2])*(C1*q1.-C2*q2)
	nΦ1 = NCF.get_nconstraints(nmcs1)
	nΦ2 = NCF.get_nconstraints(nmcs2)
	inprods = transpose(X2)*X1
	order = find_order(nmcs1,nmcs2,q1,q2,X1,X2,nΦ1,nΦ2)
	values[3:5] .= inprods[order]
	CylindricalJoint(id,nΦ,values,e2e,frame,order)
end

function make_Φ(cst::CylindricalJoint,indexed,numbered)
	(;nconstraints,values,e2e,frame,order) = cst
	(;num2sys,mem2num) = numbered
	(;mem2sysfull) = indexed
	(;end1,end2) = e2e
	pid1 = end1.pid
	pid2 = end2.pid
	nmcs1 = end1.rbsig.state.cache.funcs.nmcs
	nmcs2 = end2.rbsig.state.cache.funcs.nmcs
	rbid1 = end1.rbsig.prop.id
	rbid2 = end2.rbsig.prop.id
	function inner_Φ(q)
		ret = zeros(eltype(q),nconstraints)
		q1 = @view q[mem2sysfull[rbid1]]
		q2 = @view q[mem2sysfull[rbid2]]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		C1 = end1.rbsig.state.cache.Cps[pid1]
		C2 = end2.rbsig.state.cache.Cps[pid2]
		a1t1 = X1*frame.t1
		a1t2 = X1*frame.t2
		ret[1:2] .= transpose([a1t1 a1t2])*(C1*q1.-C2*q2) .- values[1:2]
		inprods = transpose(X2)*X1
		ret[3:5] .= inprods[order] .- values[3:5]
		ret
	end
	function inner_Φ(q,d,c)
		ret = zeros(eltype(q),nconstraints)
		q1 = @view q[mem2sysfull[rbid1]]
		q2 = @view q[mem2sysfull[rbid2]]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		a1t1 = X1*frame.t1
		a1t2 = X1*frame.t2
		c1 = c[num2sys[mem2num[rbid1]][pid1]]
		c2 = c[num2sys[mem2num[rbid2]][pid2]]
		C1 = end1.rbsig.state.cache.funcs.C(c1)*q1
		C2 = end2.rbsig.state.cache.funcs.C(c2)*q2
		ret[1:2] .= transpose([a1t1 a1t2])*(C1*q1.-C2*q2) .- values[1:2]
		inprods = transpose(X2)*X1
		ret[3:5] .= inprods[order] .- values[3:5]
		ret
	end
	inner_Φ
end

function make_A(cst::CylindricalJoint,indexed,numbered)
	(;nconstraints,e2e,frame,order) = cst
	(;mem2sysfree,mem2sysfull,nfree) = indexed
	(;num2sys,mem2num) = numbered
	(;end1,end2) = e2e
	pid1 = end1.pid
	pid2 = end2.pid
	uci1 =  end1.rbsig.state.cache.free_idx
	uci2 =  end2.rbsig.state.cache.free_idx
	rbid1 = end1.rbsig.prop.id
	rbid2 = end2.rbsig.prop.id
	memfree1 = mem2sysfree[rbid1]
	memfree2 = mem2sysfree[rbid2]
	nmcs1 = end1.rbsig.state.cache.funcs.nmcs
	nmcs2 = end2.rbsig.state.cache.funcs.nmcs
	function inner_A(q)
		ret = zeros(eltype(q),nconstraints,nfree)
		q1 = @view q[mem2sysfull[rbid1]]
		q2 = @view q[mem2sysfull[rbid2]]
		C1 = end1.rbsig.state.cache.Cps[pid1]
		C2 = end2.rbsig.state.cache.Cps[pid2]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		a1t1 = X1*frame.t1
		a1t2 = X1*frame.t2
		ret[1:2,memfree1] .= kron(
				[
					0 transpose(frame.t1); 
					0 transpose(frame.t2)
				],
				transpose(C1*q1.-C2*q2)
			)[:,uci1]
		ret[1:2,memfree1] .+= transpose([a1t1 a1t2])*C1[:,uci1]
		ret[1:2,memfree2] .-= transpose([a1t1 a1t2])*C2[:,uci2]
		k = 2
		ret_rest = @view ret[k+1:end,:]
		rot_jac!(ret_rest,order,q1,q2,X1,X2,memfree1,memfree2,uci1,uci2,)
		ret
	end
	function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
		q1 = @view q[mem2sysfull[rbid1]]
		q2 = @view q[mem2sysfull[rbid2]]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		a1t1 = X1*frame.t1
		a1t2 = X1*frame.t2
		c1 = c[num2sys[mem2num[rbid1]][pid1]]
		c2 = c[num2sys[mem2num[rbid2]][pid2]]
		C1 = end1.rbsig.state.cache.funcs.C(c1)
		C2 = end2.rbsig.state.cache.funcs.C(c2)
		ret[1:2,memfree1] .= kron(
				[
					0 transpose(frame.t1); 
					0 transpose(frame.t2)
				],
				transpose(C1*q1.-C2*q2)
			)[:,uci1]
		ret[1:2,memfree1] .+= transpose([a1t1 a1t2])*C1[:,uci1]
		ret[1:2,memfree2] .-= transpose([a1t1 a1t2])*C2[:,uci2]
		k = 2
		ret_rest = @view ret[k+1:end,:]
		rot_jac!(ret_rest,order,q1,q2,X1,X2,memfree1,memfree2,uci1,uci2,)
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
	(;end1,end2) = e2e
	pid1 = end1.pid
	pid2 = end2.pid
	rb1 = end1.rbsig
    nΦ = 6
	T = get_numbertype(rb1)
	nmcs1 = end1.rbsig.state.cache.funcs.nmcs
	nmcs2 = end2.rbsig.state.cache.funcs.nmcs
	state1 = end1.rbsig.state
	state2 = end2.rbsig.state
	q1 = NCF.rigidstate2naturalcoords(nmcs1,state1.ro,state1.R)
	q2 = NCF.rigidstate2naturalcoords(nmcs2,state2.ro,state2.R)
	X1 = NCF.make_X(nmcs1,q1)
	X2 = NCF.make_X(nmcs2,q2)
	C1 = end1.rbsig.state.cache.Cps[pid1]
	C2 = end2.rbsig.state.cache.Cps[pid2]
	values = zeros(T,nΦ)
	values[1:3] .= C1*q1.-C2*q2
	nΦ1 = NCF.get_nconstraints(nmcs1)
	nΦ2 = NCF.get_nconstraints(nmcs2)
	order = find_order(nmcs1,nmcs2,q1,q2,X1,X2,nΦ1,nΦ2)
	inprods = transpose(X2)*X1
	values[4:6] .= inprods[order]
	FixedJoint(id,nΦ,values,e2e,order)
end

function make_Φ(cst::FixedJoint,indexed,numbered)
	(;nconstraints,values,e2e,order) = cst
	(;mem2num,num2sys) = numbered
	(;mem2sysfull) = indexed
	(;end1,end2) = e2e
	nmcs1 = end1.rbsig.state.cache.funcs.nmcs
	nmcs2 = end2.rbsig.state.cache.funcs.nmcs
	pid1 = end1.pid
	pid2 = end2.pid
	C1 = end1.rbsig.state.cache.Cps[pid1]
	C2 = end2.rbsig.state.cache.Cps[pid2]
	function inner_Φ(q)
		ret = zeros(eltype(q),nconstraints)
		q1 = @view q[mem2sysfull[end1.rbsig.prop.id]]
		q2 = @view q[mem2sysfull[end2.rbsig.prop.id]]
		ret[1:3] .= C1*q1.-C2*q2 .- values[1:3]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		inprods = transpose(X2)*X1
		ret[4:6] .= inprods[order] .- values[4:6]
		ret
	end
	function inner_Φ(q,d,c)
		ret = zeros(eltype(q),nconstraints)
		rbid1 = end1.rbsig.prop.id
		rbid2 = end2.rbsig.prop.id
		q1 = @view q[mem2sysfull[rbid1]]
		q2 = @view q[mem2sysfull[rbid2]]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		c1 = c[num2sys[mem2num[rbid1][pid1]]]
		c2 = c[num2sys[mem2num[rbid2][pid2]]]
		C1 = end1.rbsig.state.cache.funcs.C(c1)*q1
		C2 = end2.rbsig.state.cache.funcs.C(c2)*q2
		ret[1:3] .= C1*q1.-C2*q2 .- values[1:3]
		inprods = transpose(X2)*X1
		ret[4:6] .= inprods[order] .- values[4:6]
		ret
	end
	inner_Φ
end

function make_A(cst::FixedJoint,indexed,numbered)
	(;nconstraints,e2e,order) = cst
	(;mem2sysfree,mem2sysfull,nfree) = indexed
	(;mem2num,num2sys) = numbered
	(;end1,end2) = e2e
	nmcs1 = end1.rbsig.state.cache.funcs.nmcs
	nmcs2 = end2.rbsig.state.cache.funcs.nmcs
	pid1 = end1.pid
	pid2 = end2.pid
	C1 = end1.rbsig.state.cache.Cps[pid1]
	C2 = end2.rbsig.state.cache.Cps[pid2]
	uci1 =  end1.rbsig.state.cache.free_idx
	uci2 =  end2.rbsig.state.cache.free_idx
	rbid1 = end1.rbsig.prop.id
	rbid2 = end2.rbsig.prop.id
	memfree1 = mem2sysfree[rbid1]
	memfree2 = mem2sysfree[rbid2]
	function inner_A(q)
		ret = zeros(eltype(q),nconstraints,nfree)
		q1 = @view q[mem2sysfull[rbid1]]
		q2 = @view q[mem2sysfull[rbid2]]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		ret[1:3,memfree1] .=  C1[:,uci1]
		ret[1:3,memfree2] .= -C2[:,uci2]
		k = 3
		ret_rest = @view ret[k+1:end,:]
		rot_jac!(ret_rest,order,q1,q2,X1,X2,memfree1,memfree2,uci1,uci2,)
		ret
	end
	function inner_A(q,c)
        ret = zeros(eltype(q),nconstraints,nfree)
		q1 = @view q[mem2sysfull[rbid1]]
		q2 = @view q[mem2sysfull[rbid2]]
		X1 = NCF.make_X(nmcs1,q1)
		X2 = NCF.make_X(nmcs2,q2)
		c1 = c[num2sys[mem2num[rbid1][pid1]]]
		c2 = c[num2sys[mem2num[rbid2][pid2]]]
        ret[1:3,memfree1] .=  end1.rbsig.state.cache.funcs.C(c1)[:,uci1]
        ret[1:3,memfree2] .= -end2.rbsig.state.cache.funcs.C(c2)[:,uci2]
		k = 3
		ret_rest = @view ret[k+1:end,:]
		rot_jac!(ret_rest,order,q1,q2,X1,X2,memfree1,memfree2,uci1,uci2,)
        ret
    end
	inner_A
end


function get_jointed_free_idx(cst)
	(;nconstraints,e2e,) = cst
	(;end1,end2) = e2e
	uci1 = end1.rbsig.state.cache.free_idx
	uci2 = end2.rbsig.state.cache.free_idx
	ncoords1 = NCF.get_ncoords(end1.rbsig.state.cache.funcs.nmcs)
	# ncoords2 = NCF.get_ncoords(end2.rbsig.state.cache.nmcs)
	uci = vcat(
		uci1,
		uci2 .+ ncoords1
	)
end

function get_jointed_free(cst,indexed)
	(;nconstraints,e2e,) = cst
	(;mem2sysfree,mem2sysfull,nfree) = indexed
	(;end1,end2) = e2e
	rbid1 = end1.rbsig.prop.id
	rbid2 = end2.rbsig.prop.id
	mem2sysfree1 = mem2sysfree[rbid1]
	mem2sysfree2 = mem2sysfree[rbid2]
	cst2sysfree = vcat(
		mem2sysfree1,
		mem2sysfree2
	)
end

function make_Φqᵀq(cst,numbered,c)
	(;e2e,) = cst
	(;mem2num,num2sys) = numbered
	(;end1,end2) = e2e
	rbid1 = end1.rbsig.prop.id
	rbid2 = end2.rbsig.prop.id
	pid1 = end1.pid
	pid2 = end2.pid
	aid1 = end1.aid
	nld1 = NCF.get_nlocaldim(end1.rbsig.state.cache.funcs.nmcs)
	nld2 = NCF.get_nlocaldim(end2.rbsig.state.cache.funcs.nmcs)
	Y1 = NCF.get_conversion(end1.rbsig.state.cache.funcs.nmcs)
	Y2 = NCF.get_conversion(end2.rbsig.state.cache.funcs.nmcs)
	cv = BlockDiagonal([Y1,Y2])
    ndim = NCF.get_ndim(end1.rbsig.state.cache.funcs.nmcs)
	T = get_numbertype(cst)
    I_Int = NCF.make_I(Int,ndim)
	c1 = c[num2sys[mem2num[rbid1][pid1]]]
	c2 = c[num2sys[mem2num[rbid2][pid2]]]
	ā1 = end1.rbsig.prop.ās[aid1]
	frame = SpatialFrame(ā1)
	t̄1 = frame.t1
	t̄2 = frame.t1
	D̄ = hcat(ā1,t̄1,t̄2)
	Φqᵀq_trans = [
		begin
			d̄ = D̄[:,id] 
			ret_raw = zeros(
				T,
				(nld1+1)+(nld2+1),
				(nld1+1)+(nld2+1),
			)
			ret_raw[1,2:nld2+1] = ret_raw[2:nld2+1,1] = -d̄
			ret_raw[2:nld2+1,2:nld2+1] = -c1*transpose(d̄)-d̄*transpose(c1)
			ret_raw[(nld1+1)+1,2:nld2+1] = ret_raw[2:nld2+1,(nld1+1)+1] = d̄
			ret_raw[(nld1+1)+2:end,2:nld2+1] = c2*transpose(d̄)
			ret_raw[2:nld2+1,(nld1+1)+2:end] = d̄*transpose(c2)
			transpose(cv)*kron(ret_raw,I_Int)*cv
		end
		for id in 1:3
	]

    Φqᵀq_rot = [
		begin
			ret_raw = zeros(
				T,
				(nld1+1)+(nld2+1),
				(nld1+1)+(nld2+1),
			)
			ret_raw[(nld1+1)+1:(nld1+1)+nld2,(nld1+1)+id] = ā1
			ret_raw[(nld1+1)+id,(nld1+1)+1:(nld1+1)+nld2] = ā1
			transpose(cv)*kron(ret_raw,I_Int)*cv
		end
		for id in 1:3
	]
	
    Φqᵀq_rot_fix = [
		begin
			ret_raw = zeros(
				T,
				(nld1+1)+(nld2+1),
				(nld1+1)+(nld2+1),
			)
			ret_raw[       1+i,(nld1+1)+j] = 1
			ret_raw[(nld1+1)+j,       1+j] = 1
			transpose(cv)*kron(ret_raw,I_Int)*cv
		end
		for i = 1:3 for j = 1:3
	]
	if cst isa FixedJoint
		Φqᵀq = vcat(zero.(Φqᵀq_trans),Φqᵀq_rot_fix[cst.order])
	elseif cst isa RevoluteJoint
		Φqᵀq = vcat(zero.(Φqᵀq_trans),Φqᵀq_rot[cst.mask])
	elseif cst isa PrismaticJoint
		Φqᵀq = vcat(Φqᵀq_trans[[2,3]],Φqᵀq_rot)
	elseif cst isa PinJoint
		Φqᵀq = vcat(zero.(Φqᵀq_trans),zero.(Φqᵀq_rot))
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

function rot_jac!(ret,order,q1,q2,X1,X2,memfree1,memfree2,uci1,uci2,)
	k = 0
	for i = 1:3
		for j = 1:3
			if 3(i-1)+j in order
				k += 1
				tpl1 = zero(q1)
				tpl2 = zero(q2) 		
				tpl1[3+3(i-1)+1:3+3i] .= X2[:,j]
				tpl2[3+3(j-1)+1:3+3j] .= X1[:,i]
				ret[k, memfree1] .= tpl1[uci1] #todo cv
				ret[k, memfree2] .= tpl2[uci2] #todo cv
			end
		end
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
	(;mem2sysfree,nfull) = indexed
	(;nconstraints,values,A) = cst
	function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nfull)
        ret .= A
        ret
    end
end

#todo 
# (joint_type == :Fixed)            && (return 0, 0) #rot order #linear trans
# (joint_type == :Prismatic)        && (return 1, 0) #rot order
# (joint_type == :Planar)           && (return 2, 0) #rot order
# (joint_type == :FixedOrientation) && (return 3, 0) #rot order
# (joint_type == :Revolute)         && (return 0, 1) #linear trans
# (joint_type == :Cylindrical)      && (return 1, 1)
# (joint_type == :PlanarAxis)       && (return 2, 1)
# (joint_type == :FreeRevolute)     && (return 3, 1)
# (joint_type == :Orbital)          && (return 0, 2) #linear trans
# (joint_type == :PrismaticOrbital) && (return 1, 2)
# (joint_type == :PlanarOrbital)    && (return 2, 2)
# (joint_type == :FreeOrbital)      && (return 3, 2)
# (joint_type == :Spherical)        && (return 0, 3) #linear trans
# (joint_type == :CylindricalFree)  && (return 1, 3)
# (joint_type == :PlanarFree)       && (return 2, 3)
# (joint_type == :Floating)         && (return 3, 3)
