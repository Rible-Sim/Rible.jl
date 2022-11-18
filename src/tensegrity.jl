abstract type AbstractTensegrityStructure end

function check_rbid_sanity(rbs)
	ids,nb = get_rbids(rbs)
    @assert minimum(ids) == 1
    @assert maximum(ids) == nb
    @assert allunique(ids)
    ids,nb
end

struct NumberedPoints
	mem2num::Vector{Vector{Int}}
	num2ID::Vector{ID{Int,Int}}
	num2sys::Vector{Vector{Int}}
    mem2sys::Vector{Vector{Int}}
	nc::Int
end

function number(rbs)
    _,nb = check_rbid_sanity(rbs)
    nr̄ps_by_member = zeros(Int,nb)
    nld_by_member = zeros(Int,nb)
	foreach(rbs) do rb
        i = rb.prop.id
        nr̄ps_by_member[i] = rb.prop.nr̄ps
        nld_by_member[i] = get_nlocaldim(rb)
	end
    mem2num = Vector{Vector{Int}}()
    num2ID = Vector{ID{Int,Int}}()
	num2sys = Vector{Vector{Int}}()
    is = 1
	js = 0
    for rbid in 1:nb
        push!(mem2num,Vector{Int}())
		nld = nld_by_member[rbid]
        for pid in 1:nr̄ps_by_member[rbid]
            push!(num2ID,ID(rbid,pid))
            push!(mem2num[rbid],is)
            is += 1
			push!(num2sys,collect(1:nld).+js)
			js += nld
        end
    end
    mem2sys = [
        reduce(vcat,num2sys[mem2num[rbid]])
        for rbid = 1:nb
    ]
    # @show mem2sys
    NumberedPoints(mem2num,num2ID,num2sys,mem2sys,js)
end

struct IndexedMemberCoords{mem2sysType,sysType}
    nfull::Int
    nfree::Int
    npres::Int
    nmem::Int
    mem2sysfull::mem2sysType
    mem2sysfree::mem2sysType
    mem2syspres::mem2sysType
    sysfree::sysType
    syspres::sysType
	ninconstraints::Int
	mem2sysincst::mem2sysType
end

function index_inconstraints(rbs)
	ids,nmem = check_rbid_sanity(rbs)
	nincst_by_member = zeros(Int,nmem)
	foreach(rbs) do rb
		nincst_by_member[rb.prop.id] = rb.state.cache.nΦ
	end
	ninconstraints = sum(nincst_by_member)
	mem2sysincst = Vector{Vector{Int}}()
	ilast = 0
	for rbid = 1:nmem
		nincst = nincst_by_member[rbid]
		push!(mem2sysincst,collect(ilast+1:ilast+nincst))
		ilast += nincst
	end
	ninconstraints,mem2sysincst
end

function index(rbs,sharing_input=Matrix{Float64}(undef,0,0))
    ids,nmem = check_rbid_sanity(rbs)
	if size(sharing_input,2) > nmem
		@warn "Cropping the sharing matrix."
		sharing = sharing_input[:,1:nmem]
	else
		sharing = sharing_input[:,:]
	end
    sysfull = Vector{Int}()
    syspres = Vector{Int}()
    sysfree = Vector{Int}()
    mem2sysfull = Vector{Vector{Int}}()
    mem2syspres = Vector{Vector{Int}}()
    mem2sysfree = Vector{Vector{Int}}()
    ntotal_by_member = zeros(Int,nmem)
    constraineds_by_member = Vector{Vector{Int}}(undef,nmem)
    unconstraineds_by_member = Vector{Vector{Int}}(undef,nmem)
    foreach(rbs) do rb
        rbid = rb.prop.id
        ntotal_by_member[rbid] = get_nbodycoords(rb)
        constraineds_by_member[rbid] = rb.state.cache.constrained_index
        unconstraineds_by_member[rbid] = rb.state.cache.unconstrained_index
    end
    for rbid = 1:nmem
        ntotal = ntotal_by_member[rbid]
        constraineds = constraineds_by_member[rbid]
        unconstraineds = unconstraineds_by_member[rbid]
        nci = length(constraineds)
        nuci = ntotal - nci
        push!(mem2sysfull,fill(-1,ntotal))
        push!(mem2syspres,Vector{Int}())
        push!(mem2sysfree,Vector{Int}())
        unshareds = collect(1:ntotal)
        shared_indices = Vector{Int}()
        for row in eachrow(sharing)
            rbids = findall(!iszero,row)
            if rbid in rbids[begin+1:end]
                myindex = row[rbid]
                formerid = first(rbids)
                formerindex = row[formerid]
                mem2sysfull[rbid][myindex] = mem2sysfull[formerid][formerindex]
                push!(shared_indices,myindex)
            end
        end
		deleteat!(unshareds,shared_indices)
        nusi = length(unshareds)
        mem2sysfull[rbid][unshareds] = collect(length(sysfull)+1:length(sysfull)+nusi)
        append!(sysfull,mem2sysfull[rbid][unshareds])
        for usi in unshareds
            if usi in constraineds
                # pres
                push!(syspres,mem2sysfull[rbid][usi])
            else
                # free
                push!(sysfree,mem2sysfull[rbid][usi])
            end
        end
        for uci in unconstraineds
            freei = findfirst((x)->x==mem2sysfull[rbid][uci],sysfree)
            push!(mem2sysfree[rbid],freei)
        end
        for ci in constraineds
            presi = findfirst((x)->x==mem2sysfull[rbid][ci],syspres)
            push!(mem2syspres[rbid],presi)
        end
    end
	ninconstraints,mem2sysincst = index_inconstraints(rbs)
    IndexedMemberCoords(
		length(sysfull),length(sysfree),length(syspres),nmem,
		mem2sysfull,mem2sysfree,mem2syspres,
		sysfree,syspres,
		ninconstraints,mem2sysincst
	)
end

struct JointedMembers{JType}
    njoints::Int
	nexconstraints::Int
    joints::JType
end

function unjoin()
	njoints = 0
	joints = Vector{Int}()
	nexconstraints = 0
	JointedMembers(njoints,nexconstraints,joints)
end

function join(joints,indexed)
	nexconstraints = mapreduce((joint)->joint.nconstraints,+,joints,init=0)
    njoints = length(joints)
    JointedMembers(njoints,nexconstraints,joints)
end

function Base.isless(rb1::AbstractRigidBody,rb2::AbstractRigidBody)
    isless(rb1.prop.id,rb2.prop.id)
end

function sort_rigidbodies(rbs::AbstractVector{<:AbstractRigidBody})
    sort!(rbs)
end

function sort_rigidbodies(rbs::TypeSortedCollection)
    sort!(reduce(vcat,rbs.data))
end

function connect(rbs,cm_input)
    _,nb = check_rbid_sanity(rbs)
	if size(cm_input,2) > nb
		@warn "Cropping the connecting matrix."
		cm = cm_input[:,1:nb]
	else
		cm = cm_input[:,:]
	end
    rbs_sorted = sort_rigidbodies(rbs)
    ret_raw = []
    is = 0
    for row in eachrow(cm)
        rbids = findall(!iszero,row)
		if isempty(rbids)
			continue
		end
        @assert length(rbids) == 2
        @assert reduce(*,row[rbids]) < 0
        rbid1,rbid2 = ifelse(row[rbids[1]]>0,rbids,reverse(rbids))
        pid1,pid2 = Int64.(abs.(row[[rbid1,rbid2]]))
        is += 1
        push!(ret_raw,End2End(is,ID(rbs_sorted[rbid1],pid1),ID(rbs_sorted[rbid2],pid2)))
    end
    ret = TypeSortedCollection(ret_raw)
end

function cluster(rbs, cm2_input)
	rbs_sorted = sort_rigidbodies(rbs)
	ret_raw = []
	cm = cm2_input
	for row in eachrow(cm)
		iret = Vector{End2End}()
		is = 0
		rbids = findall(!iszero,row)
		if isempty(rbids)
			continue
		end
		nrbid = length(rbids)
		for i in 1:nrbid-1
			is += 1
			rbid1 = rbids[i]; rbid2 = rbids[i+1]
			pid1 = Int64(row[rbids[i]]); pid2 = Int64(row[rbids[i+1]])
			push!(iret,End2End(is,ID(rbs_sorted[rbid1],pid1),ID(rbs_sorted[rbid2],pid2)))
		end
		push!(ret_raw, iret)
	end
	ret2 = TypeSortedCollection(ret_raw)
end

function connect_and_cluster(rbs, cm_input, cm2_input)
    ret1 = connect(rbs, cm_input)
	ret2 = cluster(rbs, cm2_input)
    return (connected=ret1, clustered=ret2)
end

"""
刚体连接性类。
$(TYPEDEF)
"""
struct Connectivity{numberType,indexType,tensionType,jointType}
    numbered::numberType
    indexed::indexType
    tensioned::tensionType
    jointed::jointType
end

"""
连接性构造子。
$(TYPEDSIGNATURES)
"""
function Connectivity(numbered,indexed,tensioned)
	Connectivity(numbered,indexed,tensioned,unjoin())
end

function get_nconstraints(rbs::TypeSortedCollection)
	ninconstraints = mapreduce(get_ninconstraints,+,rbs,init=0)
end

"""
刚体自然坐标状态类。
$(TYPEDEF)
"""
mutable struct NonminimalCoordinatesState{T,qT,qviewT}
	t::T
	q::qT
	q̇::qT
	q̈::qT
    F::qT
	λ::qT
	q̌::qviewT
	q̌̇::qviewT
	q̌̈::qviewT
	q̃::qviewT
	q̃̇::qviewT
	q̃̈::qviewT
	F̌::qviewT
	c::qT
end

mutable struct ClusterNonminimalCoordinatesState{T,qT,qviewT}
	t::T
	q::qT
	q̇::qT
	q̈::qT
    F::qT
	λ::qT
    s::qT
	q̌::qviewT
	q̌̇::qviewT
	q̌̈::qviewT
	q̃::qviewT
	q̃̇::qviewT
	q̃̈::qviewT
	F̌::qviewT
end

"""
自然坐标状态构造子。
$(TYPEDSIGNATURES)
"""
function NonminimalCoordinatesState(t,q,q̇,q̈,F,λ,freei,presi,c)
	t = zero(eltype(q))
	q̌ = @view q[freei]
	q̌̇ = @view q̇[freei]
	q̌̈ = @view q̈[freei]
	q̃ = @view q[presi]
	q̃̇ = @view q̇[presi]
	q̃̈ = @view q̈[presi]
	F̌ = @view F[freei]
	NonminimalCoordinatesState(t,q,q̇,q̈,F,λ,q̌,q̌̇,q̌̈,q̃,q̃̇,q̃̈,F̌,c)
end

function ClusterNonminimalCoordinatesState(t,q,q̇,q̈,F,λ,freei,presi,s)
	t = zero(eltype(q))
	q̌ = @view q[freei]
	q̌̇ = @view q̇[freei]
	q̌̈ = @view q̈[freei]
	q̃ = @view q[presi]
	q̃̇ = @view q̇[presi]
	q̃̈ = @view q̈[presi]
	F̌ = @view F[freei]
	ClusterNonminimalCoordinatesState(t,q,q̇,q̈,F,λ,s,q̌,q̌̇,q̌̈,q̃,q̃̇,q̃̈,F̌)
end

"""
张拉整体状态类。
$(TYPEDEF)
"""
struct TensegrityState{sysT, msT}
	system::sysT
	rigids::msT
end


function rigidState2coordinates(nmcs::NaturalCoordinates.LNC,ro,R,ṙo,ω)
    NaturalCoordinates.rigidstate2naturalcoords(nmcs,ro,R,ṙo,ω)
end

function rigidState2coordinates(nmcs::QuaternionCoordinates.QC,ro,R,ṙo,ω)
    QuaternionCoordinates.rigidState2coordinates(ro,R,ṙo,ω)
end
"""
张拉整体状态构造子。
$(TYPEDSIGNATURES)
"""
function TensegrityState(rigidbodies,tensiles,cnt::Connectivity{<:Any,<:Any,<:NamedTuple{(:connected, )},<:Any})
	(;numbered,indexed,jointed) = cnt
	(;mem2num,num2sys) = numbered
	(;nfull,ninconstraints,sysfree,syspres) = indexed
	(;mem2sysincst,mem2sysfull,mem2sysfree,mem2syspres) = indexed
	(;nexconstraints) = jointed
	nconstraints = ninconstraints + nexconstraints
	nb = length(rigidbodies)
	ci_by_member = Vector{Vector{Int}}(undef,nb)
	uci_by_member = Vector{Vector{Int}}(undef,nb)
	foreach(rigidbodies) do rb
		rbid = rb.prop.id
		ci_by_member[rbid] = rb.state.cache.constrained_index
		uci_by_member[rbid] = rb.state.cache.unconstrained_index
	end
	T = get_numbertype(rigidbodies)
	t = zero(T)
	q = Vector{T}(undef,nfull)
	q̇ = zero(q)
	q̈ = zero(q)
	F = zero(q)
	λ = Vector{T}(undef,nconstraints)
	c = get_c(rigidbodies,numbered)
	system = NonminimalCoordinatesState(t,q,q̇,q̈,F,λ,sysfree,syspres,c)
	rigids = [
		begin
			qmem = @view q[mem2sysfull[rbid]]
			q̇mem = @view q̇[mem2sysfull[rbid]]
			q̈mem = @view q̈[mem2sysfull[rbid]]
			Fmem = @view F[mem2sysfull[rbid]]
			λmem = @view λ[mem2sysincst[rbid]]
			cmem = @view c[reduce(vcat,num2sys[mem2num[rbid]])]
			NonminimalCoordinatesState(t,qmem,q̇mem,q̈mem,Fmem,λmem,
									uci_by_member[rbid],ci_by_member[rbid],
									cmem)
		end
		for rbid = 1:nb
	]
	foreach(rigidbodies) do rb
		(;ro,R,ṙo,ω,cache) = rb.state
		q,q̇ = rigidState2coordinates(cache.funcs.nmcs,ro,R,ṙo,ω)
		rigids[rb.prop.id].q .= q
		rigids[rb.prop.id].q̇ .= q̇
	end
	TensegrityState(system,rigids)
end

function TensegrityState(rigidbodies,tensiles,cnt::Connectivity{<:Any,<:Any,<:NamedTuple{(:connected, :clustered)},<:Any})
	(;clustercables) = tensiles
	(;indexed,jointed) = cnt
	(;nfull,ninconstraints,sysfree,syspres) = indexed
	(;mem2sysincst,mem2sysfull,mem2sysfree,mem2syspres) = indexed
	(;nexconstraints) = jointed
	nclustercables = length(clustercables)
	ns = sum([length(clustercables[i].sps) for i in 1:nclustercables])
	nconstraints = ninconstraints + nexconstraints
	nb = length(rigidbodies)
	ci_by_member = Vector{Vector{Int}}(undef,nb)
	uci_by_member = Vector{Vector{Int}}(undef,nb)
	foreach(rigidbodies) do rb
		rbid = rb.prop.id
		ci_by_member[rbid] = rb.state.cache.constrained_index
		uci_by_member[rbid] = rb.state.cache.unconstrained_index
	end
	T = get_numbertype(rigidbodies)
	t = zero(T)
	q = Vector{T}(undef,nfull)
	q̇ = zero(q)
	q̈ = zero(q)
	F = zero(q)
    s = zeros(T, 2ns)
	λ = Vector{T}(undef,nconstraints)
	system = ClusterNonminimalCoordinatesState(t,q,q̇,q̈,F,λ,sysfree,syspres,s)
	rigids = [
		begin
			qmem = @view q[mem2sysfull[rbid]]
			q̇mem = @view q̇[mem2sysfull[rbid]]
			q̈mem = @view q̈[mem2sysfull[rbid]]
			Fmem = @view F[mem2sysfull[rbid]]
			λmem = @view λ[mem2sysincst[rbid]]
			NonminimalCoordinatesState(t,qmem,q̇mem,q̈mem,Fmem,λmem,
									uci_by_member[rbid],ci_by_member[rbid])
		end
		for rbid = 1:nb
	]
	foreach(rigidbodies) do rb
		(;ro,R,ṙo,ω,cache) = rb.state
		q,q̇ = NaturalCoordinates.rigidstate2naturalcoords(cache.funcs.lncs,ro,R,ṙo,ω)
		rigids[rb.prop.id].q .= q
		rigids[rb.prop.id].q̇ .= q̇
	end
	TensegrityState(system,rigids)
end

"""
张拉整体结构类。
$(TYPEDEF)
"""
struct TensegrityStructure{RigidType,TenType,CntType,StateType,ContactsType} <: AbstractTensegrityStructure
    ndim::Int
	ndof::Int
	nconstraints::Int
    nrigids::Int
    ntensiles::Int
    # nprespoints::Int
    rigidbodies::RigidType
    tensiles::TenType
    connectivity::CntType
	state::StateType
	contacts::ContactsType
end

"""
张拉整体结构构造子。
$(TYPEDSIGNATURES)
"""
function TensegrityStructure(rbs,tensiles,cnt::Connectivity,contacts=nothing)
    ndim = get_ndim(rbs)
    nrigids = length(rbs)
    ntensiles = sum(map(length,tensiles))
	(;nfree,ninconstraints) = cnt.indexed
	(;nexconstraints) = cnt.jointed
	nconstraints = ninconstraints + nexconstraints
	ndof = nfree - nconstraints
	if ndof <= 0
		@warn "Non positive degree of freedom: $ndof."
	end
	# nprespoints = 0
	# prespoints = nothing
	state = TensegrityState(rbs,tensiles,cnt)
    tg = TensegrityStructure(
		ndim,ndof,nconstraints,
        nrigids,ntensiles,
        rbs,tensiles,
        cnt,state,
		contacts
	)
    # check_jacobian_singularity(tg)
    tg
end

"""
张拉整体机器人类。
$(TYPEDEF)
"""
struct TensegrityRobot{tgT,hubT,trajT,contacts_trajT}
    tg::tgT
    hub::hubT
    traj::trajT
	contacts_traj::contacts_trajT
end

"""
清除刚体所受作用力和力矩。
$(TYPEDSIGNATURES)
"""
function clear_forces!(tg::AbstractTensegrityStructure)
    tg.state.system.F .= 0
	clear_forces!(tg.rigidbodies)
end
function clear_forces!(rigidbodies::AbstractVector)
    foreach(clear_forces!,rigidbodies)
end
function clear_forces!(rigidbodies::TypeSortedCollection)
    foreach(clear_forces!,rigidbodies)
end
function clear_forces!(rb::AbstractRigidBody)
	(;state) = rb
	state.f .= 0
	foreach(state.fps) do fp
		fp .= 0
	end
	state.τ .= 0
	foreach(state.τps) do τp
	  	τp .= 0
	end
end

"""
更新绳索拉力
$(TYPEDSIGNATURES)
"""
function update_tensiles!(tg::AbstractTensegrityStructure)
    update_tensiles!(tg, tg.connectivity.tensioned)
end

function update_tensiles!(tg, @eponymargs(connected,))
    (;cables) = tg.tensiles
    foreach(connected) do scnt
        scable = cables[scnt.id]
        state1 = scnt.end1.rbsig.state
        state2 = scnt.end2.rbsig.state
        pid1 = scnt.end1.pid
        pid2 = scnt.end2.pid
        p1 = state1.rps[pid1]
        ṗ1 = state1.ṙps[pid1]
        f1 = state1.fps[pid1]
        p2 = state2.rps[pid2]
        ṗ2 = state2.ṙps[pid2]
        f2 = state2.fps[pid2]
		update!(scable,p1,p2,ṗ1,ṗ2)
		f1 .+=  scable.state.force
		f2 .-=  scable.state.force
    end
end

function lengthdir(v)
    l = norm(v)
    τ = v/l
    l,τ
end

function update_tensiles!(tg, @eponymargs(clustered,))
    (;clustercables) = tg.tensiles
    id = 0
    foreach(clustered) do scnt
        id += 1
        clustercable = clustercables[id]
        (;s) = clustercable.sps
        for (segid, seg) in enumerate(clustercable.segs)
            (;state, k, c, original_restlen, prestress) = seg
            (;restlen) = state
            state1 = scnt[segid].end1.rbsig.state
            state2 = scnt[segid].end2.rbsig.state
            pid1 = scnt[segid].end1.pid
            pid2 = scnt[segid].end2.pid
            p1 = state1.rps[pid1]
            ṗ1 = state1.ṙps[pid1]
            f1 = state1.fps[pid1]
            p2 = state2.rps[pid2]
            ṗ2 = state2.ṙps[pid2]
            f2 = state2.fps[pid2]
            Δr = p2 - p1
            Δṙ = ṗ2 - ṗ1
            state.length,state.direction = lengthdir(p2-p1)
            l = state.length
            τ = state.direction
            state.lengthdot = (transpose(Δr)*Δṙ)/l
            if segid == 1
                u = restlen + s[segid]
            elseif segid == length(clustercable.segs)
                u = restlen - s[segid-1]
            else
                u = restlen + s[segid] - s[segid-1]
            end
            state.tension = k*(l-u) + prestress
            state.tension = max(state.tension, 0)
            f1 .+=  τ*state.tension
            f2 .-=  τ*state.tension
            # state.force = τ*state.tension
        end
    end
end

function update_tensiles!(tg, @eponymargs(connected, clustered))
    update_tensiles!(tg, @eponymtuple(connected))
    update_tensiles!(tg, @eponymtuple(clustered))
end

function stretch_rigids!(tg,c)
	tg.state.system.c .= c
	stretch_rigids!(tg)
end

function stretch_rigids!(tg)
	(;rigidbodies,state) = tg
    (;mem2sys) = tg.connectivity.numbered
    (;c) = state.system
    foreach(rigidbodies) do rb
        rbid = rb.prop.id
        stretch_rigid!(rb,c[mem2sys[rbid]])
    end
end

function move_rigids!(tg,q,q̇=zero(q))
	tg.state.system.q .= q
	tg.state.system.q̇ .= q̇
	move_rigids!(tg)
end

function move_rigids!(tg)
    (;rigidbodies,state) = tg
    foreach(rigidbodies) do rb
        rbid = rb.prop.id
		(;q,q̇) = state.rigids[rbid]
        move_rigid!(rb,q,q̇)
    end
end

"""
更新刚体状态。
$(TYPEDSIGNATURES)
"""
function update_rigids!(tg,q)
	tg.state.system.q .= q
    tg.state.system.q̇ .= 0.0
    update_rigids!(tg)
end

function update_rigids!(tg,q,q̇)
	tg.state.system.q .= q
    tg.state.system.q̇ .= q̇
    update_rigids!(tg)
end

function update_rigids!(tg)
	(;rigidbodies,state) = tg
	foreach(rigidbodies) do rb
		rbid = rb.prop.id
		(;q, q̇) = state.rigids[rbid]
        update_rigid!(rb,q,q̇)
        update_transformations!(rb,q)
        move_rigid!(rb,q,q̇)
	end
end

"""
计算系统力的大小。
$(TYPEDSIGNATURES)
"""
function generate_forces!(tg::TensegrityStructure)
	(;rigidbodies,state) = tg
	(;system,rigids) = state
	system.F .= 0.0
    foreach(rigidbodies) do rb
        (;F) = rigids[rb.prop.id]
		generalize_force!(F,rb.state)
    end
	system.F̌
end

function get_force(tg::AbstractTensegrityStructure)
	tg.state.system.F̌
end

function get_force!(F,tg::AbstractTensegrityStructure)
	F .= get_force(tg)
end

"""
施加重力。
$(TYPEDSIGNATURES)
"""
function apply_gravity!(tg;factor=1)
    (;rigidbodies) = tg
    gravity_acceleration = factor*get_gravity(tg)
    foreach(rigidbodies) do rb
        rb.state.f .+= gravity_acceleration*rb.prop.mass
    end
end

function update!(tg::AbstractTensegrityStructure; gravity=false)
    clear_forces!(tg)
    stretch_rigids!(tg)
    # move_rigids!(tg)
    update_rigids!(tg)
    update_tensiles!(tg)
    # update_clustercables_apply_forces!(tg)
	if gravity
		apply_gravity!(tg)
	end
	generate_forces!(tg)
end

function update!(tg::TensegrityStructure,q,q̇=zero(q))
    tg.state.system.q .= q
    tg.state.system.q̇ .= q̇
	update!(tg)
end

function build_M(tg::AbstractTensegrityStructure)
    (;nfull,mem2sysfull) = tg.connectivity.indexed
	T = get_numbertype(tg)
    M = spzeros(T,nfull,nfull)
    foreach(tg.rigidbodies) do rb
        memfull = mem2sysfull[rb.prop.id]
        M[memfull,memfull] .+= rb.state.cache.M
    end
    # @assert issymmetric(M)
	M
end

function build_M⁻¹(tg::AbstractTensegrityStructure)
    (;nfull,mem2sysfull) = tg.connectivity.indexed
	T = get_numbertype(tg)
    M⁻¹ = spzeros(T,nfull,nfull)
    foreach(tg.rigidbodies) do rb
        memfull = mem2sysfull[rb.prop.id]
        M⁻¹[memfull,memfull] .+= rb.state.cache.M⁻¹
    end
    # @assert issymmetric(M⁻¹)
	M⁻¹
end

function build_M̌(tg::AbstractTensegrityStructure)
	(;sysfree) = tg.connectivity.indexed
	M = build_M(tg)
	M̌ = Symmetric(M[sysfree,sysfree])
end

function build_∂Mq̇∂q(tg::AbstractTensegrityStructure)
    (;nfull,mem2sysfull) = tg.connectivity.indexed
	T = get_numbertype(tg)
    ∂Mq̇∂q = spzeros(T,nfull,nfull)
    foreach(tg.rigidbodies) do rb
        memfull = mem2sysfull[rb.prop.id]
        ∂Mq̇∂q[memfull,memfull] .+= rb.state.cache.∂Mq̇∂q
    end
	∂Mq̇∂q
	# symsparsecsr(M;symmetrize=true)
end

function build_∂M⁻¹p∂q(tg::AbstractTensegrityStructure)
    (;nfull,mem2sysfull) = tg.connectivity.indexed
	T = get_numbertype(tg)
    ∂M⁻¹p∂q = spzeros(T,nfull,nfull)
    foreach(tg.rigidbodies) do rb
        memfull = mem2sysfull[rb.prop.id]
        ∂M⁻¹p∂q[memfull,memfull] .+= rb.state.cache.∂M⁻¹p∂q
    end
	∂M⁻¹p∂q
	# symsparsecsr(M;symmetrize=true)
end

function make_M!(tg)
    function inner_M!(M,q)
        update_rigids!(tg,q)
        M .= build_M(tg)
    end
end

function make_M⁻¹!(tg)
    function inner_M⁻¹!(M⁻¹,q)
        update_rigids!(tg,q)
        M⁻¹ .= build_M⁻¹(tg)
    end
end

function make_Jac_M!(tg)
    function Jac_M!(∂Mq̇∂q,q,q̇)
        update_rigids!(tg,q,q̇)
        ∂Mq̇∂q .= build_∂Mq̇∂q(tg)
    end
end

function make_Jac_M⁻¹!(tg)
    function Jac_M⁻¹!(∂M⁻¹p∂q,q,q̇)
        update_rigids!(tg,q,q̇)
        ∂M⁻¹p∂q .= build_∂M⁻¹p∂q(tg)
    end
end

"""
返回系统质量矩阵。
$(TYPEDSIGNATURES)
"""
function build_MassMatrices(bot::TensegrityRobot)
	(;tg) = bot
	(;nfree,npres,sysfree,syspres) = tg.connectivity.indexed
	M = build_M(tg)
	Ḿ = M[sysfree,:]
	M̌ = Symmetric(M[sysfree,sysfree])
	M̄ =           M[sysfree,syspres]
    invM̌_raw = inv(Matrix(M̌))
    invM̌ = Symmetric(sparse(invM̌_raw))
	@eponymtuple(Ḿ,M̌,M̄,invM̌)
end

function build_∂T∂qᵀ(tg::AbstractTensegrityStructure)
    (;nfull,mem2sysfull) = tg.connectivity.indexed
    T = get_numbertype(tg)
    ∂T∂qᵀ = zeros(T,nfull)
    foreach(tg.rigidbodies) do rb
        memfull = mem2sysfull[rb.prop.id]
        ∂T∂qᵀ[memfull] .+= rb.state.cache.∂T∂qᵀ
    end
    ∂T∂qᵀ
end

make_Φ(bot::TensegrityRobot) = make_Φ(bot.tg)

function make_Φ(tg::AbstractTensegrityStructure,q0::AbstractVector)
    (;rigidbodies,nconstraints) = tg
    (;numbered,indexed,jointed) = tg.connectivity
	(;nfree,nfull,syspres,sysfree,mem2sysfull,mem2syspres,mem2sysfree,ninconstraints,mem2sysincst) = indexed

    function _inner_Φ(q̌,d,c)
		q = Vector{eltype(q̌)}(undef,nfull)
		q[syspres] .= q0[syspres]
		q[sysfree] .= q̌
        ret = Vector{eltype(q̌)}(undef,nconstraints)
		foreach(rigidbodies) do rb
			rbid = rb.prop.id
			memfull = mem2sysfull[rbid]
			memfree = mem2sysfree[rbid]
			memincst = mem2sysincst[rbid]
			if !isempty(memfree)
				ret[memincst] .= rb.state.cache.funcs.Φ(q[memfull],d[memincst])
			end
		end
        is = Ref(ninconstraints)
        foreach(jointed.joints) do joint
            nc = joint.nconstraints
            ret[is[]+1:is[]+nc] .= make_Φ(joint,indexed,numbered)(q,d[is[]+1:is[]+nc],c)
            is[] += nc
        end
        ret
    end

	function inner_Φ(q̌)
		_inner_Φ(q̌,get_d(tg),get_c(tg))
	end
	function inner_Φ(q̌,d)
		_inner_Φ(q̌,d,get_c(tg))
	end
	function inner_Φ(q̌,d,c)
		_inner_Φ(q̌,d,c)
	end

    inner_Φ
end

function make_Φ(tg::AbstractTensegrityStructure)
    (;rigidbodies,nconstraints) = tg
    (;indexed,jointed,numbered) = tg.connectivity
	(;nfree,mem2sysfull,mem2sysfree,ninconstraints,mem2sysincst) = indexed
    @inline @inbounds function inner_Φ(q)
        ret = Vector{eltype(q)}(undef,nconstraints)
        is = Ref(ninconstraints)
        foreach(rigidbodies) do rb
            rbid = rb.prop.id
			memfull = mem2sysfull[rbid]
			memfree = mem2sysfree[rbid]
			memincst = mem2sysincst[rbid]
            if !isempty(memfree)
                ret[memincst] .= rb.state.cache.funcs.Φ(q[memfull])
            end
        end
		foreach(jointed.joints) do joint
            nc = joint.nconstraints
            ret[is[]+1:is[]+nc] .= make_Φ(joint,indexed,numbered)(q)
            is[] += nc
        end
        ret
    end
    inner_Φ
end

make_A(bot::TensegrityRobot) = make_A(bot.tg)

function make_A(tg::AbstractTensegrityStructure,q0::AbstractVector)
    (;rigidbodies,nconstraints) = tg
    (;numbered,indexed,jointed) = tg.connectivity
	(;nfull,nfree,syspres,sysfree,mem2sysfull,mem2sysfree,ninconstraints,mem2sysincst) = indexed

    function _inner_A(q̌,c)
		q = Vector{eltype(q̌)}(undef,nfull)
		q[syspres] .= q0[syspres]
		q[sysfree] .= q̌
        ret = zeros(eltype(q̌),nconstraints,nfree)
        foreach(rigidbodies) do rb
            rbid = rb.prop.id
			memfull = mem2sysfull[rbid]
			memfree = mem2sysfree[rbid]
			memincst = mem2sysincst[rbid]
            if !isempty(memfree)
                ret[memincst,memfree] .= rb.state.cache.funcs.Φq(q[memfull])
            end
        end
        is = Ref(ninconstraints)
        foreach(jointed.joints) do joint
            nc = joint.nconstraints
            ret[is[]+1:is[]+nc,:] .= make_A(joint,indexed,numbered)(q,c)
            is[] += nc
        end
        ret
    end
	function inner_A(q̌)
		_inner_A(q̌,get_c(tg))
	end
	function inner_A(q̌,c)
		_inner_A(q̌,c)
	end
	inner_A
end

function make_A(tg::AbstractTensegrityStructure)
    (;rigidbodies,nconstraints) = tg
    (;numbered,indexed,jointed) = tg.connectivity
	(;nfree,mem2sysfull,mem2sysfree,ninconstraints,mem2sysincst) = indexed
    @inline @inbounds function inner_A(q)
        ret = zeros(eltype(q),nconstraints,nfree)
        is = Ref(ninconstraints)
        foreach(rigidbodies) do rb
            rbid = rb.prop.id
			memfull = mem2sysfull[rbid]
			memfree = mem2sysfree[rbid]
			memincst = mem2sysincst[rbid]
            if !isempty(memfree)
                ret[memincst,memfree] .= rb.state.cache.funcs.Φq(q[memfull])
            end
        end
        foreach(jointed.joints) do joint
            nc = joint.nconstraints
            ret[is[]+1:is[]+nc,:] .= make_A(joint,indexed,numbered)(q)
            is[] += nc
        end
        ret
    end
end

function build_F̌(tg,rbid,pid,f)
	T = get_numbertype(tg)
	(;nfree,mem2sysfree) = tg.connectivity.indexed
	F̌ = zeros(T,nfree)
	foreach(tg.rigidbodies) do rb
		if rb.prop.id == rbid
	    	C = rb.state.cache.Cps[pid]
			memfree = mem2sysfree[rbid]
	        uci = rb.state.cache.unconstrained_index
			F̌[memfree] = (transpose(C)*f)[uci,:]
	    end
	end
    reshape(F̌,:,1)
end

get_q(tg) = copy(tg.state.system.q)
get_q̇(tg) = copy(tg.state.system.q̇)
get_q̌(tg) = copy(tg.state.system.q̌)
get_q̌̇(tg) = copy(tg.state.system.q̌̇)

function get_λ(tg)
	tg.state.system.λ
end

"""
返回系统初始状态。
$(TYPEDSIGNATURES)
"""
function get_initial(tg)
    _,λ = check_static_equilibrium_output_multipliers(tg)
	q̌ = get_q̌(tg)
	q = get_q(tg)
	ℓ = get_cables_len(tg)
	s = 1 ./ℓ
	d = get_d(tg)
	c = get_c(tg)
	k = get_cables_stiffness(tg)
	μ = get_cables_restlen(tg)
	@eponymtuple(q̌,q,s,λ,d,c,k,μ,)
end

function get_polyvar(tg)
    (;nconstraints,connectivity,tensiles) = tg
    (;cables) = tensiles
    (;indexed,numbered,) = connectivity
    (;nc) = numbered
    (;ninconstraints) = indexed
    (;nfree) = indexed
    ncables = length(cables)
	# state variables
    @polyvar q̌[1:nfree]
    @polyvar s[1:ncables]
    @polyvar λ[1:nconstraints]
    # parameters
	@polyvar d[1:nconstraints]
    @polyvar c[1:nc]
    @polyvar k[1:ncables]
    @polyvar μ[1:ncables]
	@eponymtuple(q̌,s,λ,d,c,k,μ,)
end

function lucompletepiv!(A)
    n=size(A, 1)
    rowpiv=zeros(Int, n)
    colpiv=zeros(Int, n)
    for k=1:n
        Asub = abs.(A[k:n, k:n])#Search for next pivot
        _, index_max = findmax(Asub)
        μ,λ = index_max.I
        μ += k-1; λ += k-1
        rowpiv[k] = μ
        A[[k,μ], 1:n] = A[[μ,k], 1:n]
        colpiv[k] = λ
        A[1:n, [k,λ]] = A[1:n, [λ,k]]
        if A[k,k]≠0
            ρ = k+1:n
            A[ρ,k] = A[ρ,k]./A[k,k]
            A[ρ,ρ] = A[ρ,ρ] - A[ρ,k:k]*A[k:k,ρ]
        end
    end
    return (rowpiv, colpiv)
end

"""
检查雅可比矩阵奇异性
$(TYPEDSIGNATURES)
"""
function check_jacobian_singularity(tg)
	(;rigidbodies,state) = tg
    q = get_q(tg)
    A = make_A(tg)
    Aq = A(q)
    sys_rank = rank(Aq)
    if sys_rank < minimum(size(Aq))
        @warn "System's Jacobian is singular: rank(A(q))=$(sys_rank)<$(minimum(size(Aq)))"
    end
    foreach(rigidbodies) do rb
        rbid = rb.prop.id
        uci = rb.state.cache.unconstrained_index
        q_rb = state.rigids[rbid].q
        Aq_rb = rb.state.cache.funcs.Φq(q_rb)
        rb_rank = rank(Aq_rb)
        if rb_rank < minimum(size(Aq_rb))
            @warn "The $(rbid)th rigid body's Jacobian is singular: rank(A(q))=$(rb_rank)<$(minimum(size(Aq_rb)))"
        end
    end
end

get_s(bot::TensegrityRobot) = get_s(bot.tg)

function get_s(tg::TensegrityStructure)
    1 ./get_cables_len(tg)
end

function get_c(tg,rbid,pid)
    (;mem2num,num2sys) = tg.connectivity.numbered
    (;c) = tg.state.system
    cidx = mem2num[rbid][pid]
    c[num2sys[cidx]]
end

get_c(bot::TensegrityRobot) = get_c(bot.tg)
get_c(tg::TensegrityStructure) = copy(tg.state.system.c)

function get_c(rigidbodies,numbered::NumberedPoints)
    ndim = get_ndim(rigidbodies)
    T = get_numbertype(rigidbodies)
	(;mem2num,num2ID,num2sys,nc) = numbered
    ret = zeros(T,nc)
    foreach(rigidbodies) do rb
        rbid = rb.prop.id
        for i in 1:rb.prop.nr̄ps
            ip = mem2num[rbid][i]
            ret[num2sys[ip]] = rb.state.cache.funcs.c(rb.prop.r̄ps[i])
        end
    end
    ret

end
function set_C!(tg,c)
	T = get_numbertype(tg)
	(;numbered,indexed) = tg.connectivity
	(;mem2num,num2ID,num2sys,nc) = numbered
    foreach(tg.rigidbodies) do rb
        rbid = rb.prop.id
        for i in 1:rb.prop.nr̄ps
            ip = mem2num[rbid][i]
            rb.state.cache.Cps[i] = rb.state.cache.funcs.C(c[num2sys[ip]])
        end
    end
end

function get_d(tg)
    (;nconstraints,rigidbodies) = tg
	(;indexed,jointed) = tg.connectivity
	(;mem2sysincst,ninconstraints) = indexed
    T = get_numbertype(tg)
    d = Vector{T}(undef,nconstraints)
	is = Ref(ninconstraints)
    foreach(rigidbodies) do rb
        rbid = rb.prop.id
		memincst = mem2sysincst[rbid]
		if !isempty(memincst)
        	d[memincst] .= NaturalCoordinates.get_deform(rb.state.cache.funcs.lncs)
		end
    end
    foreach(jointed.joints) do joint
        nc = joint.nconstraints
        d[is[]+1:is[]+nc] .= joint.values
        is[] += nc
    end
    d
end

"""
返回系统维度。
$(TYPEDSIGNATURES)
"""
get_ndim(bot::TensegrityRobot) = get_ndim(bot.tg)
get_ndim(tg::AbstractTensegrityStructure) = get_ndim(tg.rigidbodies)
get_ndim(rbs::AbstractVector{<:AbstractRigidBody}) = get_ndim(eltype(rbs))
get_ndim(rbs::TypeSortedCollection) = get_ndim(eltype(rbs.data[1]))
get_ndim(rb::AbstractRigidBody) = get_ndim(typeof(rb))
get_ndim(::Type{<:AbstractRigidBody{N,T}}) where {N,T} = N
get_ndim(::RigidBodyProperty{N}) where {N} = N

get_numbertype(bot::TensegrityRobot) = get_numbertype(bot.tg)
get_numbertype(tg::AbstractTensegrityStructure) = get_numbertype(tg.rigidbodies)
get_numbertype(rbs::AbstractVector{<:AbstractRigidBody}) = get_numbertype(eltype(rbs))
get_numbertype(rbs::TypeSortedCollection) = get_numbertype(eltype(rbs.data[1]))
get_numbertype(rb::AbstractRigidBody) = get_numbertype(typeof(rb))
get_numbertype(::Type{<:AbstractRigidBody{N,T}}) where {N,T} = T
get_numbertype(::RigidBodyProperty{N,T}) where {N,T} = T

"""
返回系统约束数量。
$(TYPEDSIGNATURES)
"""
get_nconstraints(tg::TensegrityStructure) = tg.nconstraints

get_ninconstraints(rb::AbstractRigidBody) = get_nconstraints(rb.state.cache.funcs.nmcs)
get_nbodycoords(rb::AbstractRigidBody) = get_nbodycoords(rb.state.cache.funcs.nmcs)
get_ndof(rb::AbstractRigidBody) = get_ndof(rb.state.cache.funcs.nmcs)
get_nlocaldim(rb::AbstractRigidBody) = get_nlocaldim(rb.state.cache)
get_nlocaldim(cache::NonminimalCoordinatesCache) = get_nlocaldim(cache.funcs.nmcs)

get_nconstraints(nmcs::NaturalCoordinates.LNC) = NaturalCoordinates.get_nconstraints(nmcs)
get_nbodycoords(nmcs::NaturalCoordinates.LNC) = NaturalCoordinates.get_ncoords(nmcs)
get_ndof(nmcs::NaturalCoordinates.LNC) = NaturalCoordinates.get_ndof(nmcs)
get_nlocaldim(nmcs::NaturalCoordinates.LNC) = NaturalCoordinates.get_nlocaldim(nmcs)

get_nconstraints(nmcs::QuaternionCoordinates.QC) = QuaternionCoordinates.get_nconstraints(nmcs)
get_nbodycoords(nmcs::QuaternionCoordinates.QC) = QuaternionCoordinates.get_ncoords(nmcs)
get_ndof(nmcs::QuaternionCoordinates.QC) = QuaternionCoordinates.get_ndof(nmcs)
get_nlocaldim(nmcs::QuaternionCoordinates.QC) = QuaternionCoordinates.get_nlocaldim(nmcs)

"""
返回系统重力。
$(TYPEDSIGNATURES)
"""
get_gravity(bot::TensegrityRobot) = get_gravity(bot.tg)
get_gravity(tg::TensegrityStructure) = get_gravity(tg.rigidbodies)
get_gravity(rbs::AbstractVector{<:AbstractRigidBody}) = get_gravity(eltype(rbs))
get_gravity(rbs::TypeSortedCollection) = get_gravity(eltype(rbs.data[1]))
get_gravity(rb::AbstractRigidBody) = get_gravity(typeof(rb))
get_gravity(::Type{<:AbstractRigidBody{2,T}}) where {T} = SVector{2}(zero(T),-9.81*one(T))
get_gravity(::Type{<:AbstractRigidBody{3,T}}) where {T} = SVector{3}(zero(T),zero(T),-9.81*one(T))

get_cables_len(bot::TensegrityRobot) = get_cables_len(bot.tg)
get_cables_deform(bot::TensegrityRobot) = get_cables_deform(bot.tg)
get_cables_restlen(bot::TensegrityRobot) = get_cables_restlen(bot.tg)
get_cables_len_dot(bot::TensegrityRobot) = get_cables_len_dot(bot.tg)
get_cables_tension(bot::TensegrityRobot) = get_cables_tension(bot.tg)
get_cables_stiffness(bot::TensegrityRobot) = get_cables_stiffness(bot.tg)
get_cables_force_density(bot::TensegrityRobot) = get_cables_force_density(bot.tg)

get_rigidbodies(bot::TensegrityRobot) = get_rigidbodies(bot.tg)

function get_rigidbodies(tg::AbstractTensegrityStructure)
	sort_rigidbodies(tg.rigidbodies)
end

get_rigidbars(bot::TensegrityRobot) = get_rigidbars(bot.tg)

function get_rigidbars(tg::AbstractTensegrityStructure)
	rbs = get_rigidbodies(tg)
	[rb for rb in rbs
	if rb.state.cache.funcs.nmcs isa Union{NaturalCoordinates.LNC2D4C,NaturalCoordinates.LNC3D6C}]
end

function get_cables_len!(tg::TensegrityStructure,q)
    update_rigids!(tg,q,zero(q))
    update_cables_apply_forces!(tg)
    get_cables_len(tg)
end

"""
返回系统绳索刚度。
$(TYPEDSIGNATURES)
"""
function get_cables_stiffness(tg::TensegrityStructure)
    [s.k for s in tg.tensiles.cables]
end

"""
返回系统绳索当前长度。
$(TYPEDSIGNATURES)
"""
function get_cables_len(tg::TensegrityStructure)
    [s.state.length for s in tg.tensiles.cables]
end

function get_cables_len_dot(tg::TensegrityStructure)
    [s.state.lengthdot for s in tg.tensiles.cables]
end

"""
返回系统绳索变形量。
$(TYPEDSIGNATURES)
"""
function get_cables_deform(tg::TensegrityStructure)
    [s.state.length - s.state.restlen for s in tg.tensiles.cables]
end

"""
返回系统绳索静止长度。
$(TYPEDSIGNATURES)
"""
function get_cables_restlen(tg::TensegrityStructure)
    [s.state.restlen for s in tg.tensiles.cables]
end

"""
返回系统绳索拉力。
$(TYPEDSIGNATURES)
"""
function get_cables_tension(tg::TensegrityStructure)
    [s.state.tension for s in tg.tensiles.cables]
end

"""
返回系统绳索力密度。
$(TYPEDSIGNATURES)
"""
function get_cables_force_density(tg::TensegrityStructure)
    [s.state.tension/s.state.length for s in tg.tensiles.cables]
end

"""
返回系统绳索初始长度。
$(TYPEDSIGNATURES)
"""
function get_original_restlen(botinput::TensegrityRobot)
    bot = deepcopy(botinput)
    T = get_numbertype(bot)
    actuate!(bot,zeros(T,length(bot.hub.actuators)))
    u0 = get_cables_restlen(bot.tg)
end

function force_densities_to_restlen(tg::TensegrityStructure,γs)
    [
    begin
        l = s.state.length
        l̇ = s.state.lengthdot
        k = s.k
        c = s.c
        u = l-(γ*l-c*l̇)/k
    end
        for (γ,s) in zip(γs,tg.tensiles.cables)]
end

function build_Y(bot)
	(;tg, hub) = bot
	(;actuators) = hub
    (;cables) = tg.tensiles
	ncables = length(cables)
    nact = length(actuators)
    ret = spzeros(Int,ncables,nact)
    foreach(actuators) do act
		(;id,coupler,reg) = act
		if coupler isa Serial
	        ret[act.reg.ids,id] .= 1
		elseif coupler isa Ganged
	        is1,is2 = act.reg.ids
	        ret[is1,id] =  1
	        ret[is2,id] = -1
		else
			error("Unknown actuator type")
		end
    end
    ret
end

"""
张拉整体机器人类构造子。
$(TYPEDSIGNATURES)
"""
function TensegrityRobot(tg,hub=nothing)
	update!(tg)
	traj = StructArray([deepcopy(tg.state.system)])
	contacts_traj = [deepcopy(tg.contacts)]
    TensegrityRobot(tg,hub,traj,contacts_traj)
end

"""
重置系统状态。
$(TYPEDSIGNATURES)
"""
function reset!(bot::TensegrityRobot)
    (;tg, traj) = bot
    (;q, q̇) = traj
    clear_forces!(tg)
    update_points!(tg)
    update_rigids!(tg,q[begin],q̇[begin])
    update_tensiles!(tg)
    resize!(traj,1)
end

"""
更改初始条件。
$(TYPEDSIGNATURES)
"""
function set_new_initial!(bot::TensegrityRobot,q,q̇=zero(q))
    (;tg, traj) = bot
    traj.q[begin] .= q
    traj.q̇[begin] .= q̇
    reset!(bot)
end

"""
更新系统到指定时间步状态。
$(TYPEDSIGNATURES)
"""
function goto_step!(bot::TensegrityRobot,that_step;actuate=false)
	(;tg, traj) = bot
    tg.state.system.c .= traj.c[that_step]
    tg.state.system.q .= traj.q[that_step]
    tg.state.system.q̇ .= traj.q̇[that_step]
	if actuate
		actuate!(bot,[traj.t[that_step]])
	end
	update!(tg)
	bot
end

function analyse_slack(tg::AbstractTensegrityStructure,verbose=false)
	(;cables) = tg.tensiles
	slackcases = [cable.id for cable in cables if cable.slack && (cable.state.length <= cable.state.restlen)]
	if verbose && !isempty(slackcases)
		@show slackcases
	end
	slackcases
end

"""
返回系统动能。
$(TYPEDSIGNATURES)
"""
function kinetic_energy(tg::TensegrityStructure)
	M = build_M(tg)
	(;q̇) = tg.state.system
	T = 1/2*transpose(q̇)*M*q̇
end

"""
返回系统势能。
$(TYPEDSIGNATURES)
"""
function potential_energy_gravity(tg::TensegrityStructure)
	V = Ref(zero(get_numbertype(tg)))
	foreach(tg.rigidbodies) do rb
		V[] += potential_energy_gravity(rb)
	end
	V[]
end
function potential_energy(tg::TensegrityStructure;gravity=false)
	V = zero(get_numbertype(tg))
	if gravity
		V += potential_energy_gravity(tg)
	end
	if !isempty(tg.tensiles.cables)
		V += sum(potential_energy.(tg.tensiles.cables))
	end
	V
end
"""
返回系统机械能。
$(TYPEDSIGNATURES)
"""
function mechanical_energy(tg::TensegrityStructure;gravity=false)
	T = kinetic_energy(tg)
	V = potential_energy(tg;gravity)
	E = T+V
	@eponymtuple(T,V,E)
end

function mechanical_energy!(tg::TensegrityStructure)
	update!(tg)
	mechanical_energy(tg)
end

function mechanical_energy!(bot::TensegrityRobot;actuate=false,gravity=true)
	(;tg,traj) = bot
	StructArray([
		begin
			tg.state.system.q .= trajstate.q
	        tg.state.system.q̇ .= trajstate.q̇
			if actuate
				actuate!(bot,[trajstate.t])
			end
	        update!(tg)
			mechanical_energy(tg;gravity)
		end
		for trajstate in traj
	])
end
