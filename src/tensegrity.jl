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
    NumberedPoints(mem2num,num2ID,num2sys,js)
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
        pid1,pid2 = abs.(row[[rbid1,rbid2]])
        is += 1
        push!(ret_raw,Point2Point(is,ID(rbs_sorted[rbid1],pid1),ID(rbs_sorted[rbid2],pid2)))
    end
    ret = TypeSortedCollection(ret_raw)
end

"""
刚体连接性类。
$(TYPEDEF)
"""
struct Connectivity{numberType,indexType,connectType,jointType,cType}
    numbered::numberType
    indexed::indexType
    connected::connectType
    jointed::jointType
    contacts::cType
end

"""
连接性构造子。
$(TYPEDSIGNATURES)
"""
function Connectivity(numbered,indexed,connected,jointed=unjoin())
	Connectivity(numbered,indexed,connected,jointed,nothing)
end

function get_nconstraints(rbs::TypeSortedCollection)
	ninconstraints = mapreduce(get_ninconstraints,+,rbs,init=0)
end

"""
刚体自然坐标状态类。
$(TYPEDEF)
"""
mutable struct NaturalCoordinatesState{T,qT,qviewT}
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
end

"""
自然坐标状态构造子。
$(TYPEDSIGNATURES)
"""
function NaturalCoordinatesState(t,q,q̇,q̈,F,λ,freei,presi)
	t = zero(eltype(q))
	q̌ = @view q[freei]
	q̌̇ = @view q̇[freei]
	q̌̈ = @view q̈[freei]
	q̃ = @view q[presi]
	q̃̇ = @view q̇[presi]
	q̃̈ = @view q̈[presi]
	F̌ = @view F[freei]
	NaturalCoordinatesState(t,q,q̇,q̈,F,λ,q̌,q̌̇,q̌̈,q̃,q̃̇,q̃̈,F̌)
end

"""
张拉整体状态类。
$(TYPEDEF)
"""
struct TensegrityState{sysT<:NaturalCoordinatesState,msT}
	system::sysT
	rigids::msT
end

"""
张拉整体状态构造子。
$(TYPEDSIGNATURES)
"""
function TensegrityState(rbs,cnt::Connectivity)
	(;indexed,jointed) = cnt
	(;nfull,ninconstraints,sysfree,syspres) = indexed
	(;mem2sysincst,mem2sysfull,mem2sysfree,mem2syspres) = indexed
	(;nexconstraints) = jointed
	nconstraints = ninconstraints + nexconstraints
	nb = length(rbs)
	ci_by_member = Vector{Vector{Int}}(undef,nb)
	uci_by_member = Vector{Vector{Int}}(undef,nb)
	foreach(rbs) do rb
		rbid = rb.prop.id
		ci_by_member[rbid] = rb.state.cache.constrained_index
		uci_by_member[rbid] = rb.state.cache.unconstrained_index
	end
	T = get_numbertype(rbs)
	t = zero(T)
	q = Vector{T}(undef,nfull)
	q̇ = zero(q)
	q̈ = zero(q)
	F = zero(q)
	λ = Vector{T}(undef,nconstraints)
	system = NaturalCoordinatesState(t,q,q̇,q̈,F,λ,sysfree,syspres)
	rigids = [
		begin
			qmem = @view q[mem2sysfull[rbid]]
			q̇mem = @view q̇[mem2sysfull[rbid]]
			q̈mem = @view q̈[mem2sysfull[rbid]]
			Fmem = @view F[mem2sysfull[rbid]]
			λmem = @view λ[mem2sysincst[rbid]]
			NaturalCoordinatesState(t,qmem,q̇mem,q̈mem,Fmem,λmem,
									uci_by_member[rbid],ci_by_member[rbid])
		end
		for rbid = 1:nb
	]
	foreach(rbs) do rb
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
struct TensegrityStructure{BodyType,StrType,TenType,CntType,StateType} <: AbstractTensegrityStructure
    ndim::Int
	ndof::Int
	nconstraints::Int
    nrigids::Int
    ncables::Int
    ntensiles::Int
    # nprespoints::Int
    rigidbodies::BodyType
    cables::StrType
    tensiles::TenType
    connectivity::CntType
	state::StateType
end

"""
张拉整体结构构造子。
$(TYPEDSIGNATURES)
"""
function TensegrityStructure(rbs,tensiles,cnt::Connectivity)
    ndim = get_ndim(rbs)
    nrigids = length(rbs)
    ncables = length(tensiles.cables)
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
    cables = tensiles.cables
	state = TensegrityState(rbs,cnt)
    tg = TensegrityStructure(
			ndim,ndof,nconstraints,
	        nrigids,ncables,ntensiles,
	        rbs,cables,tensiles,
	        cnt,state
	)
    check_jacobian_singularity(tg)
    tg
end

"""
张拉整体机器人类。
$(TYPEDEF)
"""
struct TensegrityRobot{tgT,hubT,trajT}
    tg::tgT
    hub::hubT
    traj::trajT
end

"""
清除刚体所受作用力和力矩。
$(TYPEDSIGNATURES)
"""
function clear_forces!(tg::TensegrityStructure)
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
function update_cables_apply_forces!(tg)
    (;cables) = tg
    (;connected) = tg.connectivity
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

"""
更新刚体状态。
$(TYPEDSIGNATURES)
"""
function update_rigids!(tg,q,q̇=zero(q))
	tg.state.system.q .= q
	tg.state.system.q̇ .= q̇
	update_rigids!(tg)
end
function update_rigids!(tg)
    (;rigidbodies,state) = tg
    foreach(rigidbodies) do rb
        rbid = rb.prop.id
		(;q,q̇) = state.rigids[rbid]
        (;cache,rps,ṙps,ro,ṙo,rg,ṙg) = rb.state
        (;Co,Cg,Cps) = cache
        mul!(ro, Co, q)
        mul!(ṙo, Co, q̇)
        mul!(rg, Cg, q)
        mul!(ṙg, Cg, q̇)
        for (i,(rp,ṙp)) in enumerate(zip(rps,ṙps))
            mul!(rp, Cps[i], q)
            mul!(ṙp, Cps[i], q̇)
        end
    end
end

function update_orientations!(tg)
	(;rigidbodies,state) = tg
	foreach(rigidbodies) do rb
		rbid = rb.prop.id
		(;q, q̇) = state.rigids[rbid]
		(;lncs) = rb.state.cache.funcs
		rb.state.R .= NaturalCoordinates.find_R(lncs,q)
		rb.state.ω .= NaturalCoordinates.find_ω(lncs,q,q̇)
	end
end

"""
计算系统力的大小。
$(TYPEDSIGNATURES)
"""
function generate_forces!(tg::TensegrityStructure)
	(;rigidbodies,state) = tg
	(;system,rigids) = tg.state
	system.F .= 0.0
    foreach(rigidbodies) do rb
		(;f,fps,cache) = rb.state
        (;Cps,Cg) = cache
        (;F) = rigids[rb.prop.id]
        for (pid,fp) in enumerate(fps)
            # F .+= transpose(Cps[pid])*fp
			mul!(F,transpose(Cps[pid]),fp,1,1)
        end
        # F .+= transpose(Cg)*f
		mul!(F,transpose(Cg),f,1,1)
    end
	system.F̌
end

function get_force(tg::TensegrityStructure)
	tg.state.system.F̌
end

function get_force!(F,tg::TensegrityStructure)
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

function update!(tg::TensegrityStructure;gravity=false)
    clear_forces!(tg)
    update_rigids!(tg)
    update_cables_apply_forces!(tg)
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

function build_M(tg::TensegrityStructure)
    (;nfull,mem2sysfull) = tg.connectivity.indexed
	T = get_numbertype(tg)
    M = spzeros(T,nfull,nfull)
    foreach(tg.rigidbodies) do rb
        memfull = mem2sysfull[rb.prop.id]
        M[memfull,memfull] .+= rb.state.cache.M
    end
    @assert issymmetric(M)
	M
	# symsparsecsr(M;symmetrize=true)
end

function build_M̌(tg::TensegrityStructure)
	(;sysfree) = tg.connectivity.indexed
	M = build_M(tg)
	M̌ = Symmetric(M[sysfree,sysfree])
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

make_Φ(bot::TensegrityRobot) = make_Φ(bot.tg)

function make_Φ(tg)
    (;rigidbodies,nconstraints) = tg
    (;indexed,jointed) = tg.connectivity
	(;nfree,mem2sysfull,mem2sysfree,ninconstraints,mem2sysincst) = indexed
    @inline @inbounds function inner_Φ(q)
        ret = Vector{eltype(q)}(undef,nconstraints)
        is = Ref(ninconstraints)
        #is[] += nbodydof*nfixbodies
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
            ret[is[]+1:is[]+nc,:] .= make_Φ(joint,mem2sysfull)(q)
            is[] += nc
        end
        ret
    end
    @inline @inbounds function inner_Φ(q,d)
        ret = Vector{eltype(q)}(undef,nconstraint)
        is = Ref(0)
        #is += nbodydof*nfixbodies
        for rbid in tg.mvbodyindex
            pindex = body2q[rbid]
            rb = rbs[rbid]
            # nc = rb.state.cache.nc
            # if nc > 0
            #     ret[is[]+1:is[]+nc] = rb.state.cache.cfuncs.Φ(q[pindex],d[is[]+1:is[]+nc])
            #     is[] += nc
            # end
            ret[is[]+1:is[]+nbodyc] .=rb.state.cache.funcs.Φ(q[pindex],d[is[]+1:is[]+nbodyc])
            is[] += nbodyc
        end
        foreach(csts) do cst
            nc = cst.nconstraints
            ret[is[]+1:is[]+nc] .= make_Φ(cst)(q,d[is[]+1:is[]+nc])
            is[] += nc
        end
        ret
    end
    inner_Φ
end

make_A(bot::TensegrityRobot) = make_A(bot.tg)

function make_A(tg)
    (;rigidbodies,nconstraints) = tg
    (;indexed,jointed) = tg.connectivity
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
            ret[is[]+1:is[]+nc,:] .= make_A(joint,mem2sysfree,nfree)(q)
            is[] += nc
        end
        ret
    end
end

function build_F(tg,rbid,pid,f)
    rbs = tg.rigidbodies
    Ti = build_Ti(tg,rbid)
    C = rbs[rbid].state.cache.Cps[pid]
    F = transpose(C*Ti)*f
    reshape(F,:,1)
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
    q0,q̇0 = get_q(tg)
    λ0 = get_λ(tg)
    q0,q̇0,λ0
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

get_c(bot::TensegrityRobot) = get_c(bot.tg)
function get_c(tg::TensegrityStructure)
    ndim = get_ndim(tg)
    T = get_numbertype(tg)
    (;numbered,indexed) = tg.connectivity
	(;mem2num,num2ID,num2sys,nc) = numbered
    ret = zeros(T,nc)
    foreach(tg.rigidbodies) do rb
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
    @unpack nconstraint = tg
    rbs = tg.rigidbodies
    csts = tg.constraints
    T = get_numbertype(tg)
    d = Vector{T}(undef,nconstraint)
    nbodyc = get_nbodyconstraint(tg)
    is = Ref(0)
    for rbid in tg.mvbodyindex
        rb = rbs[rbid]
        d[is[]+1:is[]+nbodyc] .= NaturalCoordinates.get_deform(rb.state.cache.funcs.lncs)
        is[] += nbodyc
    end
    foreach(csts) do cst
        nc = cst.nconstraints
        d[is[]+1:is[]+nc] .= cst.values
        is[] += nc
    end
    d
end

"""
返回系统维度。
$(TYPEDSIGNATURES)
"""
get_ndim(bot::TensegrityRobot) = get_ndim(bot.tg)
get_ndim(tg::TensegrityStructure) = get_ndim(tg.rigidbodies)
get_ndim(rbs::AbstractVector{<:AbstractRigidBody}) = get_ndim(eltype(rbs))
get_ndim(rbs::TypeSortedCollection) = get_ndim(eltype(rbs.data[1]))
get_ndim(rb::AbstractRigidBody) = get_ndim(typeof(rb))
get_ndim(::Type{<:AbstractRigidBody{N,T}}) where {N,T} = N

get_numbertype(bot::TensegrityRobot) = get_numbertype(bot.tg)
get_numbertype(tg::TensegrityStructure) = get_numbertype(tg.rigidbodies)
get_numbertype(rbs::AbstractVector{<:AbstractRigidBody}) = get_numbertype(eltype(rbs))
get_numbertype(rbs::TypeSortedCollection) = get_numbertype(eltype(rbs.data[1]))
get_numbertype(rb::AbstractRigidBody) = get_numbertype(typeof(rb))
get_numbertype(::Type{<:AbstractRigidBody{N,T}}) where {N,T} = T

"""
返回系统约束数量。
$(TYPEDSIGNATURES)
"""
get_nconstraints(tg::TensegrityStructure) = tg.nconstraints

get_ninconstraints(rb::AbstractRigidBody) = NaturalCoordinates.get_nconstraints(rb.state.cache.funcs.lncs)
get_nbodycoords(rb::AbstractRigidBody) = NaturalCoordinates.get_ncoords(rb.state.cache.funcs.lncs)
get_ndof(rb::AbstractRigidBody) = NaturalCoordinates.get_nlocaldim(rb.state.cache.funcs.lncs)
get_nlocaldim(rb::AbstractRigidBody) = NaturalCoordinates.get_nlocaldim(rb.state.cache.funcs.lncs)

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

function get_rigidbodies(tg::TensegrityStructure)
	sort_rigidbodies(tg.rigidbodies)
end

function get_rigidbars(tg::TensegrityStructure)
	rbs = get_rigidbodies(tg)
	[rb for rb in rbs
	if rb.state.cache.funcs.lncs isa Union{NaturalCoordinates.LNC2D4C,NaturalCoordinates.LNC3D6C}]
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
    [s.k for s in tg.cables]
end

"""
返回系统绳索当前长度。
$(TYPEDSIGNATURES)
"""
function get_cables_len(tg::TensegrityStructure)
    [s.state.length for s in tg.cables]
end

function get_cables_len_dot(tg::TensegrityStructure)
    [s.state.lengthdot for s in tg.cables]
end

"""
返回系统绳索变形量。
$(TYPEDSIGNATURES)
"""
function get_cables_deform(tg::TensegrityStructure)
    [s.state.length - s.state.restlen for s in tg.cables]
end

"""
返回系统绳索静止长度。
$(TYPEDSIGNATURES)
"""
function get_cables_restlen(tg::TensegrityStructure)
    [s.state.restlen for s in tg.cables]
end

"""
返回系统绳索拉力。
$(TYPEDSIGNATURES)
"""
function get_cables_tension(tg::TensegrityStructure)
    [s.state.tension for s in tg.cables]
end

"""
返回系统绳索力密度。
$(TYPEDSIGNATURES)
"""
function get_cables_force_density(tg::TensegrityStructure)
    [s.state.tension/s.state.length for s in tg.cables]
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
        for (γ,s) in zip(γs,tg.cables)]
end

function build_Y(bot)
	(;tg, hub) = bot
	(;actuators) = hub
    (;ncables,cables) = tg
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
function TensegrityRobot(tg,hub)
	update!(tg)
	# check_jacobian_singularity(tg)
	# check_stability(tg)
	traj = StructArray([deepcopy(tg.state.system)])
    TensegrityRobot(tg,hub,traj)
end

"""
重置系统状态。
$(TYPEDSIGNATURES)
"""
function reset!(bot::TensegrityRobot)
    (;tg, traj) = bot
    (;q, q̇) = traj
    clear_forces!(tg)
    update_rigids!(tg,q[begin],q̇[begin])
    update_cables_apply_forces!(tg)
    reset!(traj,1)
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
function goto_step!(bot::TensegrityRobot,that_step)
	(;tg, traj) = bot
    tg.state.system.q .= traj.q[that_step]
    tg.state.system.q̇ .= traj.q̇[that_step]
	update!(tg)
	bot
end

function analyse_slack(tg::TensegrityStructure,verbose=false)
	(;cables) = tg
	slackcases = [cable.id for cable in cables if cable.state.length <= cable.state.restlen]
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

"""
返回系统机械能。
$(TYPEDSIGNATURES)
"""
function mechanical_energy(tg::TensegrityStructure;gravity=false)
	T = kinetic_energy(tg)
	V = zero(T)
	if gravity
		V += potential_energy_gravity(tg)
	end
	if !isempty(tg.cables)
		V += sum(potential_energy.(tg.cables))
	end
	E = T+V
	@eponymtuple(T,V,E)
end

function mechanical_energy!(tg::TensegrityStructure)
	update!(tg)
	mechanical_energy(tg)
end

function mechanical_energy!(bot::TensegrityRobot;actuate=false,gravity=false)
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
