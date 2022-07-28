struct TensegrityRobot{tgT,hubT,trajT}
    tg::tgT
    hub::hubT
    traj::trajT
end
abstract type AbstractTensegrity end
abstract type TensegrityRobotTrajectory{T} end

struct ConstrainedCoordinatesTrajectory{T} <: TensegrityRobotTrajectory{T}
    ts::Vector{T}
    qs::Vector{Vector{T}}
    q̇s::Vector{Vector{T}}
    λs::Vector{Vector{T}}
end

struct SlidingConstrainedCoordinatesTrajectory{T} <: TensegrityRobotTrajectory{T}
    ts::Vector{T}
    qs::Vector{Vector{T}}
    q̇s::Vector{Vector{T}}
    λs::Vector{Vector{T}}
    s̄s::Vector{Vector{T}}
end

struct SlidingConstrainedCoordinatesTrajectoryRecordData{T,R,M} <: TensegrityRobotTrajectory{T}
    ts::Vector{T}
    qs::Vector{Vector{T}}
    q̇s::Vector{Vector{T}}
    λs::Vector{Vector{T}}
    s̄s::Vector{Vector{T}}
    iterations::Vector{R}
    OtherData::Vector{M}
end

struct ID
    rbid::Int
    apid::Int
end

# struct Connectivity{bType,sType,cType}
#     body2q::bType
#     string2ap::sType
#     contacts::cType
# end
#
# struct Cluster_Connectivity{bType,sType,cType}
#     body2q::bType
#     string2ap::sType
#     clusterstring2ap::Vector{sType}
#     contacts::cType
# end

Connectivity(b) = (body2q=b,)
Connectivity(b,s) = (body2q=b,string2ap=s)
Connectivity(b,s,c) = (body2q=b,string2ap=s,clusterstring2ap=c)

struct TensegrityStructure{BodyType,StrType,TenType,CntType,CstType} <: AbstractTensegrity
    ndim::Int
    ncoords::Int
    nconstraint::Int
    ndof::Int
    nbodies::Int
    nmvbodies::Int
    mvbodyindex::Vector{Int}
    nfixbodies::Int
    fixbodyindex::Vector{Int}
    npoints::Int
    nmvpoints::Int
    nstrings::Int
    rigidbodies::Vector{BodyType}
    strings::StrType
    tensiles::TenType
    connectivity::CntType
    constraints::CstType
end

function TensegrityStructure(rbs::Vector{rbT},tensiles::TT,cnt,
                            constraints = [EmptyConstraint()]) where {TT,rbT<:AbstractRigidBody{N,T}} where {N,T}
    ndim = N
    nbodies = length(rbs)
    mvbodyindex = [i for i in eachindex(rbs) if rbs[i].prop.movable]
    nmvbodies = length(mvbodyindex)
    fixbodyindex = [i for i in eachindex(rbs) if !rbs[i].prop.movable]
    nfixbodies = length(fixbodyindex)
    npoints = 0
    for (rbid,rb) in enumerate(rbs)
        npoints += rb.prop.naps
    end
    nmvpoints = 0
    for rbid in mvbodyindex
        nmvpoints += rbs[rbid].prop.naps
    end
    ncoords = maximum(maximum.(cnt.body2q))
    strings = tensiles.strings
    nstrings = length(strings)
    nconstraint = get_nconstraint(rbs,mvbodyindex,nmvbodies,constraints)
    ndof = ncoords - nconstraint
    tg = TensegrityStructure(
                    ndim,
                    ncoords,nconstraint,ndof,
                    nbodies,nmvbodies,mvbodyindex,nfixbodies,fixbodyindex,
                    npoints,nmvpoints,
                    nstrings,
                    rbs,strings,tensiles,
                    cnt,constraints)
    check_jacobian_singularity(tg)
    tg
end

struct ClusterTensegrityStructure{BodyType,StrType,CStrType,TenType,CntType,CstType} <: AbstractTensegrity
    ndim::Int
    ncoords::Int
    nconstraint::Int
    ndof::Int
    nbodies::Int
    nmvbodies::Int
    mvbodyindex::Vector{Int}
    nfixbodies::Int
    fixbodyindex::Vector{Int}
    npoints::Int
    nmvpoints::Int
    nstrings::Int
    nclusterstrings::Int
    nslidings::Int
    rigidbodies::Vector{BodyType}
    strings::StrType
    clusterstrings::CStrType
    tensiles::TenType
    connectivity::CntType
    constraints::CstType
end

function ClusterTensegrityStructure(rbs::Vector{rbT},tensiles::TT,cnt,
                            constraints = [EmptyConstraint()]) where {TT,rbT<:AbstractRigidBody{N,T}} where {N,T}
    ndim = N
    nbodies = length(rbs)
    mvbodyindex = [i for i in eachindex(rbs) if rbs[i].prop.movable]
    nmvbodies = length(mvbodyindex)
    fixbodyindex = [i for i in eachindex(rbs) if !rbs[i].prop.movable]
    nfixbodies = length(fixbodyindex)
    npoints = 0
    for (rbid,rb) in enumerate(rbs)
        npoints += rb.prop.naps
    end
    nmvpoints = 0
    for rbid in mvbodyindex
        nmvpoints += rbs[rbid].prop.naps
    end
    ncoords = maximum(maximum.(cnt.body2q))
    nconstraint = get_nconstraint(rbs,mvbodyindex,nmvbodies,constraints)
    ndof = ncoords - nconstraint
    strings = tensiles.strings
    nstrings = length(strings)
    clusterstrings = tensiles.clusterstrings
    nclusterstrings = length(clusterstrings)
    nslidings = sum(length(cs.sps.s) for cs in clusterstrings)
    tg = ClusterTensegrityStructure(
                    ndim,
                    ncoords,nconstraint,ndof,
                    nbodies,nmvbodies,mvbodyindex,nfixbodies,fixbodyindex,
                    npoints,nmvpoints,
                    nstrings,nclusterstrings,nslidings,
                    rbs,strings,clusterstrings,tensiles,
                    cnt,constraints)
    check_jacobian_singularity(tg)
    tg
end

function lengthdir(v)
    l = norm(v)
    τ = v/l
    l,τ
end

function reset_forces!(tg::AbstractTensegrity)
    reset_forces!.(tg.rigidbodies)
end

function update_strings!(tg)
    update_strings!(tg, tg.tensiles)
end

function update_strings!(tg, @eponymargs(strings,))
    rbs = tg.rigidbodies
    cnt = tg.connectivity
    function inner_update!(rbs, cnt, sstring::TensegrityRobots.SString)
        @unpack id,k,c= sstring
        sstate = sstring.state
        a,b = cnt.string2ap[id]
        state1 = rbs[a.rbid].state
        p1 = state1.rps[a.apid]
        ṗ1 = state1.ṙps[a.apid]
        f1 = state1.Faps[a.apid]
        state2 = rbs[b.rbid].state
        p2 = state2.rps[b.apid]
        ṗ2 = state2.ṙps[b.apid]
        f2 = state2.Faps[b.apid]
        Δr = p2 - p1
        Δṙ = ṗ2 - ṗ1
        l,τ = lengthdir(p2-p1)
        sstate.length = l
        sstate.direction = τ
        sstate.lengthdot = (transpose(Δr)*Δṙ)/l
        Δl = sstate.length - sstate.restlen
        f = k*Δl + c*sstate.lengthdot
        if Δl < 0
            sstate.tension = 0.0
        elseif f < 0
            sstate.tension = 0.0
        else
            sstate.tension = f
        end
        𝐟 = τ*sstate.tension
        f1 .+=  𝐟
        f2 .+= -𝐟
    end
    function inner_update!(rbs, cnt, sstring::TensegrityRobots.PrestressString)
        @unpack id,k,c,prestress = sstring
        sstate = sstring.state
        a,b = cnt.string2ap[id]
        state1 = rbs[a.rbid].state
        p1 = state1.rps[a.apid]
        ṗ1 = state1.ṙps[a.apid]
        f1 = state1.Faps[a.apid]
        state2 = rbs[b.rbid].state
        p2 = state2.rps[b.apid]
        ṗ2 = state2.ṙps[b.apid]
        f2 = state2.Faps[b.apid]
        Δr = p2 - p1
        Δṙ = ṗ2 - ṗ1
        l,τ = lengthdir(p2-p1)
        sstate.length = l
        sstate.direction = τ
        sstate.lengthdot = (transpose(Δr)*Δṙ)/l
        Δl = sstate.length - sstate.restlen
        f = k*Δl + c*sstate.lengthdot + prestress
        #if Δl < 0
        #    sstate.tension = 0.0
        #elseif f < 0
        #    sstate.tension = 0.0
        #else
        #    sstate.tension = f
        #end
        if f < 0
            sstate.tension = 0.0
        else
            sstate.tension = f
        end
        𝐟 = τ*sstate.tension
        f1 .+=  𝐟
        f2 .+= -𝐟
    end
    for sstring in strings
        inner_update!(rbs,cnt,sstring)
    end
end

function update_strings!(tg,@eponymargs(SMA_strings,))
    rbs = tg.rigidbodies
    cnt = tg.connectivity
    for SMA_string in SMA_strings
        @unpack id,law = SMA_string
        sstate = SMA_string.state
        a,b = cnt.string2ap[id]
        state1 = rbs[a.rbid].state
        p1 = state1.rps[a.apid]
        ṗ1 = state1.ṙps[a.apid]
        f1 = state1.Faps[a.apid]
        state2 = rbs[b.rbid].state
        p2 = state2.rps[b.apid]
        ṗ2 = state2.ṙps[b.apid]
        f2 = state2.Faps[b.apid]
        Δr = p2 - p1
        Δṙ = ṗ2 - ṗ1
        l,τ = lengthdir(p2-p1)
        sstate.length = l
        sstate.direction = τ
        sstate.lengthdot = (transpose(Δr)*Δṙ)/l
        Δl = sstate.length - sstate.restlen
        f = law(Δl)
        if Δl < 0
            sstate.tension = 0.0
        elseif f < 0
            sstate.tension = 0.0
        else
            sstate.tension = f
        end
        𝐟 = τ*sstate.tension
        f1 .+=  𝐟
        f2 .+= -𝐟
    end
end

<<<<<<< Updated upstream
function update_strings!(tg, @eponymargs(clusterstrings,))
    rbs = tg.rigidbodies
    cnt = tg.connectivity
    for clusterstring in clusterstrings
        s = clusterstring.sps.s
        for (segid, seg) in enumerate(clusterstring.segs)
            @unpack k,c,prestress,original_restlen = seg
            @unpack restlen = seg.state
            u0 = restlen
            segstate = seg.state
            a,b = cnt.clusterstring2ap[clusterstring.ID][segid]
            state1 = rbs[a.rbid].state
            p1 = state1.rps[a.apid]
            ṗ1 = state1.ṙps[a.apid]
            f1 = state1.Faps[a.apid]
            state2 = rbs[b.rbid].state
            p2 = state2.rps[b.apid]
            ṗ2 = state2.ṙps[b.apid]
            f2 = state2.Faps[b.apid]
            Δr = p2 - p1
            Δṙ = ṗ2 - ṗ1
            segstate.length,segstate.direction = lengthdir(p2-p1)
            l = segstate.length
            τ = segstate.direction
            segstate.lengthdot = (transpose(Δr)*Δṙ)/l
            if segid == 1
                u = u0 + s[segid]
            elseif segid == length(clusterstring.segs)
                u = u0 - s[segid-1]
            else
                u = u0 + s[segid] - s[segid-1]
            end
            #u = u0
            #segstate.tension = k*abs(u0)/u*(l-u)
            segstate.tension = k*(l-u) + prestress
            #if u0 < 0
            #    wait()
            #end

            if segstate.tension < 0
                #@show 1
                segstate.tension = 0
            end
            𝐟 = τ*segstate.tension
            f1 .+=  𝐟
            f2 .+= -𝐟
        end
    end
=======
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
        push!(ret_raw,Point2Point(is,ID(rbs_sorted[rbid1],pid1),ID(rbs_sorted[rbid2],pid2)))
    end
    ret = TypeSortedCollection(ret_raw)
end

function connect(rbs, cm_input, cm2_input)
    ret1 = connect(rbs, cm_input)
    rbs_sorted = sort_rigidbodies(rbs)
    ret_raw = []
    cm = cm2_input
    for row in eachrow(cm)
        iret = Vector{Point2Point}()
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
            push!(iret,Point2Point(is,ID(rbs_sorted[rbid1],pid1),ID(rbs_sorted[rbid2],pid2)))
        end
        push!(ret_raw, iret)
    end
    ret2 = TypeSortedCollection(ret_raw)
    return (cables=ret1, clustercables=ret2)
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

mutable struct ClusterNaturalCoordinatesState{T,qT,sT,qviewT}
	t::T
	q::qT
	q̇::qT
	q̈::qT
    F::qT
	λ::qT
    s::sT
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

function ClusterNaturalCoordinatesState(t,q,q̇,q̈,F,λ,freei,presi,s)
	t = zero(eltype(q))
	q̌ = @view q[freei]
	q̌̇ = @view q̇[freei]
	q̌̈ = @view q̈[freei]
	q̃ = @view q[presi]
	q̃̇ = @view q̇[presi]
	q̃̈ = @view q̈[presi]
	F̌ = @view F[freei]
	ClusterNaturalCoordinatesState(t,q,q̇,q̈,F,λ,s,q̌,q̌̇,q̌̈,q̃,q̃̇,q̃̈,F̌)
end

"""
张拉整体状态类。
$(TYPEDEF)
"""
struct TensegrityState{sysT, msT}
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

function TensegrityState(rbs,cnt::Connectivity, ns)
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
    s = zeros(Float64, 2ns)
	λ = Vector{T}(undef,nconstraints)
	system = ClusterNaturalCoordinatesState(t,q,q̇,q̈,F,λ,sysfree,syspres,s)
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
    # check_jacobian_singularity(tg)
    tg
end

struct ClusterTensegrityStructure{BodyType,StrType,CStrType,TenType,CntType,StateType} <: AbstractTensegrityStructure
    ndim::Int
	ndof::Int
	nconstraints::Int
    nrigids::Int
    ncables::Int
    nclustercables::Int
    ntensiles::Int
    # nprespoints::Int
    rigidbodies::BodyType
    cables::StrType
    clustercables::CStrType
    tensiles::TenType
    connectivity::CntType
	state::StateType
end

function ClusterTensegrityStructure(rbs,tensiles,cnt::Connectivity)
    ndim = get_ndim(rbs)
    nrigids = length(rbs)
    ncables = length(tensiles.cables)
    nclustercables = length(tensiles.clustercables)
    ns = sum([length(tensiles.clustercables[i].sps) for i in 1:nclustercables])
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
    (;cables, clustercables) = tensiles
    # cables = tensiles.cables
	state = TensegrityState(rbs,cnt,ns)
    tg = ClusterTensegrityStructure(
			ndim,ndof,nconstraints,
	        nrigids,ncables,nclustercables,ntensiles,
	        rbs,cables,clustercables,tensiles,
	        cnt,state
	)
    # check_jacobian_singularity(tg)
    tg
>>>>>>> Stashed changes
end

function update_strings!(tg, @eponymargs(strings,clusterstrings))
    update_strings!(tg,@eponymtuple(strings))
    update_strings!(tg,@eponymtuple(clusterstrings))
end

<<<<<<< Updated upstream
function apply_actuation(tg, css_id::Union{Int64,Vector{Int64}}, cs_id::Union{Int64,Vector{Int64}}, apply_fun)
    @unpack clusterstrings = tg
    for css in css_id
        for cs in cs_id
            clusterstrings[css].segs[cs].state.restlen += apply_fun
        end
    end
end


distribute_q_to_rbs!(tg,globalq) = distribute_q_to_rbs!(tg,globalq,zero(globalq))
function distribute_q_to_rbs!(tg,globalq,globalq̇)
    rbs = tg.rigidbodies
    cnt = tg.connectivity
    for rbid in tg.mvbodyindex
        pindex = cnt.body2q[rbid]
        @unpack q, q̇ = rbs[rbid].state.coords
        q .= globalq[pindex]
        q̇ .= globalq̇[pindex]
        @unpack cache,rps,ṙps,ro,ṙo,rg,ṙg = rbs[rbid].state
        @unpack Co,Cg,Cp = cache
=======
"""
清除刚体所受作用力和力矩。
$(TYPEDSIGNATURES)
"""
function clear_forces!(tg::TensegrityStructure)
    tg.state.system.F .= 0
	clear_forces!(tg.rigidbodies)
end
function clear_forces!(tg::ClusterTensegrityStructure)
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
function update_cables!(tg)
    update_cables!(tg, tg.tensiles)
end

function update_cables!(tg, @eponymargs(cables,))
    (;cables) = tg
    (;connected) = tg.connectivity
    foreach(connected.cables) do scnt
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

function update_cables!(tg, @eponymargs(clustercables,))
    (;clustercables) = tg
    (;connected) = tg.connectivity
    id = 0
    foreach(connected.clustercables) do scnt
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

function update_cables!(tg, @eponymargs(cables, clustercables))
    update_cables!(tg, @eponymtuple(cables))
    update_cables!(tg, @eponymtuple(clustercables))
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
    (;mem2sysfull) = tg.connectivity.indexed
    globalq = tg.state.system.q
    foreach(rigidbodies) do rb
        rbid = rb.prop.id
		(;q,q̇) = state.rigids[rbid]
        (;cache,rps,ṙps,ro,ṙo,rg,ṙg) = rb.state
        (;Co,Cg,Cps) = cache
        pindex = mem2sysfull[rbid]
        q .= globalq[pindex]
>>>>>>> Stashed changes
        mul!(ro, Co, q)
        mul!(ṙo, Co, q̇)
        mul!(rg, Cg, q)
        mul!(ṙg, Cg, q̇)
        for (i,(rp,ṙp)) in enumerate(zip(rps,ṙps))
            mul!(rp, Cp[i], q)
            mul!(ṙp, Cp[i], q̇)
        end
    end
end

function update_rbs_states!(tg,q,q̇=zero(q))
    distribute_q_to_rbs!(tg,q,q̇)
    rbs = tg.rigidbodies
    for rbid in tg.mvbodyindex
        rb = rbs[rbid]
        lncs = rb.state.cache.funcs.lncs
        @unpack q, q̇ = rb.state.coords
        R = NaturalCoordinates.find_R(lncs,q)
        Ω = NaturalCoordinates.find_ω(lncs,q,q̇)
        rb.state.R .= R
        # @show Ω
    end
end

function generate_forces!(rbs)
    for (rbid,rb) in enumerate(rbs)
        @unpack state = rb
        @unpack Faps = state
        @unpack Cp,Cg = state.cache
        @unpack Q = state.coords
        Q .= 0.0
        for (pid,f) in enumerate(Faps)
            Q .+= transpose(Cp[pid])*f
        end
        Q .+= transpose(Cg)*state.F
    end
end

<<<<<<< Updated upstream
function assemble_forces!(F,tg;factor=1.0)
    rbs = tg.rigidbodies
    @unpack body2q = tg.connectivity
    generate_forces!(rbs)
    F .= 0.0
    for rbid in tg.mvbodyindex
        pindex = body2q[rbid]
        F[pindex] .+= factor*rbs[rbid].state.coords.Q
    end
end

function assemble_forces(tg;factor=1.0)
    T = get_numbertype(tg)
    @unpack body2q = tg.connectivity
    F = zeros(T,tg.ncoords)
    assemble_forces!(F,tg,factor=factor)
    F
=======
function generate_forces!(tg::ClusterTensegrityStructure)
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

function get_force(tg::AbstractTensegrityStructure)
	tg.state.system.F̌
end

function get_force!(F,tg::AbstractTensegrityStructure)
	F .= get_force(tg)
>>>>>>> Stashed changes
end

function apply_gravity!(tg;factor=1)
    rbs = tg.rigidbodies
    gravity_acceleration = factor*get_gravity(tg)
    for (rbid,rb) in enumerate(rbs)
        @unpack prop, state = rb
        rb.state.F .+= gravity_acceleration*prop.mass
    end
end

function apply_gravity_y!(tgstruct;factor=1)
    rbs = tgstruct.rigidbodies
    gravity_acceleration = factor*get_gravity_y(tgstruct)
    for (rbid,rb) in enumerate(rbs)
        @unpack prop, state = rb
        rb.state.F .+= gravity_acceleration*prop.mass
    end
end

function kinetic_energy_coords(rb::RigidBody)
    @unpack q̇ = rb.state.coords
    @unpack M = rb.state.cache
    ke = 1/2*transpose(q̇)*M*q̇
end

function gravity_potential_energy(rb)
    q = rb.state.coords.q
    gravity_potential_energy(rb,q)
end

function gravity_potential_energy(rb::RigidBody,q)
    @unpack Cg = rb.state.cache
    r = Cg*q
    gravity_acceleration = get_gravity(rb)
    -transpose(r)*gravity_acceleration*rb.prop.mass
end

function potential_energy(s::SString)
    pe = 0.0
    @unpack k,state = s
    Δlen = s.state.length-s.state.restlen
    if Δlen > 0.0
        pe += 1/2*k*Δlen^2
    end
    pe
end

potential_energy(rb::AbstractRigidBody) = gravity_potential_energy(rb)

function kinetic_energy_coords(tg::AbstractTensegrity,q,q̇)
    distribute_q_to_rbs!(tg,q,q̇)
    ke = sum(kinetic_energy_coords.(tg.rigidbodies))
end

function gravity_potential_energy(tg::AbstractTensegrity,q)
    distribute_q_to_rbs!(tg,q)
    sum(gravity_potential_energy.(tg.rigidbodies))
end

function elastic_potential_energy(tg::TensegrityStructure)
    reset_forces!(tg)
    update_strings_apply_forces!(tg)
    pe = sum(potential_energy.(tg.strings))
end

function elastic_potential_energy(tg::TensegrityStructure,q)
    distribute_q_to_rbs!(tg,q)
    elastic_potential_energy(tg)
end

function elastic_potential_energy(bot::TensegrityRobot,q,a)
    actuate!(bot,a)
    elastic_potential_energy(bot.tg,q)
end

<<<<<<< Updated upstream
function energy(tg,q,q̇;gravity=false)
    distribute_q_to_rbs!(tg,q,q̇)
    ke = sum(kinetic_energy_coords.(tg.rigidbodies))
    update_strings_apply_forces!(tg)
    epe = sum(potential_energy.(tg.strings))
    if gravity
        gpe = gravity_potential_energy(tg,q)
    else
        gpe = 0
    end
    ke + epe + gpe
end

function build_body2q(rbs::Vector{rbType}) where rbType<:AbstractRigidBody{N,T,CType} where {N,T,CType}
    lncs = Vector{Vector{T}}()
    bp_number = Vector{Int}()
    push!(bp_number,0)
    body2q = Vector{Vector{Int}}()
    for (rbid,rb) in enumerate(rbs)
        @unpack state = rb
        xi,yi,xj,yj = state.coords.q
        bp1 = [xi,yi]
        bp2 = [xj,yj]
        bp1_find = findfirst(x->x==bp1,lncs)
        if bp1_find === nothing
            push!(lncs,bp1)
            push!(bp_number,bp_number[end]+1)
            bp1_number = bp_number[end]
        else
            bp1_number = bp1_find
        end
        bp2_find = findfirst(x->x==bp2,lncs)
        if bp2_find === nothing
            push!(lncs,bp2)
            push!(bp_number,bp_number[end]+1)
            bp2_number = bp_number[end]
        else
            bp2_number = bp2_find
        end
        push!(body2q,[2bp1_number-1,2bp1_number,
                      2bp2_number-1,2bp2_number])
    end
    body2q
=======
function update!(tg::AbstractTensegrityStructure; gravity=false)
    clear_forces!(tg)
    update_rigids!(tg)
    update_cables!(tg)
    # update_clustercables_apply_forces!(tg)
	if gravity
		apply_gravity!(tg)
	end
	generate_forces!(tg)
>>>>>>> Stashed changes
end

function build_massmatrix(tg::AbstractTensegrity)
    body2q = tg.connectivity.body2q
    ncoords = tg.ncoords
    T = get_numbertype(tg)
    mass_matrix = zeros(T,ncoords,ncoords)
    for rbid in tg.mvbodyindex
        pindex = body2q[rbid]
        mass_matrix[pindex,pindex] .+= tg.rigidbodies[rbid].state.cache.M
    end
    mass_matrix
end

<<<<<<< Updated upstream
function get_nconstraint(rbs,mvbodyindex,nmvbodies,constraints)
    nbodyconstraint = get_nbodyconstraint(rbs)
    nbodydof = get_nbodydof(rbs)
    ninconstraint = nbodyconstraint*nmvbodies
    nexconstraint = 0  #nbodydof*nfixbodies
    foreach(constraints) do cst
        nexconstraint += cst.nconstraints
    end
    nconstraint = ninconstraint + nexconstraint
end

get_nconstraint(tg) = tg.nconstraint

function build_Φ(tg)
    rbs = tg.rigidbodies
    csts = tg.constraints
    #q0,q̇0 = get_q(tg)
    @unpack body2q = tg.connectivity
    nfixbodies = tg.nfixbodies
    nconstraint = tg.nconstraint
    nbodyc = get_nbodyconstraint(tg)
    nbodydof = get_nbodydof(tg)
=======
# function build_M(tg::TensegrityStructure)
function build_M(tg::AbstractTensegrityStructure)
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

function build_M̌(tg::AbstractTensegrityStructure)
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
>>>>>>> Stashed changes
    @inline @inbounds function inner_Φ(q)
        ret = Vector{eltype(q)}(undef,nconstraint)
        is = Ref(0)
        #is[] += nbodydof*nfixbodies
        for rbid in tg.mvbodyindex
            pindex = body2q[rbid]
            rb = rbs[rbid]
            # nc = rb.state.cache.nc
            # if nc > 0
            #     ret[is[]+1:is[]+nc] = rb.state.cache.cfuncs.Φ(q[pindex])
            #     is[] += nc
            # end
            ret[is[]+1:is[]+nbodyc] .= rb.state.cache.funcs.Φ(q[pindex])
            is[] += nbodyc
        end
        foreach(csts) do cst
            nc = cst.nconstraints
            ret[is[]+1:is[]+nc] .= make_Φ(cst)(q)
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

function build_A(tg)
    rbs = tg.rigidbodies
    csts = tg.constraints
    @unpack body2q = tg.connectivity
    nfixbodies = tg.nfixbodies
    nconstraint = tg.nconstraint
    nbodyc = get_nbodyconstraint(tg)
    nbodydof = get_nbodydof(tg)
    ncoords = tg.ncoords
    @inline @inbounds function inner_A(q)
        ret = zeros(eltype(q),nconstraint,ncoords)
        is = Ref(0)
        for rbid in tg.mvbodyindex
            pindex = body2q[rbid]
            rb = rbs[rbid]
            nc = rb.state.cache.nc
            # if nc > 0
            #     ret[is[]+1:is[]+nc,pindex] = rb.state.cache.cfuncs.Φq(q[pindex])
            #     is[] += nc
            # end
            ret[is[]+1:is[]+nbodyc,pindex] .= rb.state.cache.funcs.Φq(q[pindex])
            is[] += nbodyc
        end
        foreach(csts) do cst
            nc = cst.nconstraints
            ret[is[]+1:is[]+nc,:] .= make_A(cst)(q)
            is[] += nc
        end
        ret
    end
end

function build_Aq(tg)
    rbs = tg.rigidbodies
    csts = tg.constraints
    @unpack body2q = tg.connectivity
    nfixbodies = tg.nfixbodies
    nconstraint = tg.nconstraint
    nbodyc = get_nbodyconstraint(tg)
    nbodydof = get_nbodydof(tg)
    ncoords = tg.ncoords
    @inline @inbounds function inner_Aq(λ)
        ret = zeros(Float64,ncoords,ncoords)
        is = Ref(0)
        for rbid in tg.mvbodyindex
            pindex = body2q[rbid]
            rb = rbs[rbid]
            ret[is[]+1:is[]+nbodyc,pindex] .= 2*[1 1 -1 -1;-1 -1 1 1]*diagm(λ[pindex])
            is[] += nbodyc
        end
        ret
    end
end

function build_F(tg,rbid,pid,f)
    rbs = tg.rigidbodies
    Ti = build_Ti(tg,rbid)
    C = rbs[rbid].state.cache.Cp[pid]
    F = transpose(C*Ti)*f
    reshape(F,:,1)
end

function get_q(tg)
    rbs = tg.rigidbodies
    @unpack body2q = tg.connectivity
    ncoords = tg.ncoords
    T = get_numbertype(tg)
    q = zeros(T,ncoords)
    q̇ = zeros(T,ncoords)
    for rbid in tg.mvbodyindex
        pindex = body2q[rbid]
        q[pindex] .= rbs[rbid].state.coords.q
        q̇[pindex] .= rbs[rbid].state.coords.q̇
    end
    return q,q̇
end

function get_q(bot::TensegrityRobots.TensegrityRobot)
    get_q(bot.tg)
end

function get_force(bot::TensegrityRobots.TensegrityRobot)
    get_force(bot, bot.tg.tensiles)
end

function get_force(bot::TensegrityRobots.TensegrityRobot,@eponymargs(strings,clusterstrings))
    @unpack tensiles = bot.tg
    f_list = Vector{Vector{Float64}}()
    push!(f_list, [ss.state.tension for ss in tensiles.strings])
    for clusterstring in tensiles.clusterstrings
        push!(f_list, [cs.state.tension for cs in clusterstring.segs])
    end
    return f_list
end

function get_force(bot::TensegrityRobots.TensegrityRobot,@eponymargs(strings))
    @unpack tensiles = bot.tg
    f_list = Vector{Vector{Float64}}()
    push!(f_list, [ss.state.tension for ss in tensiles.strings])
    return f_list
end

"get restlen length lengthdot tension"
function get_state(bot::TensegrityRobots.TensegrityRobot, dataname)
    get_state(bot, bot.tg.tensiles, dataname)
end

function get_state(bot::TensegrityRobots.TensegrityRobot,@eponymargs(strings,clusterstrings), dataname)
    @unpack tensiles = bot.tg
    f_list = Vector{Vector{Float64}}()
    push!(f_list, [getproperty(ss.state,dataname) for ss in tensiles.strings])
    for clusterstring in tensiles.clusterstrings
        push!(f_list, [getproperty(cs.state,dataname) for cs in clusterstring.segs])
    end
    return f_list
end

function get_state(bot::TensegrityRobots.TensegrityRobot,@eponymargs(strings), dataname)
    @unpack tensiles = bot.tg
    f_list = Vector{Vector{Float64}}()
    push!(f_list, [getproperty(ss.state,dataname) for ss in tensiles.strings])
    return f_list
end

"get μ θ α s s⁺ s⁻"
function get_sp(bot::TensegrityRobots.TensegrityRobot, dataname)
    @unpack clusterstrings = bot.tg
    f_list = Vector{Vector{Any}}()
    for css in clusterstrings
        push!(f_list, getproperty(css.sps,dataname))
    end
    return f_list
end

function csf(bot::TensegrityRobots.TensegrityRobot)
    return bot.tg.clusterstrings
end

function ssf(bot::TensegrityRobots.TensegrityRobot)
    return bot.tg.strings
end

get_λ(tg) = zeros(get_numbertype(tg),tg.nconstraint)

function get_initial(tgstruct)
    q0,q̇0 = get_q(tgstruct)
    λ0 = get_λ(tgstruct)
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

function check_jacobian_singularity(tg)
    q,_ = get_q(tg)
    A = build_A(tg)
    Aq = A(q)
    sys_rank = rank(Aq)
    if sys_rank < minimum(size(Aq))
        @warn "System's Jacobian is singular: rank(A(q))=$(sys_rank)<$(minimum(size(Aq)))"
    end
    for (rbid,rb) in enumerate(tg.rigidbodies)
        if rb.prop.movable && rb.prop.constrained
            q_rb = rb.state.coords.q
            Aq_rb = vcat(rb.state.cache.cfuncs.Φq(q_rb),
                         rb.state.cache.funcs.Φq(q_rb))
            rb_rank = rank(Aq_rb)
            intrinsic_Aq = rb.state.cache.funcs.Φq(q_rb)
            # @show rbid,lucompletepiv!(copy(intrinsic_Aq))
            # col_index = GECP(intrinsic_Aq)
            # @show rbid,col_index
            # @show rank(intrinsic_Aq[:,col_index[1:6]])
            if rb_rank < minimum(size(Aq_rb))
                @warn "The $(rbid)th rigid body's Jacobian is singular: rank(A(q))=$(rb_rank)<$(minimum(size(Aq_rb)))"
            end
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

get_ndim(bot::TensegrityRobot) = get_ndim(bot.tg)
<<<<<<< Updated upstream
get_ndim(tg::AbstractTensegrity) = get_ndim(tg.rigidbodies)
=======
get_ndim(tg::TensegrityStructure) = get_ndim(tg.rigidbodies)
get_ndim(tg::ClusterTensegrityStructure) = get_ndim(tg.rigidbodies)
>>>>>>> Stashed changes
get_ndim(rbs::AbstractVector{<:AbstractRigidBody}) = get_ndim(eltype(rbs))
get_ndim(rb::AbstractRigidBody) = get_ndim(typeof(rb))
get_ndim(::Type{<:AbstractRigidBody{N,T,C}}) where {N,T,C} = N

get_numbertype(bot::TensegrityRobot) = get_numbertype(bot.tg)
<<<<<<< Updated upstream
get_numbertype(tg::AbstractTensegrity) = get_numbertype(tg.rigidbodies)
=======
get_numbertype(tg::TensegrityStructure) = get_numbertype(tg.rigidbodies)
get_numbertype(tg::ClusterTensegrityStructure) = get_numbertype(tg.rigidbodies)
>>>>>>> Stashed changes
get_numbertype(rbs::AbstractVector{<:AbstractRigidBody}) = get_numbertype(eltype(rbs))
get_numbertype(rb::AbstractRigidBody) = get_numbertype(typeof(rb))
get_numbertype(::Type{<:AbstractRigidBody{N,T,C}}) where {N,T,C} = T

get_nbodyconstraint(bot::TensegrityRobot) = get_nbodyconstraint(bot.tg)
get_nbodyconstraint(tg::AbstractTensegrity) = get_nbodyconstraint(tg.rigidbodies)
get_nbodyconstraint(rbs::AbstractVector{<:AbstractRigidBody}) = get_nbodyconstraint(eltype(rbs))
get_nbodyconstraint(rb::AbstractRigidBody) = get_nbodyconstraint(typeof(rb))
get_nbodyconstraint(::Type{<:RigidBody{N,T,L,C,
                <:NaturalCoordinatesCache{ArrayT,MT,
                <:NaturalCoordinates.CoordinateFunctions{lncsType},
                cfT}}}) where {N,T,L,C,ArrayT,MT,lncsType,cfT} = NaturalCoordinates.get_nconstraint(lncsType)

get_nbodycoords(bot::TensegrityRobot) = get_nbodycoords(bot.tg)
get_nbodycoords(tg::AbstractTensegrity) = get_nbodycoords(tg.rigidbodies)
get_nbodycoords(rbs::AbstractVector{<:AbstractRigidBody}) = get_nbodycoords(eltype(rbs))
get_nbodycoords(rb::AbstractRigidBody) = get_nbodycoords(typeof(rb))
get_nbodycoords(::Type{<:RigidBody{N,T,L,C,
                <:NaturalCoordinatesCache{ArrayT,MT,
                <:NaturalCoordinates.CoordinateFunctions{lncsType},
                cfT}}}) where {N,T,L,C,ArrayT,MT,lncsType,cfT} = NaturalCoordinates.get_ncoords(lncsType)

get_nbodydof(bot::TensegrityRobot) = get_nbodydof(bot.tg)
get_nbodydof(tg::AbstractTensegrity) = get_nbodydof(tg.rigidbodies)
get_nbodydof(rbs::AbstractVector{<:AbstractRigidBody}) = get_nbodydof(eltype(rbs))
get_nbodydof(rb::AbstractRigidBody) = get_nbodydof(typeof(rb))
get_nbodydof(::Type{<:AbstractRigidBody{2,T,C}}) where {T,C} = 3
get_nbodydof(::Type{<:AbstractRigidBody{3,T,C}}) where {T,C} = 6

get_gravity(bot::TensegrityRobot) = get_gravity(bot.tg)
get_gravity(tg::AbstractTensegrity) = get_gravity(tg.rigidbodies)
get_gravity(rbs::AbstractVector{<:AbstractRigidBody}) = get_gravity(eltype(rbs))
get_gravity(rb::AbstractRigidBody) = get_gravity(typeof(rb))
get_gravity(::Type{<:AbstractRigidBody{2,T,C}}) where {T,C} = [zero(T),-9.81*one(T)]
get_gravity(::Type{<:AbstractRigidBody{3,T,C}}) where {T,C} = [zero(T),zero(T),-9.81*one(T)]

<<<<<<< Updated upstream

get_gravity_y(tr::TensegrityRobot) = get_gravity_y(tr.tg)
get_gravity_y(tg::AbstractTensegrity) = get_gravity_y(tg.rigidbodies)
get_gravity_y(rbs::AbstractVector{<:AbstractRigidBody}) = get_gravity_y(eltype(rbs))
get_gravity_y(rb::AbstractRigidBody) = get_gravity_y(typeof(rb))
get_gravity_y(::Type{<:AbstractRigidBody{2,T,C}}) where {T,C} = [zero(T),9.81*one(T)]
get_gravity_y(::Type{<:AbstractRigidBody{3,T,C}}) where {T,C} = [zero(T),9.81*one(T),zero(T)]
=======
function get_rigidbodies(tg::TensegrityStructure)
	sort_rigidbodies(tg.rigidbodies)
end
function get_rigidbodies(tg::ClusterTensegrityStructure)
	sort_rigidbodies(tg.rigidbodies)
end

function get_rigidbars(tg::TensegrityStructure)
	rbs = get_rigidbodies(tg)
	[rb for rb in rbs
	if rb.state.cache.funcs.lncs isa Union{NaturalCoordinates.LNC2D4C,NaturalCoordinates.LNC3D6C}]
end
function get_rigidbars(tg::ClusterTensegrityStructure)
	rbs = get_rigidbodies(tg)
	[rb for rb in rbs
	if rb.state.cache.funcs.lncs isa Union{NaturalCoordinates.LNC2D4C,NaturalCoordinates.LNC3D6C}]
end
>>>>>>> Stashed changes

get_strings_len(bot::TensegrityRobot) = get_strings_len(bot.tg)
get_strings_deform(bot::TensegrityRobot) = get_strings_deform(bot.tg)
get_strings_restlen(bot::TensegrityRobot) = get_strings_restlen(bot.tg)
get_strings_len_dot(bot::TensegrityRobot) = get_strings_len_dot(bot.tg)
get_strings_tension(bot::TensegrityRobot) = get_strings_tension(bot.tg)
get_strings_stiffness(bot::TensegrityRobot) = get_strings_stiffness(bot.tg)

function get_strings_len!(tg::AbstractTensegrity,q)
    distribute_q_to_rbs!(tg,q,zero(q))
    update_strings_apply_forces!(tg)
    get_strings_len(tg)
end

function get_strings_stiffness(tg::AbstractTensegrity)
    [s.k for s in tg.strings]
end

function get_strings_len(tg::AbstractTensegrity)
    [s.state.length for s in tg.strings]
end

function get_strings_len_dot(tg::AbstractTensegrity)
    [s.state.lengthdot for s in tg.strings]
end

function get_strings_deform(tg::AbstractTensegrity)
    [s.state.length - s.state.restlen for s in tg.strings]
end

function get_strings_restlen(tg::AbstractTensegrity)
    [s.state.restlen for s in tg.strings]
end

function get_strings_tension(tg::AbstractTensegrity)
    [s.state.tension for s in tg.strings]
end

function get_original_restlen(botinput::TensegrityRobot)
    bot = deepcopy(botinput)
    T = get_numbertype(bot)
    actuate!(bot,zeros(T,length(bot.hub.actuators)))
    u0 = get_strings_restlen(bot.tg)
end

function force_densities_to_restlen(tg::AbstractTensegrity,γs)
    [
    begin
        l = s.state.length
        l̇ = s.state.lengthdot
        k = s.k
        c = s.c
        u = l-(γ*l-c*l̇)/k
    end
        for (γ,s) in zip(γs,tg.strings)]
end

function find_remaining_index(body2q,rbs)
    original_nq = maximum(maximum.(body2q))
    switch_index = zeros(Int,original_nq)
    for (rbid,rb) in enumerate(rbs)
        qindex = body2q[rbid]
        if rb.prop.movable
            for i in qindex
                switch_index[i] = i
            end
        end
    end
    remaining_index = findall((x)->x!=0,switch_index)
end

function filter_body2q(rbs)
    body2q_raw = build_body2q(rbs)
    body2q = filter_body2q(body2q_raw,rbs)
end

function filter_body2q(body2q,rbs)
    original_nq = maximum(maximum.(body2q))
    remaining_index = find_remaining_index(body2q,rbs)
    qpointer = collect(1:original_nq)[remaining_index]
    filtered_body2q = Vector{Vector{Int}}()
    for (rbid,rb) in enumerate(rbs)
        qindex = body2q[rbid]
        filtered_index = zero(qindex)
        if rb.prop.movable
            for (j,i) in enumerate(qindex)
                filtered_index[j] = findfirst((x)->x==i,qpointer)
            end
        end
        push!(filtered_body2q,filtered_index)
    end
    filtered_body2q
end

function build_Y(bot)
	@unpack tg, hub = bot
	@unpack actuators = hub
    @unpack nstrings,strings = tg
    nact = length(actuators)
    ret = spzeros(Int,nstrings,nact)
    for (i,iact) in enumerate(actuators)
		if typeof(iact)<:ManualActuator
			is1 = iact.reg.id_string
	        ret[is1,i] = 1
		elseif typeof(iact)<:ManualGangedActuators
	        is1, is2 = iact.regs.id_strings
	        ret[is1,i] = 1
	        ret[is2,i] = -1
		else
			error("Unknown actuator type")
		end
    end
    ret
end

function new_trajectory(tg::TensegrityStructure)
    t0 = zero(get_numbertype(tg))
    q0, q̇0  = get_initial(tg)
    λ0 = get_λ(tg)
    ConstrainedCoordinatesTrajectory([t0], [q0], [q̇0], [λ0])
end

function new_trajectory(tg::ClusterTensegrityStructure)
    t0 = zero(get_numbertype(tg))
    q0, q̇0  = get_initial(tg)
    λ0 = get_λ(tg)
    s̄0 = get_s̄(tg)
    SlidingConstrainedCoordinatesTrajectory([t0], [q0], [q̇0], [λ0], [s̄0])
end

function record_trajectory(tg::ClusterTensegrityStructure)
    t0 = zero(get_numbertype(tg))
    q0, q̇0  = get_initial(tg)
    λ0 = get_λ(tg)
    s̄0 = get_s̄(tg)
    SlidingConstrainedCoordinatesTrajectoryRecordData([t0], [q0], [q̇0], [λ0], [s̄0], Vector{Float64}(), [])
end

function TensegrityRobot(tg,hub)
	reset_forces!(tg)
    # update_strings_apply_forces!(tg)
	# check_jacobian_singularity(tg)
	# check_stability(tg)
    TensegrityRobot(tg,hub,new_trajectory(tg))
end

function TensegrityRobotRecord(tg,hub)
	reset_forces!(tg)
    # update_strings_apply_forces!(tg)
	# check_jacobian_singularity(tg)
	# check_stability(tg)
    TensegrityRobot(tg,hub,record_trajectory(tg))
end

function reset!(bot::TensegrityRobot)
    @unpack tg, traj = bot
    reset!(tg,traj)
    reset!(traj)
end

#function reset!(bot::TensegrityRobotRecord)
#    @unpack tg, traj = bot
#    reset!(tg,traj)
#    reset!(traj)
#end

function reset!(tg::TensegrityStructure,traj)
    @unpack qs,q̇s = traj
    reset_forces!(tg)
    distribute_q_to_rbs!(tg,qs[begin],q̇s[begin])
    update_strings!(tg)
end

function reset!(tg::ClusterTensegrityStructure,traj)
    @unpack qs,q̇s,s̄s = traj
    reset_forces!(tg)
    reset_restlen!(tg)
    distribute_q_to_rbs!(tg,qs[begin],q̇s[begin])
    distribute_s̄!(tg,s̄s[begin])
    update_strings!(tg)
end

function reset!(traj::ConstrainedCoordinatesTrajectory)
    @unpack ts, qs, q̇s, λs = traj
    resize!(ts,1)
    resize!(qs,1)
    resize!(q̇s,1)
    resize!(λs,1)
end

function reset!(traj::SlidingConstrainedCoordinatesTrajectory)
    @unpack ts,qs,q̇s,λs,s̄s= traj
    resize!(ts,1)
    resize!(qs,1)
    resize!(q̇s,1)
    resize!(λs,1)
    resize!(s̄s,1)
end

function reset!(traj::SlidingConstrainedCoordinatesTrajectoryRecordData)
    @unpack ts,qs,q̇s,λs,s̄s,iterations,OtherData= traj
    resize!(ts,1)
    resize!(qs,1)
    resize!(q̇s,1)
    resize!(λs,1)
    resize!(s̄s,1)
    resize!(iterations,0)
    resize!(OtherData,0)
end

function reset_restlen!(tg::ClusterTensegrityStructure)
    @unpack clusterstrings = tg
    for clusterstring in clusterstrings
        for seg in clusterstring.segs
            seg.state.restlen = seg.original_restlen
        end
    end
end

<<<<<<< Updated upstream
function set_new_initial!(bot::TensegrityRobot,q,q̇=zero(q))
    @unpack tg, traj = bot
    traj.qs[begin] .= q
    traj.q̇s[begin] .= q̇
    reset!(bot)
=======
"""
更新系统到指定时间步状态。
$(TYPEDSIGNATURES)
"""
function goto_step!(bot::TensegrityRobot,that_step;actuate=false)
	(;tg, traj) = bot
    tg.state.system.q .= traj.q[that_step]
    tg.state.system.q̇ .= traj.q̇[that_step]
	if actuate
		actuate!(bot,[traj.t[that_step]])
	end
	update!(tg)
	bot
end

function analyse_slack(tg::AbstractTensegrityStructure,verbose=false)
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
>>>>>>> Stashed changes
end
