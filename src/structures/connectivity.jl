
function check_rbid_sanity(rbs)
    ids,nb = get_bodies_ids(rbs)
    @assert minimum(ids) == 1
    @assert maximum(ids) == nb
    @assert allunique(ids)
    ids,nb
end

struct NumberedPoints
    "Member's points' indices to System's points' indices"
    mem2num::Vector{Vector{Int}}
    "System's points' indices to System's points' coordinates' indices"
    num2sys::Vector{Vector{Int}}
    "Member's points' indices to System's points' coordinates' indices"
    mem2sys::Vector{Vector{Int}}
    "Number of the System's points' coodinates"
    nc::Int
end

function number(rbs)
    _,nb = check_rbid_sanity(rbs)
    nnodes_by_mem = zeros(Int,nb)
    nld_by_mem = zeros(Int,nb)
    foreach(rbs) do body
        i = body.prop.id
        nnodes_by_mem[i] = body.prop.loci |> length
        nld_by_mem[i] = get_num_of_local_dims(body)
    end
    mem2num = Vector{Int}[]
    num2ID = Vector{ID{Int,Int}}()
    num2sys = Vector{Int}[]
    is = 1
    js = 0
    for bodyid in 1:nb
        push!(mem2num,Int[])
        nld = nld_by_mem[bodyid]
        for pid in 1:nnodes_by_mem[bodyid]
            push!(num2ID,ID(bodyid,pid))
            push!(mem2num[bodyid],is)
            is += 1
            push!(num2sys,collect(1:nld).+js)
            js += nld
        end
    end
    mem2sys = [
        reduce(vcat,num2sys[mem2num[bodyid]])
        for bodyid = 1:nb
    ]
    # @show mem2sys
    NumberedPoints(mem2num,num2sys,mem2sys,js)
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
    sysndof::Int
    mem2sysndof::mem2sysType
end

function index_inconstraints(rbs)
    ids,nmem = check_rbid_sanity(rbs)
    nincst_by_mem = zeros(Int,nmem)
    ndof_by_mem = zeros(Int,nmem)
    foreach(rbs) do body
        nincst_by_mem[body.prop.id] = body.state.cache.nΦ
        ndof_by_mem[body.prop.id] = get_num_of_dof(body)
    end
    ninconstraints = sum(nincst_by_mem)
    sysndof = sum(ndof_by_mem)
    mem2sysincst = Vector{Int}[]
    mem2sysndof = Vector{Int}[]
    ilast = 0
    jlast = 0
    for bodyid = 1:nmem
        nincst = nincst_by_mem[bodyid]
        push!(mem2sysincst,collect(ilast+1:ilast+nincst))
        ilast += nincst
        ndof = ndof_by_mem[bodyid]
        push!(mem2sysndof,collect(jlast+1:jlast+ndof))
        jlast += ndof 
    end
    ninconstraints,mem2sysincst,sysndof,mem2sysndof
end

function index(rbs,sharing_input=Int[;;])
    ids,nmem = check_rbid_sanity(rbs)
    if size(sharing_input,2) > nmem
        @warn "Cropping the sharing matrix."
        sharing = sharing_input[:,1:nmem]
    else
        sharing = sharing_input[:,:]
    end
    sysfull = Int[]
    syspres = Int[]
    sysfree = Int[]
    mem2sysfull = Vector{Int}[]
    mem2syspres = Vector{Int}[]
    mem2sysfree = Vector{Int}[]
    ntotal_by_mem = zeros(Int,nmem)
    pres_idx_by_mem = Vector{Vector{Int}}(undef,nmem)
    free_idx_by_mem = Vector{Vector{Int}}(undef,nmem)
    foreach(rbs) do body
        bodyid = body.prop.id
        ntotal_by_mem[bodyid] = get_nbodycoords(body)
        pres_idx_by_mem[bodyid] = body.state.cache.pres_idx
        free_idx_by_mem[bodyid] = body.state.cache.free_idx
    end
    for bodyid = 1:nmem
        ntotal = ntotal_by_mem[bodyid]
        pres = pres_idx_by_mem[bodyid]
        free = free_idx_by_mem[bodyid]
        npres = length(pres)
        nfree = ntotal - npres
        push!(mem2sysfull,fill(-1,ntotal))
        push!(mem2syspres,Int[])
        push!(mem2sysfree,Int[])
        unshareds = collect(1:ntotal)
        shared_indices = Int[]
        for row in eachrow(sharing)
            rbids = findall(!iszero,row)
            if bodyid in rbids[begin+1:end]
                myindex = row[bodyid]
                formerid = first(rbids)
                formerindex = row[formerid]
                mem2sysfull[bodyid][myindex] = mem2sysfull[formerid][formerindex]
                push!(shared_indices,myindex)
            end
        end
        deleteat!(unshareds,shared_indices)
        nusi = length(unshareds)
        mem2sysfull[bodyid][unshareds] = collect(length(sysfull)+1:length(sysfull)+nusi)
        append!(sysfull,mem2sysfull[bodyid][unshareds])
        for i in unshareds
            if i in pres
                # pres
                push!(syspres,mem2sysfull[bodyid][i])
            else
                # free
                push!(sysfree,mem2sysfull[bodyid][i])
            end
        end
        for i in free
            freei = findfirst((x)->x==mem2sysfull[bodyid][i],sysfree)
            push!(mem2sysfree[bodyid],freei)
        end
        for i in pres
            presi = findfirst((x)->x==mem2sysfull[bodyid][i],syspres)
            push!(mem2syspres[bodyid],presi)
        end
    end
    ninconstraints,mem2sysincst,sysndof,mem2sysndof = index_inconstraints(rbs)
    IndexedMemberCoords(
        length(sysfull),length(sysfree),length(syspres),nmem,
        mem2sysfull,mem2sysfree,mem2syspres,
        sysfree,syspres,
        ninconstraints,mem2sysincst,
        sysndof,mem2sysndof
    )
end

struct JointedMembers{JType,joint2sysexcstType}
    njoints::Int
    nexconstraints::Int
    joints::JType
    joint2sysexcst::joint2sysexcstType
end

function unjoin()
    njoints = 0
    joints = Int[]
    nexconstraints = 0
    joint2sysexcst = Vector{Int}[]
    JointedMembers(njoints,nexconstraints,joints,joint2sysexcst)
end

function join(joints,indexed)
    # nexconstraints = mapreduce((joint)->joint.nconstraints,+,joints,init=0)
    njoints = length(joints)
    joint2sysexcst = Vector{Int}[]
    nexcst_by_joint = zeros(Int,njoints)
    foreach(joints) do joint
        nexcst_by_joint[joint.id] = joint.nconstraints
    end
    nexconstraints = sum(nexcst_by_joint)
    joint2sysexcst = Vector{Int}[]
    ilast = 0
    for jointid = 1:njoints
        nexcst = nexcst_by_joint[jointid]
        push!(joint2sysexcst,collect(ilast+1:ilast+nexcst))
        ilast += nexcst
    end
    JointedMembers(njoints,nexconstraints,joints,joint2sysexcst)
end

function Base.isless(rb1::AbstractBody,rb2::AbstractBody)
    isless(rb1.prop.id,rb2.prop.id)
end

function Base.isless(e2e1::End2End,e2e2::End2End)
    isless(e2e1.id,e2e2.id)
end

function Base.sort(c::TypeSortedCollection)
    sort!(reduce(vcat,c.data))
end

function connect(rbs,cm_input=Int[;;])
    _,nb = check_rbid_sanity(rbs)
    if size(cm_input,2) > nb
        @warn "Cropping the connecting matrix."
        cm = cm_input[:,1:nb]
    else
        cm = cm_input[:,:]
    end
    rbs_sorted = sort(rbs)
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
    rbs_sorted = sort(rbs)
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
Rigid Body Connectivity Type.
$(TYPEDEF)
"""
struct Connectivity{numberType,indexType,tensionType,jointType}
    numbered::numberType
    indexed::indexType
    tensioned::tensionType
    jointed::jointType
end

"""
Connectivity Constructor.
$(TYPEDSIGNATURES)
"""
function Connectivity(numbered,indexed,tensioned)
    Connectivity(numbered,indexed,tensioned,unjoin())
end

