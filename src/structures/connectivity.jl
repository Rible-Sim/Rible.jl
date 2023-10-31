
function check_rbid_sanity(rbs)
    ids,nb = get_bodies_ids(rbs)
    @assert minimum(ids) == 1
    @assert maximum(ids) == nb
    @assert allunique(ids)
    ids,nb
end

struct Numbered
    "body's loci' idx to System's loci' idx"
    bodyid2sys_loci_idx::Vector{Vector{Int}}
    "System's loci' idx to System's loci' coords' idx"
    sys_loci2coords_idx::Vector{Vector{Int}}
    "body's loci' idx to System's loci' coords' idx"
    bodyid2sys_loci_coords_idx::Vector{Vector{Int}}
    "Number of the System's loci' coords"
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
    bodyid2sys_loci_idx = Vector{Int}[]
    num2ID = Vector{ID{Int,Int}}()
    sys_loci2coords_idx = Vector{Int}[]
    is = 1
    js = 0
    for bodyid in 1:nb
        push!(bodyid2sys_loci_idx,Int[])
        nld = nld_by_mem[bodyid]
        for pid in 1:nnodes_by_mem[bodyid]
            push!(num2ID,ID(bodyid,pid))
            push!(bodyid2sys_loci_idx[bodyid],is)
            is += 1
            push!(sys_loci2coords_idx,collect(1:nld).+js)
            js += nld
        end
    end
    bodyid2sys_loci_coords_idx = [
        reduce(vcat,sys_loci2coords_idx[bodyid2sys_loci_idx[bodyid]])
        for bodyid = 1:nb
    ]
    # @show bodyid2sys_loci_coords_idx
    Numbered(
        bodyid2sys_loci_idx,
        sys_loci2coords_idx,
        bodyid2sys_loci_coords_idx,
        js
    )
end

struct Indexed{mem2sysType,sysType}
    num_of_full_coords::Int
    num_of_free_coords::Int
    num_of_pres_coords::Int
    num_of_bodies::Int
    bodyid2sys_full_coords::mem2sysType
    bodyid2sys_free_coords::mem2sysType
    bodyid2sys_pres_coords::mem2sysType
    sys_free_coords_idx::sysType
    sys_pres_coords_idx::sysType
    num_of_intrinsic_cstr::Int
    bodyid2sys_intrinsic_cstr_idx::mem2sysType
    sys_num_of_dof::Int
    bodyid2sys_dof_idx::mem2sysType
end

function index_incstr(rbs)
    ids,num_of_bodies = check_rbid_sanity(rbs)
    nincst_by_mem = zeros(Int,num_of_bodies)
    ndof_by_mem = zeros(Int,num_of_bodies)
    foreach(rbs) do body
        nincst_by_mem[body.prop.id] = body.state.cache.num_of_cstr
        ndof_by_mem[body.prop.id] = get_num_of_dof(body)
    end
    num_of_intrinsic_cstr = sum(nincst_by_mem)
    sys_num_of_dof = sum(ndof_by_mem)
    bodyid2sys_intrinsic_cstr_idx = Vector{Int}[]
    bodyid2sys_dof_idx = Vector{Int}[]
    ilast = 0
    jlast = 0
    for bodyid = 1:num_of_bodies
        nincst = nincst_by_mem[bodyid]
        push!(bodyid2sys_intrinsic_cstr_idx,collect(ilast+1:ilast+nincst))
        ilast += nincst
        num_of_dof = ndof_by_mem[bodyid]
        push!(bodyid2sys_dof_idx,collect(jlast+1:jlast+num_of_dof))
        jlast += num_of_dof 
    end
    num_of_intrinsic_cstr,bodyid2sys_intrinsic_cstr_idx,sys_num_of_dof,bodyid2sys_dof_idx
end

function index(rbs,sharing_input=Int[;;])
    ids,num_of_bodies = check_rbid_sanity(rbs)
    if size(sharing_input,2) > num_of_bodies
        @warn "Cropping the sharing matrix."
        sharing = sharing_input[:,1:num_of_bodies]
    else
        sharing = sharing_input[:,:]
    end
    sysfull = Int[]
    sys_pres_coords_idx = Int[]
    sys_free_coords_idx = Int[]
    bodyid2sys_full_coords = Vector{Int}[]
    bodyid2sys_pres_coords = Vector{Int}[]
    bodyid2sys_free_coords = Vector{Int}[]
    ntotal_by_mem = zeros(Int,num_of_bodies)
    pres_idx_by_mem = Vector{Vector{Int}}(undef,num_of_bodies)
    free_idx_by_mem = Vector{Vector{Int}}(undef,num_of_bodies)
    foreach(rbs) do body
        bodyid = body.prop.id
        ntotal_by_mem[bodyid] = get_num_of_coords(body)
        pres_idx_by_mem[bodyid] = body.state.cache.pres_idx
        free_idx_by_mem[bodyid] = body.state.cache.free_idx
    end
    for bodyid = 1:num_of_bodies
        ntotal = ntotal_by_mem[bodyid]
        pres = pres_idx_by_mem[bodyid]
        free = free_idx_by_mem[bodyid]
        num_of_pres_coords = length(pres)
        num_of_free_coords = ntotal - num_of_pres_coords
        push!(bodyid2sys_full_coords,fill(-1,ntotal))
        push!(bodyid2sys_pres_coords,Int[])
        push!(bodyid2sys_free_coords,Int[])
        unshareds = collect(1:ntotal)
        shared_idx = Int[]
        for row in eachrow(sharing)
            rbids = findall(!iszero,row)
            if bodyid in rbids[begin+1:end]
                myindex = row[bodyid]
                formerid = first(rbids)
                formerindex = row[formerid]
                bodyid2sys_full_coords[bodyid][myindex] = bodyid2sys_full_coords[formerid][formerindex]
                push!(shared_idx,myindex)
            end
        end
        deleteat!(unshareds,shared_idx)
        nusi = length(unshareds)
        bodyid2sys_full_coords[bodyid][unshareds] = collect(length(sysfull)+1:length(sysfull)+nusi)
        append!(sysfull,bodyid2sys_full_coords[bodyid][unshareds])
        for i in unshareds
            if i in pres
                # pres
                push!(sys_pres_coords_idx,bodyid2sys_full_coords[bodyid][i])
            else
                # free
                push!(sys_free_coords_idx,bodyid2sys_full_coords[bodyid][i])
            end
        end
        for i in free
            freei = findfirst((x)->x==bodyid2sys_full_coords[bodyid][i],sys_free_coords_idx)
            push!(bodyid2sys_free_coords[bodyid],freei)
        end
        for i in pres
            presi = findfirst((x)->x==bodyid2sys_full_coords[bodyid][i],sys_pres_coords_idx)
            push!(bodyid2sys_pres_coords[bodyid],presi)
        end
    end
    num_of_intrinsic_cstr,bodyid2sys_intrinsic_cstr_idx,sys_num_of_dof,bodyid2sys_dof_idx = index_incstr(rbs)
    Indexed(
        length(sysfull),length(sys_free_coords_idx),length(sys_pres_coords_idx),num_of_bodies,
        bodyid2sys_full_coords,bodyid2sys_free_coords,bodyid2sys_pres_coords,
        sys_free_coords_idx,sys_pres_coords_idx,
        num_of_intrinsic_cstr,bodyid2sys_intrinsic_cstr_idx,
        sys_num_of_dof,bodyid2sys_dof_idx
    )
end

struct JointedMembers{JType,joint2sysexcstType}
    njoints::Int
    num_of_extrinsic_cstr::Int
    joints::JType
    joint2sysexcst::joint2sysexcstType
end

function unjoin()
    njoints = 0
    joints = Int[]
    num_of_extrinsic_cstr = 0
    joint2sysexcst = Vector{Int}[]
    JointedMembers(njoints,num_of_extrinsic_cstr,joints,joint2sysexcst)
end

function join(joints,indexed)
    # num_of_extrinsic_cstr = mapreduce((joint)->joint.num_of_cstr,+,joints,init=0)
    njoints = length(joints)
    joint2sysexcst = Vector{Int}[]
    nexcst_by_joint = zeros(Int,njoints)
    foreach(joints) do joint
        nexcst_by_joint[joint.id] = joint.num_of_cstr
    end
    num_of_extrinsic_cstr = sum(nexcst_by_joint)
    joint2sysexcst = Vector{Int}[]
    ilast = 0
    for jointid = 1:njoints
        nexcst = nexcst_by_joint[jointid]
        push!(joint2sysexcst,collect(ilast+1:ilast+nexcst))
        ilast += nexcst
    end
    JointedMembers(njoints,num_of_extrinsic_cstr,joints,joint2sysexcst)
end

function Base.isless(rb1::AbstractBody,rb2::AbstractBody)
    isless(rb1.prop.id,rb2.prop.id)
end

function Base.isless(e2e1::Hen2Egg,e2e2::Hen2Egg)
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
        push!(ret_raw,Hen2Egg(is,ID(rbs_sorted[rbid1],pid1),ID(rbs_sorted[rbid2],pid2)))
    end
    ret = TypeSortedCollection(ret_raw)
end

function cluster(rbs, cm2_input)
    rbs_sorted = sort(rbs)
    ret_raw = []
    cm = cm2_input
    for row in eachrow(cm)
        iret = Vector{Hen2Egg}()
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
            push!(iret,Hen2Egg(is,ID(rbs_sorted[rbid1],pid1),ID(rbs_sorted[rbid2],pid2)))
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

