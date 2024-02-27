function connect_spring(bodies, spring_dampers; cm=Int[;;], istart=0)
    _, nb = check_id_sanity(bodies)
    @assert size(cm, 2) == 4
    rbs_sorted = sort(bodies)
    ret = []
    j = 0
    for row in eachrow(cm)
        rbid1, pid1, rbid2, pid2 = row
        j += 1
        joint = CableJoint(
            Hen2Egg(Signifier(rbs_sorted[rbid1], pid1), Signifier(rbs_sorted[rbid2], pid2)),
            0,
        )
        full_coords_idx, free_coords_idx = get_joint_idx(joint)
        force = spring_dampers[j]
        cable = Apparatus(
            istart + j,
            joint,
            force,
            0,
            full_coords_idx,
            free_coords_idx
        )
        push!(ret, cable)
    end
    ret
end

function connect_clusters(bodies, cluster_sps, cluster_segs, cms)
    _, nb = check_id_sanity(bodies)
    rets = []
    clusterID = 0
    num_segs = 0
    for (cluster_seg, cm) in zip(cluster_segs, cms)
        @assert size(cm, 2) == 4
        rbs_sorted = sort(bodies)
        #TODO 这里ret的类型是Vector{Any}， 导致下面的Apparatus不是concrete， 会影响运算性能， 建议修改
        ret = []
        j = 0
        clusterID += 1
        for row in eachrow(cm)
            rbid1, pid1, rbid2, pid2 = row
            j += 1
            joint = CableJoint(
                Hen2Egg(Signifier(rbs_sorted[rbid1], pid1), Signifier(rbs_sorted[rbid2], pid2)),
                0,
            )
            segs = Apparatus(
                j,
                joint,
                cluster_seg[j],
                num_segs,
                Int[],
                Int[]
            )
            push!(ret, segs)
        end
        num_segs += size(cm, 1)
        push!(rets, Apparatus(clusterID, ClusterJoint(cluster_sps[clusterID], 0), ret, 2j-2, Int[], Int[]))
    end
    rets
end

function connect_spring_and_clusters(bodies, spring_dampers, cluster_sps, cluster_segs, connecting_matrix, connecting_cluster_matrix; istart=0)
    ret1 = connect_spring(bodies, spring_dampers; cm=connecting_matrix, istart)
    ret2 = connect_clusters(bodies, cluster_sps, cluster_segs, connecting_cluster_matrix)
    ret1, ret2
end

function connect(bodies, spring_dampers; connecting_matrix=Int[;;], istart=0)
    _, nb = check_id_sanity(bodies)
    if size(connecting_matrix, 2) > nb
        @warn "Cropping the connecting matrix."
        cm = connecting_matrix[:, 1:nb]
    else
        cm = connecting_matrix[:, :]
    end
    rbs_sorted = sort(bodies)
    ret = []
    j = 0
    for row in eachrow(cm)
        rbids = findall(!iszero, row)
        if isempty(rbids)
            continue
        end
        @assert length(rbids) == 2
        @assert reduce(*, row[rbids]) < 0
        rbid1, rbid2 = ifelse(row[rbids[1]] > 0, rbids, reverse(rbids))
        pid1, pid2 = Int64.(abs.(row[[rbid1, rbid2]]))
        j += 1
        joint = CableJoint(
            Hen2Egg(Signifier(rbs_sorted[rbid1], pid1), Signifier(rbs_sorted[rbid2], pid2)),
            0,
        )
        full_coords_idx, free_coords_idx = get_joint_idx(joint)
        force = spring_dampers[j]
        cable = Apparatus(
            istart + j,
            joint,
            force,
            0,
            full_coords_idx,
            free_coords_idx
        )
        push!(ret, cable)
    end
    ret
end



struct Indexed{id2idxType,idxType}
    num_of_bodies::Int
    num_of_apparatuses::Int
    num_of_joint_apparatuses::Int
    num_of_force_apparatuses::Int
    num_of_full_coords::Int
    num_of_free_coords::Int
    num_of_pres_coords::Int
    num_of_intrinsic_cstr::Int
    num_of_extrinsic_cstr::Int
    num_of_cstr::Int
    num_of_dof_unconstrained::Int
    num_of_dof::Int
    num_of_add_var::Int
    sys_free_coords_idx::idxType
    sys_pres_coords_idx::idxType
    bodyid2sys_full_coords::id2idxType
    bodyid2sys_free_coords::id2idxType
    bodyid2sys_pres_coords::id2idxType
    bodyid2sys_intrinsic_cstr_idx::id2idxType
    bodyid2sys_dof_idx::id2idxType
    apparid2sys_extrinsic_cstr_idx::id2idxType
    apparid2sys_free_coords_idx::id2idxType
    apparid2sys_add_var_idx::id2idxType
end

function index_incstr(bodies)
    ids,num_of_bodies = check_id_sanity(bodies)
    nincst_by_body = zeros(Int,num_of_bodies)
    ndof_by_body = zeros(Int,num_of_bodies)
    foreach(bodies) do body
        nincst_by_body[body.prop.id] = length(body.coords.cstr_idx)
        ndof_by_body[body.prop.id] = get_num_of_dof(body)
    end
    num_of_intrinsic_cstr = sum(nincst_by_body)
    num_of_dof_unconstrained = sum(ndof_by_body)
    bodyid2sys_intrinsic_cstr_idx = Vector{Int}[]
    bodyid2sys_dof_idx = Vector{Int}[]
    ilast = 0
    jlast = 0
    for bodyid = 1:num_of_bodies
        nincst = nincst_by_body[bodyid]
        push!(bodyid2sys_intrinsic_cstr_idx,collect(ilast+1:ilast+nincst))
        ilast += nincst
        num_of_dof = ndof_by_body[bodyid]
        push!(bodyid2sys_dof_idx,collect(jlast+1:jlast+num_of_dof))
        jlast += num_of_dof 
    end
    num_of_intrinsic_cstr,bodyid2sys_intrinsic_cstr_idx,num_of_dof_unconstrained,bodyid2sys_dof_idx
end

function index(bodies,apparatuses=Int[];sharing_matrix::AbstractMatrix=Int[;;])
    bodies_ids,num_of_bodies = check_id_sanity(bodies)
    apparatuses_ids,num_of_apparatuses = check_id_sanity(apparatuses)
    if size(sharing_matrix,2) > num_of_bodies
        @warn "Cropping the sharing matrix."
        sharing = sharing_matrix[:,1:num_of_bodies]
    else
        sharing = sharing_matrix[:,:]
    end
    sys_full_coords_idx = Int[]
    sys_pres_coords_idx = Int[]
    sys_free_coords_idx = Int[]
    bodyid2sys_full_coords = Vector{Int}[]
    bodyid2sys_pres_coords = Vector{Int}[]
    bodyid2sys_free_coords = Vector{Int}[]
    ntotal_by_body = zeros(Int,num_of_bodies)
    pres_idx_by_body = Vector{Vector{Int}}(undef,num_of_bodies)
    free_idx_by_body = Vector{Vector{Int}}(undef,num_of_bodies)
    foreach(bodies) do body
        bodyid = body.prop.id
        ntotal_by_body[bodyid] = get_num_of_coords(body)
        pres_idx_by_body[bodyid] = body.coords.pres_idx
        free_idx_by_body[bodyid] = body.coords.free_idx
    end
    for bodyid = 1:num_of_bodies
        ntotal = ntotal_by_body[bodyid]
        pres = pres_idx_by_body[bodyid]
        free = free_idx_by_body[bodyid]
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
        bodyid2sys_full_coords[bodyid][unshareds] = collect(length(sys_full_coords_idx)+1:length(sys_full_coords_idx)+nusi)
        append!(sys_full_coords_idx,bodyid2sys_full_coords[bodyid][unshareds])
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
            free_idx = findfirst((x)->x==bodyid2sys_full_coords[bodyid][i],sys_free_coords_idx)
            push!(bodyid2sys_free_coords[bodyid],free_idx)
        end
        for i in pres
            pres_idx = findfirst((x)->x==bodyid2sys_full_coords[bodyid][i],sys_pres_coords_idx)
            push!(bodyid2sys_pres_coords[bodyid],pres_idx)
        end
    end
    num_of_intrinsic_cstr,bodyid2sys_intrinsic_cstr_idx,num_of_dof_unconstrained,bodyid2sys_dof_idx = index_incstr(bodies)
    num_of_free_coords = length(sys_free_coords_idx)
    num_of_pres_coords = length(sys_pres_coords_idx)
    num_of_full_coords = length(sys_full_coords_idx)
    # num_of_extrinsic_cstr = mapreduce((apparatus)->apparatus.num_of_cstr,+,apparatuses,init=0)
    num_of_excst_by_joint = zeros(Int,num_of_apparatuses)
    num_of_add_var_by_appar = zeros(Int,num_of_apparatuses)
    num_of_joint_apparatuses = 0
    num_of_force_apparatuses = 0
    apparid2sys_free_coords_idx = Vector{Vector{Int}}(undef,num_of_apparatuses)
    apparid2sys_add_var_idx = Vector{Vector{Int}}(undef,num_of_apparatuses)
    apparid2sys_extrinsic_cstr_idx = Vector{Vector{Int}}(undef,num_of_apparatuses)
    num_of_add_var_last = 0
    num_excst_last = 0
    if num_of_apparatuses > 0
        foreach(apparatuses) do apparatus
            id = apparatus.id
            num_of_add_var_by_appar[id] = apparatus.num_of_add_var
            apparid2sys_free_coords_idx[id] = get_appar_idx(apparatus,bodyid2sys_free_coords)
            apparid2sys_add_var_idx[id] = collect(num_of_add_var_last+1:num_of_add_var_last+num_of_add_var_by_appar[id])
            num_of_add_var_last += num_of_add_var_by_appar[id]
            num_of_excst_by_joint[id] = apparatus.joint.num_of_cstr
            apparid2sys_extrinsic_cstr_idx[id] = collect(num_excst_last+1:num_excst_last+num_of_excst_by_joint[id])
            num_excst_last += num_of_excst_by_joint[id]
            if !(apparatus.joint isa Nothing)
                num_of_joint_apparatuses += 1
            end
            if !(apparatus.force isa Nothing)
                num_of_force_apparatuses += 1
            end
        end
    end
    num_of_extrinsic_cstr = sum(num_of_excst_by_joint)
    num_of_add_var = sum(num_of_add_var_by_appar)
    num_of_cstr = num_of_intrinsic_cstr + num_of_extrinsic_cstr
    num_of_dof = num_of_free_coords - num_of_cstr
    if num_of_dof <= 0
        @warn "Non positive degree of freedom: $num_of_dof."
    end

    Indexed(
        num_of_bodies,
        num_of_apparatuses,
        num_of_joint_apparatuses,
        num_of_force_apparatuses,
        num_of_full_coords,
        num_of_free_coords,
        num_of_pres_coords,
        num_of_intrinsic_cstr,
        num_of_extrinsic_cstr,
        num_of_cstr,
        num_of_dof_unconstrained,
        num_of_dof,
        num_of_add_var,
        sys_free_coords_idx,
        sys_pres_coords_idx,
        bodyid2sys_full_coords,
        bodyid2sys_free_coords,
        bodyid2sys_pres_coords,
        bodyid2sys_intrinsic_cstr_idx,
        bodyid2sys_dof_idx,
        apparid2sys_extrinsic_cstr_idx,
        apparid2sys_free_coords_idx,
        apparid2sys_add_var_idx
    )
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

function number(bodies,apparatuses=nothing)
    _,nb = check_id_sanity(bodies)
    nnodes_by_body = zeros(Int,nb)
    nld_by_body = zeros(Int,nb)
    foreach(bodies) do body
        i = body.prop.id
        nnodes_by_body[i] = body.prop.loci |> length
        nld_by_body[i] = get_num_of_local_dims(body)
    end
    bodyid2sys_loci_idx = Vector{Int}[]
    num2ID = Vector{Signifier{Int,Int}}()
    sys_loci2coords_idx = Vector{Int}[]
    is = 1
    js = 0
    for bodyid in 1:nb
        push!(bodyid2sys_loci_idx,Int[])
        nld = nld_by_body[bodyid]
        for pid in 1:nnodes_by_body[bodyid]
            push!(num2ID,Signifier(bodyid,pid))
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

"""
Rigid Body Connectivity Type.
$(TYPEDEF)
"""
struct Connectivity{indexType,numberType}
    indexed::indexType
    numbered::numberType
end

"""
Connectivity Constructor.
$(TYPEDSIGNATURES)
"""
function Connectivity(numbered::Numbered,indexed::Indexed)
    Connectivity(indexed,numbered,)
end


