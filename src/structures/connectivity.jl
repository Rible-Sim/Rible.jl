struct Apparatus{jointType,forceType}
    id::Int
    joint::jointType
    force::forceType
end

function Base.isless(a::Apparatus,b::Apparatus)
    isless(a.id,b.id)
end

get_id(appar::Apparatus) = appar.id

function connect(bodies, spring_dampers, cm_input=Int[;;], istart = 0)
    _,nb = check_id_sanity(bodies)
    if size(cm_input,2) > nb
        @warn "Cropping the connecting matrix."
        cm = cm_input[:,1:nb]
    else
        cm = cm_input[:,:]
    end
    rbs_sorted = sort(bodies)
    ret = []
    j = 0
    for row in eachrow(cm)
        rbids = findall(!iszero,row)
        if isempty(rbids)
            continue
        end
        @assert length(rbids) == 2
        @assert reduce(*,row[rbids]) < 0
        rbid1,rbid2 = ifelse(row[rbids[1]]>0,rbids,reverse(rbids))
        pid1,pid2 = Int64.(abs.(row[[rbid1,rbid2]]))
        j += 1
        joint = CableJoint(
            Hen2Egg(ID(rbs_sorted[rbid1],pid1),ID(rbs_sorted[rbid2],pid2)),
            0,
        )
        force = spring_dampers[j]
        cable = Apparatus(
            istart+j,
            joint,
            force
        )
        push!(ret,cable)
    end
    ret
end

function cluster(bodies, cm2_input)
    rbs_sorted = sort(bodies)
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

function connect_and_cluster(bodies, cm_input, cm2_input)
    ret1 = connect(bodies, cm_input)
    ret2 = cluster(bodies, cm2_input)
    return (connected=ret1, clustered=ret2)
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
    sys_free_coords_idx::idxType
    sys_pres_coords_idx::idxType
    bodyid2sys_full_coords::id2idxType
    bodyid2sys_free_coords::id2idxType
    bodyid2sys_pres_coords::id2idxType
    bodyid2sys_intrinsic_cstr_idx::id2idxType
    bodyid2sys_dof_idx::id2idxType
    apparid2sys_extrinsic_cstr_idx::id2idxType
    apparid2full_idx::id2idxType
    apparid2free_idx::id2idxType
    apparid2sys_free_coords_idx::id2idxType
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

function index(bodies,apparatuses,sharing_input::AbstractMatrix=Int[;;])
    bodies_ids,num_of_bodies = check_id_sanity(bodies)
    apparatuses_ids,num_of_apparatuses = check_id_sanity(apparatuses)
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
    # num_of_extrinsic_cstr = mapreduce((apparatus)->apparatus.num_of_cstr,+,apparatuses,init=0)
    apparid2sys_extrinsic_cstr_idx = Vector{Int}[]
    nexcst_by_joint = zeros(Int,num_of_apparatuses)
    num_of_joint_apparatuses = 0
    num_of_force_apparatuses = 0
    apparid2full_idx = Vector{Vector{Int}}(undef,num_of_apparatuses)
    apparid2free_idx = Vector{Vector{Int}}(undef,num_of_apparatuses)
    apparid2sys_free_coords_idx = Vector{Vector{Int}}(undef,num_of_apparatuses)
    apparid2sys_extrinsic_cstr_idx = Vector{Int}[]
    ilast = 0
    if num_of_apparatuses > 0
        foreach(apparatuses) do apparatus
            nexcst_by_joint[apparatus.id] = apparatus.joint.num_of_cstr
            (;id,joint) = apparatus
            apparid2full_idx[id], apparid2free_idx[id], apparid2sys_free_coords_idx[id] = 
            get_joint_idx(joint,bodyid2sys_free_coords)
            if !(apparatus.joint isa Nothing)
                num_of_joint_apparatuses += 1
            end
            if !(apparatus.force isa Nothing)
                num_of_force_apparatuses += 1
            end
        end
    end
    num_of_extrinsic_cstr = sum(nexcst_by_joint)
    num_of_cstr = num_of_intrinsic_cstr + num_of_extrinsic_cstr
    num_of_dof = num_of_free_coords - num_of_cstr
    if num_of_dof <= 0
        @warn "Non positive degree of freedom: $num_of_dof."
    end
    for apparid = 1:num_of_apparatuses
        nexcst = nexcst_by_joint[apparid]
        push!(apparid2sys_extrinsic_cstr_idx,collect(ilast+1:ilast+nexcst))
        ilast += nexcst
    end

    Indexed(
        num_of_bodies,
        num_of_apparatuses,
        num_of_joint_apparatuses,
        num_of_force_apparatuses,
        length(sysfull),
        num_of_free_coords,
        num_of_pres_coords,
        num_of_intrinsic_cstr,
        num_of_extrinsic_cstr,
        num_of_cstr,
        num_of_dof_unconstrained,
        num_of_dof,
        sys_free_coords_idx,
        sys_pres_coords_idx,
        bodyid2sys_full_coords,
        bodyid2sys_free_coords,
        bodyid2sys_pres_coords,
        bodyid2sys_intrinsic_cstr_idx,
        bodyid2sys_dof_idx,
        apparid2sys_extrinsic_cstr_idx,
        apparid2full_idx,
        apparid2free_idx,
        apparid2sys_free_coords_idx
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

function number(bodies,apparatuses)
    _,nb = check_id_sanity(bodies)
    nnodes_by_body = zeros(Int,nb)
    nld_by_body = zeros(Int,nb)
    foreach(bodies) do body
        i = body.prop.id
        nnodes_by_body[i] = body.prop.loci |> length
        nld_by_body[i] = get_num_of_local_dims(body)
    end
    bodyid2sys_loci_idx = Vector{Int}[]
    num2ID = Vector{ID{Int,Int}}()
    sys_loci2coords_idx = Vector{Int}[]
    is = 1
    js = 0
    for bodyid in 1:nb
        push!(bodyid2sys_loci_idx,Int[])
        nld = nld_by_body[bodyid]
        for pid in 1:nnodes_by_body[bodyid]
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


