function index_incstr(bodies)
    ids,num_of_bodies = check_id_sanity(bodies)
    nincst_by_body = zeros(Int,num_of_bodies)
    nintrinsic_by_body = zeros(Int,num_of_bodies)
    ndof_by_body = zeros(Int,num_of_bodies)
    foreach(bodies) do body
        bodyid = body.prop.id
        # Active intrinsic constraints participate in system constraint assembly.
        nincst_by_body[bodyid] = get_num_of_cstr(body.coords)
        # Intrinsic nullspace DOF uses the full intrinsic model (ignores cstr_idx customization).
        nintrinsic_by_body[bodyid] = get_num_of_intrinsic_cstr(body)
        ndof_by_body[bodyid] = get_num_of_coords(body) - nintrinsic_by_body[bodyid]
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

"""
Standard Rigid Body Connectivity Type.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Connectivity{id2idxType} <: AbstractConnectivity
    num_of_bodies::Int
    num_of_apparatuses::Int
    num_of_joint_apparatuses::Int
    num_of_force_apparatuses::Int
    num_of_full_coords::Int
    num_of_intrinsic_cstr::Int
    num_of_extrinsic_cstr::Int
    num_of_cstr::Int
    num_of_dof_unconstrained::Int
    num_of_dof::Int
    num_of_aux_var::Int
    sys_full_coords_idx::Vector{Int}
    bodyid2sys_full_coords::id2idxType
    bodyid2sys_intrinsic_cstr_idx::id2idxType
    bodyid2sys_dof_idx::id2idxType
    apparid2sys_extrinsic_cstr_idx::id2idxType
    apparid2sys_aux_var_idx::id2idxType
    "body's loci' idx to System's loci' idx"
    bodyid2sys_locus_id::Vector{Vector{Int}}
    "System's loci' idx to Signifier"
    sys_locus_id2sig::Vector{Signifier{Int64}}
    "System's loci' idx to System's loci' coords' idx"
    sys_locus_id2coords_idx::Vector{Vector{Int}}
    "body's loci' idx to System's loci' coords' idx"
    bodyid2sys_loci_coords_idx::Vector{Vector{Int}}
    "Number of the System's loci"
    num_of_sys_loci::Int
    "Number of the System's loci' coords"
    num_of_sys_loci_coords::Int
    num_of_body_params::Int
    num_of_appar_params::Int
    num_of_params::Int
    apparid2params_idx::Vector{Vector{Int}}
    apparid2sys_full_coords_idx::Vector{Vector{Int}}
end

function process_apparatus_indices(apparatuses, bodyid2sys_full_coords, num_of_intrinsic_cstr)
    _, num_of_apparatuses = check_id_sanity(apparatuses)
    num_of_excst_by_joint = zeros(Int,num_of_apparatuses)
    num_of_aux_var_by_appar = zeros(Int,num_of_apparatuses)
    num_of_joint_apparatuses = 0
    num_of_force_apparatuses = 0
    apparid2sys_full_coords_idx = Vector{Vector{Int}}(undef,num_of_apparatuses)
    apparid2sys_aux_var_idx = Vector{Vector{Int}}(undef,num_of_apparatuses)
    apparid2sys_extrinsic_cstr_idx = Vector{Vector{Int}}(undef,num_of_apparatuses)
    num_of_aux_var_last = 0
    num_excst_last = 0
    if num_of_apparatuses > 0
        foreach(apparatuses) do apparatus
            id = apparatus.id
            num_of_aux_var_by_appar[id] = apparatus.num_of_aux_var
            apparid2sys_full_coords_idx[id] = get_appar_idx(apparatus,bodyid2sys_full_coords)
            apparid2sys_aux_var_idx[id] = collect(num_of_aux_var_last+1:num_of_aux_var_last+num_of_aux_var_by_appar[id])
            num_of_aux_var_last += num_of_aux_var_by_appar[id]
            num_of_excst_by_joint[id] = apparatus.joint.num_of_cstr
            apparid2sys_extrinsic_cstr_idx[id] = num_of_intrinsic_cstr .+ collect(num_excst_last+1:num_excst_last+num_of_excst_by_joint[id])
            num_excst_last += num_of_excst_by_joint[id]
            if !(apparatus.joint isa NoJoint)
                num_of_joint_apparatuses += 1
            end
            if !(apparatus.force isa NoForce)
                num_of_force_apparatuses += 1
            end
        end
    end
    num_of_extrinsic_cstr = sum(num_of_excst_by_joint)
    num_of_aux_var = sum(num_of_aux_var_by_appar)
    @eponymtuple(
        num_of_joint_apparatuses,
        num_of_force_apparatuses,
        num_of_extrinsic_cstr,
        num_of_aux_var,
        apparid2sys_full_coords_idx,
        apparid2sys_aux_var_idx,
        apparid2sys_extrinsic_cstr_idx
    )
end

function process_loci_indices(bodies)
    _,nb = check_id_sanity(bodies)
    nloci_by_body = zeros(Int,nb)
    nld_by_body = zeros(Int,nb)
    foreach(bodies) do body
        i = body.prop.id
        nloci_by_body[i] = body.prop.loci |> length
        nld_by_body[i] = get_num_of_local_dims(body)
    end
    bodyid2sys_locus_id = Vector{Int}[]
    sys_locus_id2sig = Vector{Signifier{Int}}()
    sys_locus_id2coords_idx = Vector{Int}[]
    num_of_sys_loci = 1
    num_of_sys_loci_coords = 0
    for bodyid in 1:nb
        push!(bodyid2sys_locus_id,Int[])
        nld = nld_by_body[bodyid]
        for pid in 1:nloci_by_body[bodyid]
            push!(sys_locus_id2sig,Signifier(bodyid,pid))
            push!(bodyid2sys_locus_id[bodyid],num_of_sys_loci)
            num_of_sys_loci += 1
            push!(sys_locus_id2coords_idx,collect(1:nld).+num_of_sys_loci_coords)
            num_of_sys_loci_coords += nld
        end
    end
    bodyid2sys_loci_coords_idx = [
        reduce(vcat,sys_locus_id2coords_idx[bodyid2sys_locus_id[bodyid]])
        for bodyid = 1:nb
    ]
    @eponymtuple(
        bodyid2sys_locus_id,
        sys_locus_id2sig,
        sys_locus_id2coords_idx,
        bodyid2sys_loci_coords_idx,
        num_of_sys_loci,
        num_of_sys_loci_coords
    )
end

function process_param_indices(apparatuses, num_of_body_params)
    _, num_of_apparatuses = check_id_sanity(apparatuses)
    apparid2params_per_appar = Vector{Int}(undef,num_of_apparatuses)
    apparid2params_idx = Vector{Int}[]
    foreach(apparatuses) do appar
        nc = get_num_of_params(appar)
        apparid2params_per_appar[appar.id] = nc
    end
    num_of_appar_params = 0
    for apparid = 1:num_of_apparatuses
        push!(apparid2params_idx,collect(1:apparid2params_per_appar[apparid]).+num_of_body_params.+num_of_appar_params)
        num_of_appar_params += apparid2params_per_appar[apparid]
    end
    num_of_params = num_of_body_params + num_of_appar_params
    @eponymtuple(
        num_of_appar_params,
        num_of_params,
        apparid2params_idx
    )
end


function process_coords_idx(bodies; sharing_matrix=Int[;;])
    bodies_ids,num_of_bodies = check_id_sanity(bodies)
    if size(sharing_matrix,2) > num_of_bodies
        @warn "Cropping the sharing matrix."
        sharing = sharing_matrix[:,1:num_of_bodies]
    else
        sharing = sharing_matrix[:,:]
    end
    sys_full_coords_idx = Int[]
    bodyid2sys_full_coords = Vector{Int}[]
    ntotal_by_body = zeros(Int,num_of_bodies)
    foreach(bodies) do body
        bodyid = body.prop.id
        ntotal_by_body[bodyid] = get_num_of_coords(body)
    end
    for bodyid = 1:num_of_bodies
        ntotal = ntotal_by_body[bodyid]
        push!(bodyid2sys_full_coords,fill(-1,ntotal))
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
        deleteat!(unshareds,sort(shared_idx))
        nusi = length(unshareds)
        bodyid2sys_full_coords[bodyid][unshareds] = collect(length(sys_full_coords_idx)+1:length(sys_full_coords_idx)+nusi)
        append!(sys_full_coords_idx,bodyid2sys_full_coords[bodyid][unshareds])
    end
    num_of_full_coords = length(sys_full_coords_idx)
    @eponymtuple(
        num_of_bodies, sys_full_coords_idx, num_of_full_coords,
        bodyid2sys_full_coords
    )
end

"""
Constructor for the standard `Connectivity` type.
$(TYPEDSIGNATURES)
"""
function Connectivity(bodies,apparatuses=Int[];sharing_matrix=Int[;;])
    apparatuses_ids,num_of_apparatuses = check_id_sanity(apparatuses)
    num_of_intrinsic_cstr,bodyid2sys_intrinsic_cstr_idx,num_of_dof_unconstrained,bodyid2sys_dof_idx = index_incstr(bodies)
    (;
        num_of_bodies, sys_full_coords_idx,num_of_full_coords, bodyid2sys_full_coords 
    ) = process_coords_idx(bodies;sharing_matrix)

    (;
        num_of_joint_apparatuses,
        num_of_force_apparatuses,
        num_of_extrinsic_cstr,
        num_of_aux_var,
        apparid2sys_aux_var_idx,
        apparid2sys_extrinsic_cstr_idx,
        apparid2sys_full_coords_idx
    ) = process_apparatus_indices(apparatuses, bodyid2sys_full_coords, num_of_intrinsic_cstr)

    num_of_cstr = num_of_intrinsic_cstr + num_of_extrinsic_cstr
    num_of_dof = num_of_full_coords - num_of_cstr
    if num_of_dof <= 0
        @warn "Non positive degree of freedom: $num_of_dof."
    end

    (;
        bodyid2sys_locus_id,
        sys_locus_id2sig,
        sys_locus_id2coords_idx,
        bodyid2sys_loci_coords_idx,
        num_of_sys_loci,
        num_of_sys_loci_coords
    ) = process_loci_indices(bodies)

    num_of_body_params = num_of_sys_loci_coords

    (;
        num_of_appar_params,
        num_of_params,
        apparid2params_idx
    ) = process_param_indices(apparatuses, num_of_body_params)

    Connectivity(
        num_of_bodies,
        num_of_apparatuses,
        num_of_joint_apparatuses,
        num_of_force_apparatuses,
        num_of_full_coords,
        num_of_intrinsic_cstr,
        num_of_extrinsic_cstr,
        num_of_cstr,
        num_of_dof_unconstrained,
        num_of_dof,
        num_of_aux_var,
        sys_full_coords_idx,
        bodyid2sys_full_coords,
        bodyid2sys_intrinsic_cstr_idx,
        bodyid2sys_dof_idx,
        apparid2sys_extrinsic_cstr_idx,
        apparid2sys_aux_var_idx,
        bodyid2sys_locus_id,
        sys_locus_id2sig,
        sys_locus_id2coords_idx,
        bodyid2sys_loci_coords_idx,
        num_of_sys_loci,
        num_of_sys_loci_coords,
        num_of_body_params,
        num_of_appar_params,
        num_of_params,
        apparid2params_idx,
        apparid2sys_full_coords_idx
    )
end

get_num_of_coords(cnt::Connectivity) = cnt.num_of_full_coords
get_num_of_free_coords(cnt::Connectivity) = cnt.num_of_full_coords
get_num_of_pres_coords(cnt::Connectivity) = 0
get_free_coords_idx(cnt::Connectivity) = cnt.sys_full_coords_idx
get_pres_coords_idx(cnt::Connectivity) = Int[]
