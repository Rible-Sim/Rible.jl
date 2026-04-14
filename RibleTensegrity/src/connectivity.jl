
function change_connecting_format(bodies, connecting_matrix=Int[;;],)
    _, nb = check_id_sanity(bodies)
    if size(connecting_matrix, 2) > nb
        @warn "Cropping the connecting matrix."
        cm = connecting_matrix[:, 1:nb]
    else
        cm = connecting_matrix[:, :]
    end
    ret = Matrix{Int}(undef, size(cm, 1), 4)
    
    for (i, row) in enumerate(eachrow(cm))
        rbids = findall(!iszero, row)
        if isempty(rbids)
            continue
        end
        @assert length(rbids) == 2
        @assert reduce(*, row[rbids]) < 0
        bid1, bid2 = ifelse(row[rbids[1]] > 0, rbids, reverse(rbids))
        pid1, pid2 = Int64.(abs.(row[[bid1, bid2]]))
        ret[i,:] = [bid1, pid1, bid2, pid2]
    end
    ret
end

function connect_spring(bodies, spring_dampers; cm=Int[;;], istart=0)
    @assert size(cm, 2) == 4 "The connecting matrix should have 4 columns with the format: bid1, pid1, bid1, pid2."
    rbs_sorted = sort(bodies)
    ret = []
    j = 0
    for row in eachrow(cm)
        bid1, pid1, bid2, pid2 = row
        j += 1
        joint = CableJoint(
            Hen2Egg(Signifier(rbs_sorted[bid1], pid1), Signifier(rbs_sorted[bid2], pid2)),
            0,
        )
        full_coords_idx = get_joint_idx(joint)
        force = spring_dampers[j]
        cable = Apparatus(
            istart + j,
            joint,
            force,
            0, #num_of_aux_var,
            full_coords_idx,
            full_coords_idx,
            ApparatusCache(joint,force,full_coords_idx,)
        )
        push!(ret, cable)
    end
    ret
end

function connect_clusters(bodies, cluster_sps, cluster_segs, cms)
    _, nb = check_id_sanity(bodies)
    rets = []
    id = 0
    nsegs = 0
    for (cluster_seg, cm) in zip(cluster_segs, cms)
        @assert size(cm, 2) == 4
        rbs_sorted = sort(bodies)
        j = 0
        id += 1
        begin_id = 0
        ret = [
            begin
                bid1, pid1, bid2, pid2 = row
                j += 1
                joint = CableJoint(
                    Hen2Egg(Signifier(rbs_sorted[bid1], pid1), Signifier(rbs_sorted[bid2], pid2)),
                    0,
                )
                (;hen,egg) = joint.hen2egg
                coords_hen = hen.body.coords
                coords_egg = egg.body.coords
                ncoords_hen = get_num_of_coords(coords_hen)
                ncoords_egg = get_num_of_coords(coords_egg)
                full_coords_idx = vcat(
                    collect(1:ncoords_hen),
                    collect(1:ncoords_egg) .+ ncoords_hen
                ) .+ begin_id
                begin_id+ncoords_hen, full_coords_idx, full_coords_idx
                force = cluster_seg[j]
                segs = Apparatus(
                    j,
                    joint,
                    force,
                    nsegs,
                    full_coords_idx,
                    full_coords_idx,
                    ApparatusCache(joint,force,full_coords_idx,)
                )
                segs
            end
            for row in eachrow(cm)
        ] |> ClusterDistanceSpringDampers
        nsegs += size(cm, 1)
        joint = ClusterJoint(cluster_sps[id], 0)
        push!(rets, 
            Apparatus(
                id, 
                joint, #joint
                ret, #force
                2j-2, #num_of_aux_var
                Int[], #full_coords_idx
                Int[], #full_coords_idx
                ApparatusCache(joint,ret,Int[],)
            )
        )
    end
    rets
end

function connect_spring_and_clusters(bodies, spring_dampers, 
        cluster_sps, cluster_segs, 
        connecting_matrix, connecting_cluster_matrix; 
        istart=0
    )
    ret1 = connect_spring(bodies, spring_dampers; cm=connecting_matrix, istart)
    ret2 = connect_clusters(bodies, cluster_sps, cluster_segs, connecting_cluster_matrix)
    ret1, ret2
end
