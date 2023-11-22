function build_joint_cache(
        nmcs_hen::NC,
        nmcs_egg::NC,
        loci_position_hen,
        loci_position_egg,
        loci_axes_trl_hen,
        loci_axes_trl_egg,
        loci_axes_rot_hen,
        loci_axes_rot_egg,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg
    )
    c_hen = to_local_coords(nmcs_hen,loci_position_hen)
    c_egg = to_local_coords(nmcs_egg,loci_position_egg)
    invX̄_hen = nmcs_hen.data.invX̄
    invX̄_egg = nmcs_egg.data.invX̄
    axes_trl_hen = invX̄_hen*loci_axes_trl_hen
    axes_trl_egg = invX̄_egg*loci_axes_trl_egg
    axes_rot_egg = invX̄_egg*loci_axes_rot_egg
    num_of_coords_hen = get_num_of_coords(nmcs_hen)
    num_of_coords_egg = get_num_of_coords(nmcs_egg)
    num_of_jointed_coords = num_of_coords_hen + num_of_coords_egg
    # translate
    T = eltype(q_hen)
    heros = spzeros(T,num_of_jointed_coords,num_of_jointed_coords)
    # half 1st
    half_1st = fill(heros,3)
    # half 4th
    half_2nd = fill(heros,1)
    # half 3rd
    half_3rd = fill(heros,6)
    # half 2nd
    half_4th = fill(heros,3)
    # translate

    # translate
    C_hen = to_transformation(nmcs_hen,q_hen,c_hen)
    C_egg = to_transformation(nmcs_egg,q_egg,c_egg)
    J = [-C_hen C_egg;] |> sparse
    transformations = J[mask_1st,:]

    I3_Bool = I(3)

    half_2nd[1] = (J'*J) |> sparse

    X_hen = get_X(nmcs_hen,q_hen)
    X_egg = get_X(nmcs_egg,q_egg)
    select_uvw_hen = BlockDiagonal([nmcs_hen.conversion_to_X,zero(nmcs_egg.conversion_to_X)])
    select_uvw_egg = BlockDiagonal([zero(nmcs_hen.conversion_to_X),nmcs_egg.conversion_to_X])
    q = vcat(q_hen,q_egg)
    if (nmcs_hen isa NC3D12C) && (nmcs_egg isa NC3D12C)
        # hes 3rd on hen
        # translate on hen
        for i = 1:3
            axis_hen = axes_trl_hen.X[:,i]
            axis_zero = zero(axis_hen)
            half_3rd[i] = select_uvw_hen'*kron(vcat(0,axis_hen,0,axis_zero),I3_Bool)*J |> sparse
        end
        # hes 3rd on egg
        # translate on egg
        for i = 1:3
            axis_egg = axes_trl_egg.X[:,i]
            axis_zero = zero(axis_egg)
            half_3rd[3+i] = select_uvw_egg'*kron(vcat(0,axis_zero,0,axis_egg),I3_Bool)*J |> sparse
        end
        # hes 2nd
        # rotate of egg
        axes_rot_hen = inv(X_hen)*X_egg*axes_rot_egg
        axes_idx = [
            (2,3),
            (2,1),
            (3,1)
        ]
        for (i,(id_axis_hen,id_axis_egg)) in enumerate(axes_idx)
            axis_hen = axes_rot_hen.X[:,id_axis_hen]
            axis_egg = axes_rot_egg.X[:,id_axis_egg]
            axis_zero = zero(axis_egg)
            half_4th[i] = select_uvw_hen'*
                kron(vcat(0,axis_hen,0,axis_zero),I3_Bool)*
                kron(vcat(0,axis_zero,0,axis_egg),I3_Bool)'*
                select_uvw_egg |> sparse
        end
    else
        # rotational joint not to be used for bars
        axes_rot_hen = axes_rot_egg
    end

    halves = vcat(
        half_1st[mask_1st],
        half_2nd[mask_2nd],
        half_3rd[mask_3rd],
        half_4th[mask_4th];
    )
    # cstr values
    Refq = Ref(q)
    RefqT = Ref(q')
    values = RefqT.*halves.*Refq
    values[mask_1st] .+= transformations*q


    # valid for NC only, wrong for QC
    hess_1st = [(H .+ H') |> Symmetric for H in half_1st]
    hess_2nd = [(H .+ H') |> Symmetric for H in half_2nd]
    hess_3rd = [(H .+ H') |> Symmetric for H in half_3rd]
    hess_4th = [(H .+ H') |> Symmetric for H in half_4th]
    hessians = vcat(
        hess_1st[mask_1st],
        hess_2nd[mask_2nd],
        hess_3rd[mask_3rd],
        hess_4th[mask_4th];
    )

    cache = @eponymtuple(
        transformations,halves,hessians
    )
    cache, values
end

function get_joint_violations!(
        ret,
        nmcs_hen::NC, 
        nmcs_egg::NC,
        loci_position_hen,
        loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg,
        violations
    )
    (;transformations,halves,) = cache
    q = vcat(q_hen,q_egg)
    for (icstr,half) = enumerate(halves)
        ret[icstr] = q'*half*q
    end
    ret[mask_1st] .+= transformations*q
    ret .-= violations
end

function get_joint_jacobian!(
        ret,
        nmcs_hen::NC,nmcs_egg::NC,
        loci_position_hen,
        loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        free_idx,
        q_hen,q_egg
    )
    (;transformations,hessians) = cache
    q = vcat(q_hen,q_egg)
    num_of_coords_hen = get_num_of_coords(nmcs_hen)
    for (icstr,hess) = enumerate(hessians)
        A = q'*hess
        ret[icstr,:] .= A[1,free_idx]
    end
    ret[mask_1st,:] .+= transformations[mask_1st,free_idx]
end

function get_joint_forces_jacobian!(
        ret,
        num_of_cstr,
        nmcs_hen::NC, nmcs_egg::NC,
        loci_position_hen,
        loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        free_idx,
        q_hen,q_egg,
        λ
    )
    (;hessians) = cache
    ret .= 0.0
    for i = 1:num_of_cstr
        ret .+= -λ[i] .* hessians[i][free_idx,free_idx]
    end
end
