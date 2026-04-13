function build_joint_cache(
        coords_hen::NC{Nh,Mh,Th,Lh,NCH,NCH2},
        coords_egg::NC{Ne,Me,Te,Le,NCE,NCE2},
        loci_position_hen,
        loci_position_egg,
        loci_axes_trl_hen,
        loci_axes_trl_egg,
        loci_axes_rot_hen,
        loci_axes_rot_egg,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg
    ) where {Nh,Mh,Th,Lh,NCH,NCH2,Ne,Me,Te,Le,NCE,NCE2}
    c_hen = to_local_coords(coords_hen,loci_position_hen)
    c_egg = to_local_coords(coords_egg,loci_position_egg)
    invX̄_hen = coords_hen.data.invX̄
    invX̄_egg = coords_egg.data.invX̄
    axes_trl_hen = invX̄_hen*loci_axes_trl_hen
    axes_trl_egg = invX̄_egg*loci_axes_trl_egg
    axes_rot_egg = invX̄_egg*loci_axes_rot_egg
    T = promote_type(eltype(q_hen), eltype(q_egg))
    num_of_coords_hen = get_num_of_coords(coords_hen)
    num_of_coords_egg = get_num_of_coords(coords_egg)
    num_of_jointed_coords = num_of_coords_hen + num_of_coords_egg
    joint_q = Vector{T}(undef, num_of_jointed_coords)
    joint_work = similar(joint_q)
    trf_work = Vector{T}(undef,length(mask_1st))
    # translate
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
    C_hen = to_position_jacobian(coords_hen,q_hen,c_hen)
    C_egg = to_position_jacobian(coords_egg,q_egg,c_egg)
    J = [-C_hen C_egg;] |> sparse
    transformations = J[mask_1st,:]

    I3_Bool = I(3)

    half_2nd[1] = (J'*J) |> sparse

    X_hen = get_X(coords_hen,q_hen)
    X_egg = get_X(coords_egg,q_egg)
    select_uvw_hen = BlockDiagonal([coords_hen.conversion_to_X,zero(coords_egg.conversion_to_X)])
    select_uvw_egg = BlockDiagonal([zero(coords_hen.conversion_to_X),coords_egg.conversion_to_X])
    q = vcat(q_hen,q_egg)
    if (coords_hen isa NC3D12C) && (coords_egg isa NC3D12C)
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
            (3,1), #normal * tangent
            (3,2), #normal * bitangent
            (1,2)  #tangent * bitangent
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
    relative_core = axes_rot_egg.X*axes_rot_hen.X'

    halves = vcat(
        half_1st[mask_1st],
        half_2nd[mask_2nd],
        half_3rd[mask_3rd],
        half_4th[mask_4th];
    )
    # cstr violations
    Refq = Ref(q)
    RefqT = Ref(q')
    violations = RefqT.*halves.*Refq
    violations[mask_1st] .+= transformations*q
    hess_1st = [(H .+ H') for H in half_1st]
    hess_2nd = [(H .+ H') for H in half_2nd]
    hess_3rd = [(H .+ H') for H in half_3rd]
    hess_4th = [(H .+ H') for H in half_4th]
    hessians = vcat(
        hess_1st[mask_1st],
        hess_2nd[mask_2nd],
        hess_3rd[mask_3rd],
        hess_4th[mask_4th];
    )
    cache = @eponymtuple(
        relative_core,
        transformations,halves,hessians,
        joint_q,joint_work,trf_work,
        num_of_coords_hen,num_of_coords_egg
    )
    cache, violations
end

function assemble_joint_q!(joint_q::AbstractVector, q_hen, q_egg, num_of_coords_hen, num_of_coords_egg)
    copyto!(joint_q, 1, q_hen, 1, num_of_coords_hen)
    copyto!(joint_q, num_of_coords_hen+1, q_egg, 1, num_of_coords_egg)
    joint_q
end

function get_joint_violations!(
        ret,
        coords_hen::NC, 
        coords_egg::NC,
        loci_position_hen,
        loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg,
        violations
    )
    (;transformations,halves,joint_q,joint_work,trf_work,num_of_coords_hen,num_of_coords_egg) = cache
    assemble_joint_q!(joint_q, q_hen, q_egg, num_of_coords_hen, num_of_coords_egg)
    # for (icstr,half) = enumerate(halves)
    #     ret[icstr] = q'*half*q
    # end
    # ret[mask_1st] .+= transformations*q
    for (icstr,half) = enumerate(halves)
        mul!(joint_work, half, joint_q)
        ret[icstr] = dot(joint_q, joint_work)
    end
    mul!(trf_work, transformations, joint_q)
    for (i, idx) in enumerate(mask_1st)
        ret[idx] += trf_work[i]
    end
    ret .-= violations
end

function get_joint_jacobian!(ret, coords_hen::NC, coords_egg::NC,
        loci_position_hen, loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg
    )
    (;transformations,hessians,joint_q,joint_work,num_of_coords_hen,num_of_coords_egg) = cache
    assemble_joint_q!(joint_q, q_hen, q_egg, num_of_coords_hen, num_of_coords_egg)
    # for (icstr,hess) = enumerate(hessians)
    #     A = q'*hess
    #     ret[icstr,:] .+= A[1,:]
    # end
    # ret[mask_1st,:] .+= transformations[mask_1st,:]
    ncols = length(joint_work)
    for (icstr,hess) = enumerate(hessians)
        mul!(joint_work, hess, joint_q)
        for j = 1:ncols
            ret[icstr, j] += joint_work[j]
        end
    end
    rowvals_t = rowvals(transformations)
    nzvals_t = nonzeros(transformations)
    colptr_t = transformations.colptr
    for col = 1:length(colptr_t)-1
        for nz = colptr_t[col]:(colptr_t[col+1]-1)
            r = rowvals_t[nz]
            ret[mask_1st[r], col] += nzvals_t[nz]
        end
    end
end

function add_joint_forces_jacobian!(ret,num_of_cstr,
        coords_hen::NC, coords_egg::NC,
        loci_position_hen, loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg,
        λ
    )
    (;hessians) = cache
    for i = 1:num_of_cstr
        ret .+= -λ[i] .* hessians[i]
    end
end

function get_joint_velocity_jacobian!(ret,num_of_cstr,
        coords_hen::NC, coords_egg::NC,
        loci_position_hen, loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg,
        q̇_hen,q̇_egg,
    )
    (;hessians,joint_q,joint_work,num_of_coords_hen,num_of_coords_egg) = cache
    assemble_joint_q!(joint_q, q̇_hen, q̇_egg, num_of_coords_hen, num_of_coords_egg)
    q̇_jointed = joint_q
    # for i = 1:num_of_cstr
    #     ret[[i],:] .= transpose(q̇_jointed)*hessians[i]
    # end
    for i = 1:num_of_cstr
        mul!(joint_work, hessians[i], q̇_jointed)
        @views ret[i,:] .= joint_work
    end
end
