function build_joint_cache(
        nmcs_hen::QC,
        nmcs_egg::QC,
        loci_position_hen,
        loci_position_egg,
        loci_axes_trl_hen,
        loci_axes_trl_egg,
        loci_axes_rot_hen,
        loci_axes_rot_egg,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg
    )
    # translate
    c_hen = to_local_coords(nmcs_hen,loci_position_hen)
    c_egg = to_local_coords(nmcs_egg,loci_position_egg)
    quat_trl_rel_hen = QuatRotation(loci_axes_trl_hen.X[:,[2,3,1]]).q
    quat_trl_rel_egg = QuatRotation(loci_axes_trl_egg.X[:,[2,3,1]]).q
    # quat_rot_rel_hen = QuatRotation(loci_axes_rot_hen.X[:,[2,3,1]]).q
    quat_rot_rel_egg = QuatRotation(loci_axes_rot_egg.X[:,[2,3,1]]).q
    @show loci_axes_rot_egg.X
    quat_hen = Quaternion(q_hen[4:7]...)
    quat_egg = Quaternion(q_egg[4:7]...)
    quat_trl_hen = quat_hen*quat_trl_rel_hen
    quat_trl_egg = quat_egg*quat_trl_rel_egg
    quat_rot_egg = quat_egg*quat_rot_rel_egg
    quat_rot_rel_hen = conj(quat_hen)*quat_rot_egg
    quat_rot_hen = quat_hen*quat_rot_rel_hen
    r_hen = to_position(nmcs_hen,q_hen,c_hen)
    r_egg = to_position(nmcs_egg,q_egg,c_egg)
    d = r_egg - r_hen
    # cstr 1st
    vio_1st = d
    vio_2nd = [d'*d]
    vio_3rd = vcat(
        Rmat(quat_trl_hen)'*d,
        Rmat(quat_trl_egg)'*d
    )
    rot_quat = conj(quat_rot_egg)*quat_rot_hen
    vio_4th = vec(conj(rot_quat)*rot_quat)[2:4] # always zero
    violations = vcat(
        vio_1st[mask_1st],
        vio_2nd[mask_2nd],
        vio_3rd[mask_3rd],
        vio_4th[mask_4th];
    )

    @show r_egg, r_hen
    cache = @eponymtuple(
        quat_trl_rel_hen, quat_trl_rel_egg, quat_rot_rel_hen, quat_rot_rel_egg,
        rot_quat,
    )
    cache, violations
end

function get_joint_violations!(
        ret,
        nmcs_hen::QC, 
        nmcs_egg::QC,
        loci_position_hen,
        loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg,
        violations
    )
    c_hen = to_local_coords(nmcs_hen,loci_position_hen)
    c_egg = to_local_coords(nmcs_egg,loci_position_egg)
    (;quat_trl_rel_hen, quat_trl_rel_egg, quat_rot_rel_hen, quat_rot_rel_egg,rot_quat) = cache
    quat_hen = Quaternion(q_hen[4:7]...)
    quat_egg = Quaternion(q_egg[4:7]...)
    quat_trl_hen = quat_hen*quat_trl_rel_hen
    quat_trl_egg = quat_egg*quat_trl_rel_egg
    quat_rot_hen = quat_hen*quat_rot_rel_hen
    quat_rot_egg = quat_egg*quat_rot_rel_egg
    r_hen = to_position(nmcs_hen,q_hen,c_hen)
    r_egg = to_position(nmcs_egg,q_egg,c_egg)
    d = r_egg - r_hen
    # cstr 1st
    vio_1st = d
    vio_2nd = [d'*d]
    vio_3rd = vcat(
        Rmat(quat_trl_hen)'*d,
        Rmat(quat_trl_egg)'*d
    )
    vio_4th = vec(conj(rot_quat)*conj(quat_rot_egg)*quat_rot_hen)[2:4]
    # @show vio_1st |> norm
    # @show violations[1:3]  |> norm
    ret .= vcat(
        vio_1st[mask_1st],
        vio_2nd[mask_2nd],
        vio_3rd[mask_3rd],
        vio_4th[mask_4th];
    ) .- violations
end

function get_joint_jacobian!(
        ret,
        nmcs_hen::QC, 
        nmcs_egg::QC,
        loci_position_hen,
        loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        free_idx_hen,free_idx_egg,
        free_hen,free_egg,
        q_hen,q_egg
    )

    num_of_coords_hen = get_num_of_coords(nmcs_hen)
    num_of_coords_egg = get_num_of_coords(nmcs_egg)
    num_of_jointed_coords = num_of_coords_hen + num_of_coords_egg
    c_hen = to_local_coords(nmcs_hen,loci_position_hen)
    c_egg = to_local_coords(nmcs_egg,loci_position_egg)
    (;quat_trl_rel_hen, quat_trl_rel_egg, quat_rot_rel_hen, quat_rot_rel_egg, rot_quat) = cache
    quat_hen = Quaternion(q_hen[4:7]...)
    quat_egg = Quaternion(q_egg[4:7]...)
    quat_trl_hen = quat_hen*quat_trl_rel_hen
    quat_trl_egg = quat_egg*quat_trl_rel_egg
    quat_rot_hen = quat_hen*quat_rot_rel_hen
    quat_rot_egg = quat_egg*quat_rot_rel_egg
    T = eltype(q_hen)
    # @show quat_trl_hen
    # @show quat_trl_egg
    # @show quat_rot_egg
    C_hen = to_transformation(nmcs_hen,q_hen,c_hen)
    C_egg = to_transformation(nmcs_egg,q_egg,c_egg)
    J = [-C_hen C_egg;] 
    r_hen = to_position(nmcs_hen,q_hen,c_hen)
    r_egg = to_position(nmcs_egg,q_egg,c_egg)
    d = r_egg - r_hen
    # hes 1st
    jac_1st = zeros(T,3,num_of_jointed_coords)
    # hes 4th
    jac_2nd = zeros(T,1,num_of_jointed_coords)
    # hes 3rd
    jac_3rd = zeros(T,6,num_of_jointed_coords)
    # rotate
    jac_4th = zeros(T,3,num_of_jointed_coords)
    # jac 1st 
    jac_1st[1:3,:] .= J
    # jac 4th 
    jac_2nd[:,:]   .= 2d'*J
    # jac 3rd on hen
    # translate on hen
    jac_3rd[1:3,:]                         .= Rmat(quat_trl_hen)'*J
    jac_3rd[1:3,4:7]                      .+= ∂Rᵀf∂q(quat_trl_hen,d)*Mmat(quat_trl_rel_hen)
    # jac 3rd on egg
    # translate on egg
    jac_3rd[4:6,:]                         .= Rmat(quat_trl_egg)'*J
    jac_3rd[4:6,num_of_coords_hen.+(4:7)] .+= ∂Rᵀf∂q(quat_trl_egg,d)*Mmat(quat_trl_rel_egg)
    # jac 2nd
    # rotate of egg
    O43 = @SMatrix zeros(T,4,3)
    Pcq_egg = Pmat(conj(rot_quat)*conj(quat_rot_rel_egg))
    jac_4th[1:3,:] = (
        [
            O43 Pcq_egg*Pmat(conj(quat_egg))*Mmat(quat_rot_rel_hen) O43 Pcq_egg*Mmat(quat_rot_hen)*Inv_mat
        ]
    )[2:4,:]
    jac = vcat(
        jac_1st[mask_1st,:],
        jac_2nd[mask_2nd,:],
        jac_3rd[mask_3rd,:],
        jac_4th[mask_4th,:];
    )
    # @show jac
    ret[:,free_hen] .= jac[:,free_idx_hen]
    ret[:,free_egg] .= jac[:,(free_idx_egg.+num_of_coords_hen)]
end