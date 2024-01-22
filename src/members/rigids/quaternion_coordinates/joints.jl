function build_joint_cache(
        nmcs_hen::QC, nmcs_egg::QC,
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
    quat_trl_rel_hen = QuatRotation(loci_axes_trl_hen.X).q |> vec
    quat_trl_rel_egg = QuatRotation(loci_axes_trl_egg.X).q |> vec
    # quat_rot_rel_hen = QuatRotation(loci_axes_rot_hen.X).q |> vec
    quat_rot_rel_egg = QuatRotation(loci_axes_rot_egg.X).q |> vec
    quat_hen = @view q_hen[4:7]
    quat_egg = @view q_egg[4:7]
    quat_trl_hen = Pmat(quat_hen)*quat_trl_rel_hen
    quat_trl_egg = Pmat(quat_egg)*quat_trl_rel_egg
    quat_rot_egg = Pmat(quat_egg)*quat_rot_rel_egg
    quat_rot_rel_hen = Pmat(Inv_mat*quat_hen)*quat_rot_egg
    quat_rot_hen = Pmat(quat_hen)*quat_rot_rel_hen
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
    rot_quat = Pmat(Inv_mat*quat_rot_egg)*quat_rot_hen
    vio_4th = (Pmat(Inv_mat*rot_quat)*rot_quat)[2:4] # always zero
    violations = vcat(
        vio_1st[mask_1st],
        vio_2nd[mask_2nd],
        vio_3rd[mask_3rd],
        vio_4th[mask_4th];
    )
    cache = @eponymtuple(
        quat_trl_rel_hen, quat_trl_rel_egg, quat_rot_rel_hen, quat_rot_rel_egg,
        rot_quat,
    )
    cache, violations
end

function get_joint_violations!(
        ret,
        nmcs_hen::QC, nmcs_egg::QC,
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
    quat_hen = @view q_hen[4:7]
    quat_egg = @view q_egg[4:7]
    quat_trl_hen = Pmat(quat_hen)*quat_trl_rel_hen
    quat_trl_egg = Pmat(quat_egg)*quat_trl_rel_egg
    quat_rot_hen = Pmat(quat_hen)*quat_rot_rel_hen
    quat_rot_egg = Pmat(quat_egg)*quat_rot_rel_egg
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
    vio_4th = (Pmat(Inv_mat*rot_quat)*(Pmat(Inv_mat*quat_rot_egg)*quat_rot_hen))[2:4]
    # @show vio_1st |> norm
    # @show violations[1:3]  |> norm
    ret .= vcat(
        vio_1st[mask_1st],
        vio_2nd[mask_2nd],
        vio_3rd[mask_3rd],
        vio_4th[mask_4th];
    ) .- violations
end

function make_cstr_jacobian(
        nmcs_hen::QC, nmcs_egg::QC,
        loci_position_hen,
        loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
    )
    num_of_coords_hen = get_num_of_coords(nmcs_hen)
    num_of_coords_egg = get_num_of_coords(nmcs_egg)
    num_of_jointed_coords = num_of_coords_hen + num_of_coords_egg
    c_hen = to_local_coords(nmcs_hen,loci_position_hen)
    c_egg = to_local_coords(nmcs_egg,loci_position_egg)
    (;quat_trl_rel_hen, quat_trl_rel_egg, quat_rot_rel_hen, quat_rot_rel_egg, rot_quat) = cache
    T = eltype(loci_position_hen)
    _jac_1st = zeros(T,3,num_of_jointed_coords)
    # hes 4th
    _jac_2nd = zeros(T,1,num_of_jointed_coords)
    # hes 3rd
    _jac_3rd = zeros(T,6,num_of_jointed_coords)
    # rotate
    _jac_4th = zeros(T,3,num_of_jointed_coords)
    jac_1st = DiffCache(_jac_1st)
    jac_2nd = DiffCache(_jac_2nd)
    jac_3rd = DiffCache(_jac_3rd)
    jac_4th = DiffCache(_jac_4th)
    _jac = vcat(
        _jac_1st[mask_1st,:],
        _jac_2nd[mask_2nd,:],
        _jac_3rd[mask_3rd,:],
        _jac_4th[mask_4th,:];
    )
    jac = DiffCache(_jac)
    function inner_cstr_jacobian(q)
        inner_jac = get_tmp(jac,q)
        inner_jac_1st = get_tmp(jac_1st,q)
        inner_jac_2nd = get_tmp(jac_2nd,q)
        inner_jac_3rd = get_tmp(jac_3rd,q)
        inner_jac_4th = get_tmp(jac_4th,q)
        inner_q_hen = @view q[1:num_of_coords_hen]
        inner_q_egg = @view q[num_of_coords_hen+1:end]
        quat_hen = @view inner_q_hen[4:7]
        quat_egg = @view inner_q_egg[4:7]
        quat_trl_hen = Pmat(quat_hen)*quat_trl_rel_hen
        quat_trl_egg = Pmat(quat_egg)*quat_trl_rel_egg
        quat_rot_hen = Pmat(quat_hen)*quat_rot_rel_hen
        quat_rot_egg = Pmat(quat_egg)*quat_rot_rel_egg
        # @show quat_trl_hen
        # @show quat_trl_egg
        # @show quat_rot_egg
        C_hen = to_position_jacobian(nmcs_hen,inner_q_hen,c_hen)
        C_egg = to_position_jacobian(nmcs_egg,inner_q_egg,c_egg)
        J = [-C_hen C_egg;] 
        r_hen = to_position(nmcs_hen,inner_q_hen,c_hen)
        r_egg = to_position(nmcs_egg,inner_q_egg,c_egg)
        d = r_egg - r_hen
        # hes 1st
        # jac 1st 
        inner_jac_1st[1:3,:] .= J
        # jac 4th 
        inner_jac_2nd[:,:]   .= 2d'*J
        # jac 3rd on hen
        # translate on hen
        inner_jac_3rd[1:3,:]                         .= Rmat(quat_trl_hen)'*J
        inner_jac_3rd[1:3,4:7]                      .+= ∂Rᵀf∂q(quat_trl_hen,d)*Mmat(quat_trl_rel_hen)
        # jac 3rd on egg
        # translate on egg
        inner_jac_3rd[4:6,:]                         .= Rmat(quat_trl_egg)'*J
        inner_jac_3rd[4:6,num_of_coords_hen.+(4:7)] .+= ∂Rᵀf∂q(quat_trl_egg,d)*Mmat(quat_trl_rel_egg)
        # jac 2nd
        # rotate of egg
        O43 = @SMatrix zeros(T,4,3)
        Pcq_egg = Pmat(Inv_mat*rot_quat)*(Inv_mat*quat_rot_rel_egg)
        inner_jac_4th[1:3,:] = (
            [
                O43 Pmat(Pcq_egg)*Pmat(Inv_mat*quat_egg)*Mmat(quat_rot_rel_hen) O43 Pmat(Pcq_egg)*Mmat(quat_rot_hen)*Inv_mat
            ]
        )[2:4,:]
        n1 = length(mask_1st)
        n2 = length(mask_2nd)
        n3 = length(mask_3rd)
        n4 = length(mask_4th)
        inner_jac[            (1:n1),:] = inner_jac_1st[mask_1st,:]
        inner_jac[        n1.+(1:n2),:] = inner_jac_2nd[mask_2nd,:]
        inner_jac[    n1.+n2.+(1:n3),:] = inner_jac_3rd[mask_3rd,:]
        inner_jac[n1.+n2.+n3.+(1:n4),:] = inner_jac_4th[mask_4th,:]
        inner_jac
    end
end

function get_joint_jacobian!(
        ret,
        nmcs_hen::QC, nmcs_egg::QC,
        loci_position_hen,
        loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg
    )
    q_jointed = vcat(
        q_hen,q_egg
    )
    A = make_cstr_jacobian(
        nmcs_hen, 
        nmcs_egg,
        loci_position_hen,
        loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
    )(q_jointed)
    ret .= A
end


function get_joint_forces_jacobian!(
        ret,
        num_of_cstr,
        nmcs_hen::QC, nmcs_egg::QC,
        loci_position_hen,
        loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg,
        λ
    )
    A = make_cstr_jacobian(
        nmcs_hen, 
        nmcs_egg,
        loci_position_hen,
        loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
    )
    function cstr_forces(q)
        jac = A(q)
        transpose(jac)*λ
    end
    q_jointed = vcat(
        q_hen,q_egg,
    )
    ret .= ForwardDiff.jacobian(cstr_forces, q_jointed)
end

function get_joint_velocity_jacobian!(
        ret,
        num_of_cstr,
        nmcs_hen::QC, nmcs_egg::QC,
        loci_position_hen,
        loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg,
        q̇_hen,q̇_egg,
    )
    A = make_cstr_jacobian(
        nmcs_hen, 
        nmcs_egg,
        loci_position_hen,
        loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
    )
    q_jointed = vcat(
        q_hen,q_egg,
    )
    q̇ = vcat(
        q̇_hen,q̇_egg,
    )
    function cstr_forces(q)
        A(q)*q̇
    end
    ret .= ForwardDiff.jacobian(cstr_forces, q_jointed)
end
