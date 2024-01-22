
function cstr_function(appar::Apparatus{<:CableJoint},structure,q,c)
    T = get_numbertype(structure)
    zeros(T,0)
end

function cstr_jacobian(appar::Apparatus{<:CableJoint},structure,q)
    (;apparid2sys_free_coords_idx) = structure.connectivity.indexed
    jointed_sys_free_idx = apparid2sys_free_coords_idx[appar.id]
    T = get_numbertype(structure)
    zeros(T,0,length(jointed_sys_free_idx))
end

function cstr_forces_jacobian(appar::Apparatus{<:CableJoint},structure,q,λ)
    (;apparid2sys_free_coords_idx) = structure.connectivity.indexed
    jointed_sys_free_idx = apparid2sys_free_coords_idx[appar.id]
    T = get_numbertype(structure)
    zeros(T,length(jointed_sys_free_idx),length(jointed_sys_free_idx))

end

function cstr_velocity_jacobian(appar::Apparatus{<:CableJoint},structure,q,q̇)
    (;apparid2sys_free_coords_idx) = structure.connectivity.indexed
    jointed_sys_free_idx = apparid2sys_free_coords_idx[appar.id]
    T = get_numbertype(structure)
    zeros(T,0,length(jointed_sys_free_idx))
end

function cstr_function(appar::Apparatus{<:LinearJoint},structure,q,c)
    (;indexed) = structure.connectivity
    sys_free_coords_idx = indexed.apparid2sys_free_coords_idx[appar.id]
    (;A,violations) = appar.joint
    A*q[sys_free_coords_idx] .- violations
end

function cstr_jacobian(appar::Apparatus{<:LinearJoint},structure,q)
    (;indexed) = structure.connectivity
    appar.joint.A
end

function cstr_forces_jacobian(appar::Apparatus{<:LinearJoint},structure,q,λ)
    (;indexed) = structure.connectivity
    sys_free_coords_idx = indexed.apparid2sys_free_coords_idx[appar.id]
    n = length(sys_free_coords_idx)
    zeros(eltype(λ),n,n)
end

function cstr_velocity_jacobian(appar::Apparatus{<:LinearJoint},structure,q,q̇)
    (;indexed) = structure.connectivity
    sys_free_coords_idx = indexed.apparid2sys_free_coords_idx[appar.id]
    (;num_of_cstr,) = appar.joint
    n = length(sys_free_coords_idx)
    zeros(eltype(q̇),num_of_cstr,n)
end

function cstr_function(appar::Apparatus{<:PrototypeJoint},structure::Structure,q, c = get_local_coords(structure))
    (;indexed,numbered) = structure.connectivity
    (;bodyid2sys_full_coords) = indexed
    (;sys_loci2coords_idx,bodyid2sys_loci_idx) = numbered
    (;
        num_of_cstr,hen2egg,
        violations,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = appar.joint
    (;hen,egg) = hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    nmcs_hen = hen.body.coords.nmcs
    nmcs_egg = egg.body.coords.nmcs
    T = eltype(q)
    ret = zeros(T,num_of_cstr)
    q_hen = @view q[bodyid2sys_full_coords[id_hen]]
    q_egg = @view q[bodyid2sys_full_coords[id_egg]]
    q = vcat(q_hen,q_egg)
    # cstr violations
    # @show joint.id
    get_joint_violations!(
        ret,
        nmcs_hen, nmcs_egg,
        hen.body.prop.loci[hen.pid].position,
        egg.body.prop.loci[egg.pid].position,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg,
        violations
    )
    ret
end

function cstr_jacobian(appar::Apparatus{<:PrototypeJoint},structure::Structure,q,c = get_local_coords(structure))
    (;indexed,numbered) = structure.connectivity
    (;bodyid2sys_free_coords,
      bodyid2sys_full_coords,
      num_of_free_coords,
      num_of_full_coords,
      apparid2full_idx,
      apparid2free_idx,
    ) = indexed
    (;sys_loci2coords_idx,bodyid2sys_loci_idx) = numbered
    (;
        num_of_cstr,
        hen2egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = appar.joint
    full_idx = apparid2full_idx[appar.id]
    free_idx = apparid2free_idx[appar.id]
    (;hen,egg) = hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    nmcs_hen = hen.body.coords.nmcs
    nmcs_egg = egg.body.coords.nmcs
    T = eltype(q)
    ret = zeros(T,num_of_cstr,length(full_idx))
    q_hen = @view q[bodyid2sys_full_coords[id_hen]]
    q_egg = @view q[bodyid2sys_full_coords[id_egg]]
    # translate
    get_joint_jacobian!(
        ret,
        nmcs_hen, nmcs_egg,
        hen.body.prop.loci[hen.pid].position,
        egg.body.prop.loci[egg.pid].position,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg
    )
    @view ret[:,free_idx]
end

function cstr_forces_jacobian(appar::Apparatus{<:PrototypeJoint},structure,q,λ)
    (;indexed,numbered) = structure.connectivity
    (;bodyid2sys_free_coords,
      bodyid2sys_full_coords,
      num_of_free_coords,
      num_of_full_coords,
      apparid2full_idx,
      apparid2free_idx
    ) = indexed
    (;
        num_of_cstr,
        hen2egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = appar.joint
    full_idx = apparid2full_idx[appar.id]
    free_idx = apparid2free_idx[appar.id]
    (;hen,egg) = hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    q_hen = @view q[bodyid2sys_full_coords[id_hen]]
    q_egg = @view q[bodyid2sys_full_coords[id_egg]]
    nmcs_hen = hen.body.coords.nmcs
    nmcs_egg = egg.body.coords.nmcs
    T = eltype(λ)
    num_of_full_idx = length(full_idx)
    ret = zeros(T,num_of_full_idx,num_of_full_idx)
    get_joint_forces_jacobian!(
        ret,
        num_of_cstr,
        nmcs_hen, nmcs_egg,
        hen.body.prop.loci[hen.pid].position,
        egg.body.prop.loci[egg.pid].position,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg,
        λ
    )
    @view ret[free_idx,free_idx]
end

function cstr_velocity_jacobian(appar::Apparatus{<:PrototypeJoint},structure,q,q̇)
    (;indexed,numbered) = structure.connectivity
    (;bodyid2sys_free_coords,
      bodyid2sys_full_coords,
      num_of_free_coords,
      num_of_full_coords,
      apparid2full_idx,
      apparid2free_idx
    ) = indexed
    (;
        num_of_cstr,
        hen2egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = appar.joint
    full_idx = apparid2full_idx[appar.id]
    free_idx = apparid2free_idx[appar.id]
    (;hen,egg) = hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    q_hen = @view q[bodyid2sys_full_coords[id_hen]]
    q_egg = @view q[bodyid2sys_full_coords[id_egg]]
    q̇_hen = @view q̇[bodyid2sys_full_coords[id_hen]]
    q̇_egg = @view q̇[bodyid2sys_full_coords[id_egg]]
    nmcs_hen = hen.body.coords.nmcs
    nmcs_egg = egg.body.coords.nmcs
    T = eltype(q̇)
    num_of_full_idx = length(full_idx)
    ret = zeros(T,num_of_cstr,num_of_full_idx)
    get_joint_velocity_jacobian!(
        ret,
        num_of_cstr,
        nmcs_hen, nmcs_egg,
        hen.body.prop.loci[hen.pid].position,
        egg.body.prop.loci[egg.pid].position,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg,
        q̇_hen,q̇_egg,
    )
    @view ret[:,free_idx]
end