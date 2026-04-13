

# out of space for convinience only
function cstr_function(structure::AbstractStructure, inst_state::AbstractCoordinatesState)
    ret = zeros(eltype(inst_state.q),structure.connectivity.num_of_cstr)
    cstr_function!(ret, structure, inst_state)
    ret
end

function cstr_function!(ret, structure::AbstractStructure, inst_state::AbstractCoordinatesState)
    (;bodies,apparatuses,connectivity) = structure
    cnt = connectivity
    (;  
        num_of_cstr,
        num_of_full_coords,
        bodyid2sys_full_coords,
        num_of_intrinsic_cstr,
        bodyid2sys_intrinsic_cstr_idx,
        apparid2sys_extrinsic_cstr_idx
    ) = cnt
    (;q,) = inst_state
    fill!(ret, zero(eltype(q)))
    foreach(bodies) do body
        bodyid = body.prop.id
        memfull = bodyid2sys_full_coords[bodyid]
        memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
        cstr_function!((@view ret[memincst]), 
            body.coords, body.cache,
            (@view q[memfull])#,d[memincst]
        )
    end
    foreach(apparatuses) do appar
        jointexcst = apparid2sys_extrinsic_cstr_idx[appar.id]
        cstr_function!((@view ret[jointexcst]), appar,cnt,
            inst_state
        )
    end
end


# out of space for convinience only
function cstr_jacobian(structure::AbstractStructure, inst_state::AbstractCoordinatesState)
    (;num_of_cstr,num_of_full_coords) = structure.connectivity
    (;q,) = inst_state
    ret = zeros(eltype(q),num_of_cstr,num_of_full_coords)
    cstr_jacobian!(ret,structure,inst_state)
    ret
end

function cstr_jacobian!(ret, structure::AbstractStructure, inst_state::AbstractCoordinatesState)
    (;bodies,apparatuses,connectivity) = structure
    (;
        num_of_cstr,
        num_of_full_coords,
        bodyid2sys_full_coords,
        num_of_intrinsic_cstr,
        bodyid2sys_intrinsic_cstr_idx,
        apparid2sys_full_coords_idx,
        apparid2sys_extrinsic_cstr_idx
    ) = connectivity
    (;q,) = inst_state
    structure.state.system.q .= inst_state.q
    fill!(ret, zero(eltype(q)))
    foreach(bodies) do body
        bodyid = body.prop.id
        memfull = bodyid2sys_full_coords[bodyid]
        memincst = bodyid2sys_intrinsic_cstr_idx[bodyid]
        cstr_jacobian!(
            (@view ret[memincst,memfull]),
            body,structure.state.members[bodyid]
        )
    end
    foreach(apparatuses) do appar
        (;id,) = appar
        ext = apparid2sys_extrinsic_cstr_idx[id]
        cstr_jacobian!((@view ret[ext,apparid2sys_full_coords_idx[id]]), 
            appar,connectivity,inst_state
        )
    end
end

function cstr_function!(ret,appar::Apparatus,cnt,inst_state)
    # default: no constraints
end

function cstr_jacobian!(ret, appar::Apparatus,cnt,inst_state)
    # default: no constraints
end

function add_cstr_forces_jacobian!(ret, appar::Apparatus,cnt,q,λ)
    # default: no constraints
end

function cstr_velocity_jacobian!(ret, appar::Apparatus,cnt,q,q̇)
    # default: no constraints
end

function cstr_function!(ret, appar::Apparatus{<:LinearJoint},cnt,inst_state)
    appar_sys_full_coords_idx = cnt.apparid2sys_full_coords_idx[appar.id]
    (;A,violations) = appar.joint
    q_view = @view inst_state.q[appar_sys_full_coords_idx]
    mul!(ret, A, q_view)
    @. ret -= violations
end

function cstr_jacobian!(ret, appar::Apparatus{<:LinearJoint},cnt,inst_state)
    ret .+= appar.joint.A
end

function add_cstr_forces_jacobian!(ret, appar::Apparatus{<:LinearJoint},cnt,q,λ)
    # linear constraints vanish here
end

function cstr_velocity_jacobian!(ret, appar::Apparatus{<:LinearJoint},cnt,q,q̇)
    # linear constraints vanish here
end

function cstr_function!(ret, appar::Apparatus{<:PrototypeJoint},cnt,inst_state)
    (;bodyid2sys_full_coords) = cnt
    (;sys_locus_id2coords_idx,bodyid2sys_locus_id) = cnt
    (;
        num_of_cstr,hen2egg,
        violations,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        loci_position_hen,loci_position_egg
    ) = appar.joint
    (;hen,egg) = hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    coords_hen = hen.body.coords
    coords_egg = egg.body.coords
    (;q) = inst_state
    q_hen = @view q[bodyid2sys_full_coords[id_hen]]
    q_egg = @view q[bodyid2sys_full_coords[id_egg]]
    # cstr violations
    # @show joint.id
    get_joint_violations!(
        ret,
        coords_hen, coords_egg,
        loci_position_hen,
        loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg,
        violations
    )
end

function cstr_jacobian!(ret, appar::Apparatus{<:PrototypeJoint},cnt,inst_state)
    (;bodyid2sys_full_coords,
    ) = cnt
    (;sys_locus_id2coords_idx,bodyid2sys_locus_id) = cnt
    (;
        num_of_cstr,
        hen2egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        loci_position_hen,loci_position_egg
    ) = appar.joint
    (;q) = inst_state
    full_idx = appar.full_coords_idx
    (;hen,egg) = hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    coords_hen = hen.body.coords
    coords_egg = egg.body.coords
    q_hen = @view q[bodyid2sys_full_coords[id_hen]]
    q_egg = @view q[bodyid2sys_full_coords[id_egg]]
    # translate
    get_joint_jacobian!(
        ret,
        coords_hen, coords_egg,
        loci_position_hen,
        loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg
    )
end

function add_cstr_forces_jacobian!(ret, appar::Apparatus{<:PrototypeJoint},cnt,q,λ)
    (;bodyid2sys_full_coords,
    ) = cnt
    (;
        num_of_cstr,
        hen2egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        loci_position_hen,loci_position_egg
    ) = appar.joint
    full_idx = appar.full_coords_idx
    (;hen,egg) = hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    q_hen = @view q[bodyid2sys_full_coords[id_hen]]
    q_egg = @view q[bodyid2sys_full_coords[id_egg]]
    coords_hen = hen.body.coords
    coords_egg = egg.body.coords
    add_joint_forces_jacobian!(
        ret,
        num_of_cstr,
        coords_hen, coords_egg,
        loci_position_hen,loci_position_egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg,
        λ
    )
end

function cstr_velocity_jacobian!(ret, appar::Apparatus{<:PrototypeJoint},cnt,q,q̇)
    (;bodyid2sys_full_coords,
    ) = cnt
    (;
        num_of_cstr,
        hen2egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = appar.joint
    (;hen,egg) = hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    q_hen = @view q[bodyid2sys_full_coords[id_hen]]
    q_egg = @view q[bodyid2sys_full_coords[id_egg]]
    q̇_hen = @view q̇[bodyid2sys_full_coords[id_hen]]
    q̇_egg = @view q̇[bodyid2sys_full_coords[id_egg]]
    coords_hen = hen.body.coords
    coords_egg = egg.body.coords
    fill!(ret, zero(eltype(q̇)))
    get_joint_velocity_jacobian!(
        ret,
        num_of_cstr,
        coords_hen, coords_egg,
        hen.position,
        egg.position,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th,
        q_hen,q_egg,
        q̇_hen,q̇_egg,
    )
end

# FixedPointJoint Constraints
function cstr_function!(ret, appar::Apparatus{<:FixedPointJoint},cnt,inst_state)
    (;bodyid2sys_full_coords,
    ) = cnt
    (;sys_locus_id2coords_idx,bodyid2sys_locus_id) = cnt
    (;
        body,
        num_of_cstr,
        r̄o, ro
    ) = appar.joint
    full_idx = appar.full_coords_idx
    (;prop,coords) = body
    (;id) = prop
    c = to_local_coords(coords,r̄o)
    (;q) = inst_state
    q = @view q[bodyid2sys_full_coords[id]]
    ret .= to_position(coords,q,c) .- ro
end

function cstr_jacobian!(ret, appar::Apparatus{<:FixedPointJoint},cnt,inst_state) 
    (;bodyid2sys_full_coords,
    ) = cnt
    (;sys_locus_id2coords_idx,bodyid2sys_locus_id) = cnt
    (;
        body,
        num_of_cstr,
        r̄o,
    ) = appar.joint
    full_idx = appar.full_coords_idx
    (;prop,coords) = body
    (;id) = prop
    (;q) = inst_state
    q = @view q[bodyid2sys_full_coords[id]]
    c = to_local_coords(coords,r̄o)
    ret .+= to_position_jacobian(coords,q,c)
end

function add_cstr_forces_jacobian!(ret, appar::Apparatus{<:FixedPointJoint},cnt,q,λ) 
    (;bodyid2sys_full_coords,
    ) = cnt
    (;
        num_of_cstr,
        body,
    ) = appar.joint
    full_idx = appar.full_coords_idx
    id = body.prop.id
    q = @view q[bodyid2sys_full_coords[id]]
    T = eltype(λ)
    num_of_full_idx = length(full_idx)
    #todo 
    # add_joint_forces_jacobian!(
    #     ret,
    #     num_of_cstr,
    #     nmcs, 
    #     body.prop.loci[pid].position,
    #     cache,
    #     mask_1st,mask_2nd,mask_3rd,mask_4th,
    #     q,
    #     λ
    # )
end

function cstr_velocity_jacobian!(ret, appar::Apparatus{<:FixedPointJoint},cnt,q,q̇) 
    
end
