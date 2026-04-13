function clear_forces!(structure::AbstractStructure)
    structure.state.system.F .= 0
    clear_forces!(structure.bodies)
end
function clear_forces!(bodies::AbstractVector)
    foreach(clear_forces!,bodies)
end
function clear_forces!(bodies::TypeSortedCollection)
    foreach(clear_forces!,bodies)
end

"""
Update apparatuses of a structure
$(TYPEDSIGNATURES)
"""
function update_apparatuses!(structure::AbstractStructure, s::AbstractVector{T}) where {T}
    structure.state.system.s .= s
    update_apparatuses!(structure)
end

function update_apparatuses!(structure::AbstractStructure, s=nothing)
    (;apparatuses) = structure
    foreach(apparatuses) do appar
        appar_aux_var_idx = structure.connectivity.apparid2sys_aux_var_idx[appar.id]
        update_apparatus!(structure, appar, 
            (@view structure.state.system.s[appar_aux_var_idx])
        )
    end
end

function update_apparatus!(structure::AbstractStructure, appar::Apparatus, s)
    #default no update
end

function update_apparatus!(structure::AbstractStructure, appar::Apparatus{<:PrototypeJoint,<:TorsionalSpringDamper}, s)
    cnt = structure.connectivity
    (;system,members) = structure.state
    (;
      bodyid2sys_full_coords,
      num_of_full_coords,
      apparid2sys_full_coords_idx
    ) = cnt
    (;
        num_of_cstr,
        hen2egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = appar.joint
    (;hen,egg) = hen2egg
    full_idx = appar.full_coords_idx
    appar_sys_full_idx = apparid2sys_full_coords_idx[appar.id]
    spring_damper = appar.force
    (;state,k) = spring_damper
    state.angle = s[1]
    deformation = state.angle - state.rest_angle
    ## @show asin.(angles) .|> rad2deg
    state.torque = -k*deformation
end

function apply_force_to_bodies!(structure::AbstractStructure, appar::Apparatus{<:PrototypeJoint,<:TorsionalSpringDamper})
    cnt = structure.connectivity
    (;system,members) = structure.state
    (;
      bodyid2sys_full_coords,
      num_of_full_coords,
      apparid2sys_full_coords_idx
    ) = cnt
    (;
        num_of_cstr,
        hen2egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = appar.joint
    (;hen,egg) = hen2egg
    full_idx = appar.full_coords_idx
    appar_sys_full_idx = apparid2sys_full_coords_idx[appar.id]
    spring_damper = appar.force
    (;state,k) = spring_damper

    id_egg = egg.body.prop.id
    q_egg = @view system.q[bodyid2sys_full_coords[id_egg]]
    coords_egg = egg.body.coords
    X_egg = NCF.get_X(coords_egg,q_egg)
    invX̄_egg = coords_egg.data.invX̄
    ## angles_jacobian = ForwardDiff.jacobian(jointed2angles,q_jointed)
    ## @show q_jointed
    ## @show angles_jacobian
    # couple of forces
    loci_axes_rot_egg = egg.rot_axes
    axes_rot_egg = X_egg*invX̄_egg*loci_axes_rot_egg
    # X_hen*invX̄_hen*loci_axes_rot_hen == X_egg*invX̄_egg*loci_axes_rot_egg
    # R_hen*loci_axes_rot_hen == R_egg*loci_axes_rot_egg
    r̄_home = egg.position
    r̄_away = r̄_home + loci_axes_rot_egg.X[:,1]
    C_home = to_position_jacobian(coords_egg,q_egg,to_local_coords(coords_egg,r̄_home))
    C_away = to_position_jacobian(coords_egg,q_egg,to_local_coords(coords_egg,r̄_away))
    (;F) = members[egg.body.prop.id]
    F .+= state.torque*(-transpose(C_home)+transpose(C_away))*axes_rot_egg.X[:,2]
end

function update_apparatus!(structure::AbstractStructure, appar::Apparatus{<:BodyJoint,<:TorsionalSpringDamper}, s)
    (;joint,force) = appar
    spring_damper = appar.force
    (;state,k) = spring_damper
    
    state.angle = s[1]
    deformation = state.angle - state.rest_angle
    state.torque = -k*deformation
end

function apply_force_to_bodies!(structure::AbstractStructure, appar::Apparatus{<:BodyJoint,<:TorsionalSpringDamper})
   (;state,k) = appar.force

    body = appar.joint.body
    (;members) = structure.state
    
    (;F, q) = members[body.prop.id]
    
    num_of_dim = NCF.get_num_of_dims(body.coords)
    select_uvw = body.coords.conversion_to_X
    ∂x∂q = @view select_uvw[num_of_dim+1, :]
    ∂y∂q = @view select_uvw[num_of_dim+2, :]

    u,v = NCF.get_uv(body.coords,q)
    
    # F = torque * v (Geometric couple)
    # This avoids the singularity of implicit function theorem (1/Ss)
    
    # Map back to full coords using ∂x∂q and ∂y∂q

    F .+= state.torque .* (v[1] .* ∂x∂q .+ v[2] .* ∂y∂q)
end


function stretch!(structure::AbstractStructure,c)
    structure.state.system.c .= c
    stretch!(structure)
end

function stretch!(structure::AbstractStructure)
    (;bodies,apparatuses,state) = structure
    (;bodyid2sys_loci_coords_idx,apparid2params_idx) = structure.connectivity
    (;c) = state.system
    foreach(bodies) do body
        bodyid = body.prop.id
        @views stretch_loci!(body, c[bodyid2sys_loci_coords_idx[bodyid]])
    end
    foreach(apparatuses) do appar
        @views stretch!(appar, c[apparid2params_idx[appar.id]])
    end
end


function move_rigids!(structure)
    (;bodies,state) = structure
    foreach(bodies) do body
        bodyid = body.prop.id
        body_inst_state = state.members[bodyid]
        update_loci_states!(body,body_inst_state)
    end
end

function lazy_update_bodies!(structure::AbstractStructure,inst_state::AbstractCoordinatesState)
    structure.state.system.q .= inst_state.q
    structure.state.system.q̇ .= inst_state.q̇
    lazy_update_bodies!(structure)
end

function lazy_update_bodies!(structure::AbstractStructure)
    (;bodies,state) = structure
    foreach(bodies) do body
        bodyid = body.prop.id
        body_inst_state = state.members[bodyid]
        lazy_update_state!(body,body_inst_state)
        update_transformations!(body,body_inst_state)
        update_loci_states!(body,body_inst_state)
    end
    structure.cache.system.dirty = true
end


function update_bodies!(structure::AbstractStructure,inst_state::AbstractCoordinatesState)
    structure.state.system.q .= inst_state.q
    structure.state.system.q̇ .= inst_state.q̇
    update_bodies!(structure)
end

function update_bodies!(structure::AbstractStructure)
    (;bodies,state) = structure
    foreach(bodies) do body
        bodyid = body.prop.id
        body_inst_state = state.members[bodyid]
        update_inertia_cache!(body,body_inst_state)
        update_state!(body,body_inst_state)
        update_transformations!(body,body_inst_state)
        update_loci_states!(body,body_inst_state)
    end
    structure.cache.system.dirty = true
end


function apply_force_to_bodies!(structure::AbstractStructure, appar::Apparatus)
    # default no apply
end

"""
$(TYPEDSIGNATURES)
"""
function assemble_forces!(structure::AbstractStructure)
    (;bodies,apparatuses) = structure
    (;system,members) = structure.state

    foreach(apparatuses) do appar
        apply_force_to_bodies!(structure,appar)
    end

    foreach(bodies) do body
        (;F) = members[body.prop.id]
        (;prop,coords,state,cache) = body
        centrifugal_force!(F,state,cache)
        mass_center_force!(F,state,cache)
        concentrated_force!(F,state,cache)
        strain!(F,prop,coords,state,cache)
    end
end

function assemble_forces!(inst_state::AbstractCoordinatesState, structure::AbstractStructure)
    assemble_forces!(structure)
    inst_state.F .= structure.state.system.F
end

"""
$(TYPEDSIGNATURES)
"""
function apply_field!(structure::AbstractStructure, field::AbstractField ;)
    (;bodies,state) = structure
    foreach(bodies) do body
        (;id) = body.prop
        state.members[id].q
        apply_field!(
            state.members[id].F,
            body,
            field,
            state.members[id].q;
        )
    end
end


"""
$(TYPEDSIGNATURES)
"""
function field_jacobian!(∂F∂q, structure::AbstractStructure, field::AbstractField;)
    (;bodies,state,connectivity) = structure
    (;
        bodyid2sys_full_coords_idx,
    ) = connectivity
    foreach(bodies) do body
        (;id) = body.prop
        state.members[id].q
        field_jacobian!(
            (@view ∂F∂q[bodyid2sys_full_coords_idx[id],bodyid2sys_full_coords_idx[id]]),
            body,
            field,
            state.members[id].q;
        )
    end
end

