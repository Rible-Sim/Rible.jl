function clear_forces!(st::AbstractStructure)
    st.state.system.F .= 0
    clear_forces!(st.bodies)
end
function clear_forces!(bodies::AbstractVector)
    foreach(clear_forces!,bodies)
end
function clear_forces!(bodies::TypeSortedCollection)
    foreach(clear_forces!,bodies)
end


function rotation2angles(R::AbstractVector{T}) where {T}
    rotation2angles(reshape(R,3,:))
end

function rotation2angles(R::AbstractMatrix{T}) where {T}
    c1 = [1.0,0,0]×(R[:,1]/norm(R[:,1]))
    c2 = [0,1.0,0]×(R[:,2]/norm(R[:,2]))
    c3 = [0,0,1.0]×(R[:,3]/norm(R[:,3]))
    asin.( [0,0,c3[1]])   
end

function make_jointed2angles(hen2egg,relative_core)
    (;hen,egg) = hen2egg
    nmcs_hen = hen.body.coords.nmcs
    nmcs_egg = egg.body.coords.nmcs
    num_of_coords_hen = get_num_of_coords(nmcs_hen)
    num_of_coords_egg = get_num_of_coords(nmcs_egg)
    num_of_jointed_coords = num_of_coords_hen + num_of_coords_egg
    function inner_jointed2angles(q_jointed)
        q_hen = @view q_jointed[1:num_of_coords_hen]
        q_egg = @view q_jointed[num_of_coords_hen+1:end]
        X_hen = NCF.get_X(nmcs_hen,q_hen)
        X_egg = NCF.get_X(nmcs_egg,q_egg)
        relative_axes = X_egg*X_hen'
        angles = rotation2angles(relative_axes) 
        angles
    end
end

"""
Update DistanceSpringDamper Tension 
$(TYPEDSIGNATURES)
"""
function update_apparatuses!(st::AbstractStructure)
    (;apparatuses) = st
    foreach(apparatuses) do appar
        update_apparatus!(st, appar)
    end
end

function update_apparatus!(st::AbstractStructure, appar::Apparatus)
    #default no update
end

function update_apparatus!(st::AbstractStructure, appar::Apparatus{<:CableJoint})
    (;id,joint,force) = appar
    (;hen,egg) = joint.hen2egg
    locus_state_hen = hen.body.state.loci_states[hen.pid]
    locus_state_egg = egg.body.state.loci_states[egg.pid]
    p_hen = locus_state_hen.frame.position
    ṗ_hen = locus_state_hen.frame.velocity
    f_hen = locus_state_hen.force
    p_egg = locus_state_egg.frame.position
    ṗ_egg = locus_state_egg.frame.velocity
    f_egg = locus_state_egg.force
    update!(force,p_hen,p_egg,ṗ_hen,ṗ_egg)
    f_hen .+= force.state.force
    f_egg .-= force.state.force
end

function update_apparatus!(st::AbstractStructure, appar::Apparatus{<:ClusterJoint})
    (;id, joint, force) = appar
    (;sps) = joint
    nAppar = length(appar.force)
    for (idx, iappar) in enumerate(force)
        ijoint = iappar.joint
        iforce = iappar.force
        hen = ijoint.hen2egg.hen
        egg = ijoint.hen2egg.egg
        locus_state_hen = hen.body.state.loci_states[hen.pid]
        locus_state_egg = egg.body.state.loci_states[egg.pid]
        p_hen = locus_state_hen.frame.position
        ṗ_hen = locus_state_hen.frame.velocity
        f_hen = locus_state_hen.force
        p_egg = locus_state_egg.frame.position
        ṗ_egg = locus_state_egg.frame.velocity
        f_egg = locus_state_egg.force
        if (idx==1)
            s1 = 0
            s2 = sps[idx].s
        elseif (idx==nAppar)
            s1 = sps[nAppar-1].s
            s2 = 0
        else
            s1 = sps[idx-1].s
            s2 = sps[idx].s
        end
        update!(iforce, p_hen, p_egg, ṗ_hen, ṗ_egg, s1, s2)
        f_hen .+= iforce.state.force
        f_egg .-= iforce.state.force
    end
end

function update_apparatus!(st::AbstractStructure, appar::Apparatus{<:PrototypeJoint,<:RotationalSpringDamper})
    (;indexed,numbered) = st.connectivity
    (;system) = st.state
    (;bodyid2sys_free_coords,
      bodyid2sys_full_coords,
      num_of_free_coords,
      num_of_full_coords,
      apparid2sys_free_coords_idx
    ) = indexed
    (;
        num_of_cstr,
        hen2egg,
        cache,
        mask_1st,mask_2nd,mask_3rd,mask_4th
    ) = appar.joint
    (;
        relative_core
    ) = cache
    full_idx = appar.full_coords_idx
    free_idx = appar.free_coords_idx
    sys_free_coords_idx = apparid2sys_free_coords_idx[appar.id]
    spring_damper = appar.force
    (;mask,k) = spring_damper
    (;hen,egg) = hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    q = get_coords(st)
    q_hen = @view q[bodyid2sys_full_coords[id_hen]]
    q_egg = @view q[bodyid2sys_full_coords[id_egg]]
    q_jointed = vcat(
        q_hen,
        q_egg
    )
    jointed2angles = make_jointed2angles(hen2egg,relative_core)
    angles = jointed2angles(q_jointed)
    ## @show asin.(angles) .|> rad2deg
    torques = k.*angles
    update!(spring_damper,angles,)
    angles_jacobian = ForwardDiff.jacobian(jointed2angles,q_jointed)
    ## @show q_jointed
    ## @show angles_jacobian
    for i in mask
        torque = torques[i]
        system.F̌[sys_free_coords_idx] .-= angles_jacobian[i,free_idx]*torque
    end
end

function update_s!(st::AbstractStructure, s̄)
    (;apparatuses) = st
    st.state.system.s .= s̄
    foreach(apparatuses) do appar 
        if isa(appar, Apparatus{<:ClusterJoint})
            (;id) = appar
            (;sps) = appar.joint
            idx = st.connectivity.indexed.apparid2sys_add_var_idx[id]
            for (i, sp) in enumerate(sps)
                sp.s⁺ = s̄[idx[2i-1]]
                sp.s⁻ = s̄[idx[2i]]
                sp.s = sp.s⁺ - sp.s⁻
            end
        end
    end
end

function stretch_rigids!(st::AbstractStructure,c)
    st.state.system.c .= c
    stretch_rigids!(st)
end

function stretch_rigids!(st)
    (;bodies,state) = st
    (;bodyid2sys_loci_coords_idx) = st.connectivity.numbered
    (;c) = state.system
    foreach(bodies) do body
        bodyid = body.prop.id
        stretch_loci!(body,c[bodyid2sys_loci_coords_idx[bodyid]])
    end
end

function move_rigids!(st::AbstractStructure,q,q̇=zero(q))
    st.state.system.q .= q
    st.state.system.q̇ .= q̇
    move_rigids!(st)
end

function move_rigids!(st)
    (;bodies,state) = st
    foreach(bodies) do body
        bodyid = body.prop.id
        (;q,q̇) = state.parts[bodyid]
        update_loci_states!(body,q,q̇)
    end
end

function lazy_update_bodies!(st::AbstractStructure,q)
    st.state.system.q .= q
    st.state.system.q̇ .= 0.0
    lazy_update_bodies!(st)
end

function lazy_update_bodies!(st::AbstractStructure,q,q̇)
    st.state.system.q .= q
    st.state.system.q̇ .= q̇
    lazy_update_bodies!(st)
end

function lazy_update_bodies!(st::AbstractStructure)
    (;bodies,state) = st
    foreach(bodies) do body
        bodyid = body.prop.id
        (;q, q̇) = state.members[bodyid]
        lazy_update_state!(body,q,q̇)
        update_transformations!(body,q)
        update_loci_states!(body,q,q̇)
    end
end

function update_bodies!(st::AbstractStructure,q)
    st.state.system.q .= q
    st.state.system.q̇ .= 0.0
    update_bodies!(st)
end

function update_bodies!(st::AbstractStructure,q,q̇)
    st.state.system.q .= q
    st.state.system.q̇ .= q̇
    update_bodies!(st)
end

function update_bodies!(st::AbstractStructure)
    (;bodies,state) = st
    foreach(bodies) do body
        bodyid = body.prop.id
        (;q, q̇) = state.members[bodyid]
        update_inertia_cache!(body,q,q̇)
        update_state!(body,q,q̇)
        update_transformations!(body,q)
        update_loci_states!(body,q,q̇)
    end
end

"""
$(TYPEDSIGNATURES)
"""
function assemble_forces!(generliazed_force,st::Structure)
    generliazed_force .= assemble_forces!(st::Structure)
end

"""
$(TYPEDSIGNATURES)
"""
function assemble_forces!(st::Structure)
    (;bodies,state) = st
    (;system,members) = state
    foreach(bodies) do body
        (;F) = members[body.prop.id]
        (;state,cache) = body
        centrifugal_force!(F,state,cache)
        mass_center_force!(F,state,cache)
        concentrated_force!(F,state,cache)
        strain!(F,state,cache)
    end
    system.F̌
end

function get_force(st::AbstractStructure)
    st.state.system.F̌
end

function get_force!(F,st::AbstractStructure)
    F .= get_force(st)
end

"""
$(TYPEDSIGNATURES)
"""
function apply_gravity!(st;factor=1)
    (;bodies) = st
    gravity_acceleration = factor*get_gravity(st)
    foreach(bodies) do body
        body.state.mass_locus_state.force .+= gravity_acceleration*body.prop.mass
    end
end