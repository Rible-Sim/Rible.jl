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

"""
Update DistanceSpringDamper Tension 
$(TYPEDSIGNATURES)
"""
function update_tensiles!(st::AbstractStructure)
    update_tensiles!(st::AbstractStructure, st.connectivity.tensioned)
    update_tensiles!(st::AbstractStructure, st.connectivity.jointed)
end

function update_tensiles!(st::AbstractStructure, @eponymargs(connected,))
    (;cables) = st.tensiles
    foreach(connected) do scnt
        cable = cables[scnt.id]
        locus_state_hen = scnt.hen.bodysig.state.loci_states[scnt.hen.pid]
        locus_state_egg = scnt.egg.bodysig.state.loci_states[scnt.egg.pid]
        p_hen = locus_state_hen.frame.position
        ṗ_hen = locus_state_hen.frame.velocity
        f_hen = locus_state_hen.force
        p_egg = locus_state_egg.frame.position
        ṗ_egg = locus_state_egg.frame.velocity
        f_egg = locus_state_egg.force
        update!(cable,p_hen,p_egg,ṗ_hen,ṗ_egg)
        f_hen .+= cable.state.force
        f_egg .-= cable.state.force
    end
end

function rotation2angles(R::AbstractVector{T}) where {T}
    rotation2angles(reshape(R,3,:))
end

function rotation2angles(R::AbstractMatrix{T}) where {T}
    R31 = R[3,1]
    if R31 < 1
        if R31 > -1
            θy = asin(-R31)
            cosθy = cos(θy)
            θz = atan((R[2,1])/cosθy,(R[1,1])/cosθy)
            θx = atan((R[3,2])/cosθy,(R[3,3])/cosθy)
        else # R31 == -1
            # Not a unique solution
            θx = zero(T)
            θy = T(π/2)
            θz =  -atan(-R[2,3],R[2,2])
        end
    else # R31 == 1
        # Not a unique solution
        θx = zero(T)
        θy = -π/2
        θz = atan(-R[2,3],R[2,2])
    end
    [θz,θy,θx]
    c1 = [1.0,0,0]×(R[:,1]/norm(R[:,1]))
    c2 = [0,1.0,0]×(R[:,2]/norm(R[:,2]))
    c3 = [0,0,1.0]×(R[:,3]/norm(R[:,3]))
    asin.([0,0,c3[1]])
end

function make_jointed2angles(hen2egg,relative_core)
    (;hen,egg) = hen2egg
    nmcs_hen = hen.bodysig.coords.nmcs
    nmcs_egg = egg.bodysig.coords.nmcs
    num_of_coords_hen = get_num_of_coords(nmcs_hen)
    num_of_coords_egg = get_num_of_coords(nmcs_egg)
    num_of_jointed_coords = num_of_coords_hen + num_of_coords_egg
    function inner_jointed2angles(q_jointed)
        q_hen = @view q_jointed[1:num_of_coords_hen]
        q_egg = @view q_jointed[num_of_coords_hen+1:end]
        X_hen = NCF.get_X(nmcs_hen,q_hen)
        X_egg = NCF.get_X(nmcs_egg,q_egg)
        relative_axes = X_egg*X_hen'
        # angles = [0,0,acos(X_hen[:,2]'*X_egg[:,2]/(norm(X_hen[:,2])*norm(X_egg[:,2])))]
        angles = rotation2angles(relative_axes) 
        angles
    end
end

function update_tensiles!(st::AbstractStructure, jointed::Jointed)
    (;spring_dampers) = st.tensiles
    (;indexed,numbered) = st.connectivity
    (;system) = st.state
    (;bodyid2sys_free_coords,
      bodyid2sys_full_coords,
      num_of_free_coords,
      num_of_full_coords,
    ) = indexed
    q = get_coords(st)
    foreach(jointed.joints) do joint
        if joint isa PrototypeJoint
            (;
                num_of_cstr,
                hen2egg,
                cache,
                mask_1st,mask_2nd,mask_3rd,mask_4th
            ) = joint
            (;
                relative_core
            ) = cache
            _, free_idx, sys_free_idx = get_joint_idx(joint,indexed)
            spring_damper = spring_dampers[joint.id]
            (;mask,k) = spring_damper
            (;hen,egg) = hen2egg
            nmcs_hen = hen.bodysig.coords.nmcs
            nmcs_egg = egg.bodysig.coords.nmcs
            num_of_coords_hen = get_num_of_coords(nmcs_hen)
            num_of_coords_egg = get_num_of_coords(nmcs_egg)
            id_hen = hen.bodysig.prop.id
            id_egg = egg.bodysig.prop.id
            q_hen = @view q[bodyid2sys_full_coords[id_hen]]
            q_egg = @view q[bodyid2sys_full_coords[id_egg]]
            q_jointed = vcat(
                q_hen,
                q_egg
            )
            jointed2angles = make_jointed2angles(hen2egg,relative_core)
            angles = jointed2angles(q_jointed)
	        torques = k.*angles
            update!(spring_damper,angles,)
            angles_jacobian = ForwardDiff.jacobian(jointed2angles,q_jointed)
            for i in mask
                torque = torques[i]
                angle = angles[i]
                # @show angle
                # @show angles_jacobian[i,  4:12]*torque
                # @show angles_jacobian[i,16:end]*torque
                system.F̌[sys_free_idx] .-= angles_jacobian[i,free_idx]*torque
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
function assemble_forces!(st::Structure)
    (;bodies,state) = st
    (;system,members) = state
    foreach(bodies) do body
        (;F) = members[body.prop.id]
        (;state,cache) = body
        centrifugal_force!(F,state,cache)
        mass_center_force!(F,state,cache)
        concentrated_force!(F,state,cache)
        # strain!(F,state)
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