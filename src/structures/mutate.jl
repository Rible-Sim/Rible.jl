
"""
Clear Rigid Body 所受作用力和力矩。
$(TYPEDSIGNATURES)
"""
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
Update Cable Tension 
$(TYPEDSIGNATURES)
"""
function update_tensiles!(st::AbstractStructure)
    update_tensiles!(st, st.connectivity.tensioned)
end

function update_tensiles!(st::Structure, @eponymargs(connected,))
    (;cables) = st.tensiles
    foreach(connected) do scnt
        scable = cables[scnt.id]
        locus_state_hen = scnt.hen.rbsig.state.loci_states[scnt.hen.pid]
        locus_state_egg = scnt.egg.rbsig.state.loci_states[scnt.egg.pid]
        p_hen = locus_state_hen.position
        ṗ_hen = locus_state_hen.velocity
        f_hen = locus_state_hen.force
        p_egg = locus_state_egg.position
        ṗ_egg = locus_state_egg.velocity
        f_egg = locus_state_egg.force
        update!(scable,p_hen,p_egg,ṗ_hen,ṗ_egg)
        f_hen .+= scable.state.force
        f_egg .-= scable.state.force
    end
end

function stretch_rigids!(st,c)
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

function move_rigids!(st,q,q̇=zero(q))
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

"""
Update Rigid Body State.
$(TYPEDSIGNATURES)
"""
function lazy_update_bodies!(st,q)
    st.state.system.q .= q
    st.state.system.q̇ .= 0.0
    lazy_update_bodies!(st)
end

function lazy_update_bodies!(st,q,q̇)
    st.state.system.q .= q
    st.state.system.q̇ .= q̇
    lazy_update_bodies!(st)
end

function lazy_update_bodies!(st)
    (;bodies,state) = st
    foreach(bodies) do body
        bodyid = body.prop.id
        (;q, q̇) = state.members[bodyid]
        lazy_update_state!(body,q,q̇)
        update_transformations!(body,q)
        update_loci_states!(body,q,q̇)
    end
end

function update_bodies!(st,q)
    st.state.system.q .= q
    st.state.system.q̇ .= 0.0
    update_bodies!(st)
end

function update_bodies!(st,q,q̇)
    st.state.system.q .= q
    st.state.system.q̇ .= q̇
    update_bodies!(st)
end

function update_bodies!(st)
    (;bodies,state) = st
    foreach(bodies) do body
        bodyid = body.prop.id
        (;q, q̇) = state.members[bodyid]
        update_cache!(body,q,q̇)
        update_state!(body,q,q̇)
        update_transformations!(body,q)
        update_loci_states!(body,q,q̇)
    end
end

"""
计算System 力的大小。
$(TYPEDSIGNATURES)
"""
function assemble_force!(st::Structure)
    (;bodies,state) = st
    (;system,members) = state
    system.F .= 0.0
    foreach(bodies) do body
        (;F) = members[body.prop.id]
        generalize_force!(F,body.state,body.cache)
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
施加重力。
$(TYPEDSIGNATURES)
"""
function apply_gravity!(st;factor=1)
    (;bodies) = st
    gravity_acceleration = factor*get_gravity(st)
    foreach(bodies) do body
        body.state.mass_locus_state.force .+= gravity_acceleration*body.prop.mass
    end
end