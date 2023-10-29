
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

function update_tensiles!(st, @eponymargs(connected,))
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
    (;mem2sys) = st.connectivity.numbered
    (;c) = state.system
    foreach(bodies) do body
        bodyid = body.prop.id
        stretch_body!(body,c[mem2sys[bodyid]])
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
        move_body!(body,q,q̇)
    end
end

"""
Update Rigid Body State.
$(TYPEDSIGNATURES)
"""
function update_rigids!(st,q)
    st.state.system.q .= q
    st.state.system.q̇ .= 0.0
    update_rigids!(st)
end

function update_rigids!(st,q,q̇)
    st.state.system.q .= q
    st.state.system.q̇ .= q̇
    update_rigids!(st)
end

function update_rigids!(st)
    (;bodies,state) = st
    foreach(bodies) do body
        bodyid = body.prop.id
        (;q, q̇) = state.parts[bodyid]
        update_body!(body,q,q̇)
        update_transformations!(body,q)
        move_body!(body,q,q̇)
    end
end

"""
计算System 力的大小。
$(TYPEDSIGNATURES)
"""
function generate_forces!(st::Structure)
    (;bodies,state) = st
    (;system,parts) = state
    system.F .= 0.0
    foreach(bodies) do body
        (;F) = parts[body.prop.id]
        generalize_force!(F,body.state)
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