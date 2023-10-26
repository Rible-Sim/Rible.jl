
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
        state1 = scnt.hen.rbsig.state
        state2 = scnt.egg.rbsig.state
        pid1 = scnt.hen.pid
        pid2 = scnt.egg.pid
        p1 = state1.rps[pid1]
        ṗ1 = state1.ṙps[pid1]
        f1 = state1.fps[pid1]
        p2 = state2.rps[pid2]
        ṗ2 = state2.ṙps[pid2]
        f2 = state2.fps[pid2]
        update!(scable,p1,p2,ṗ1,ṗ2)
        f1 .+= scable.state.force
        f2 .-= scable.state.force
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
    foreach(bodies) do rb
        rbid = rb.prop.id
        stretch_rigid!(rb,c[mem2sys[rbid]])
    end
end

function move_rigids!(st,q,q̇=zero(q))
    st.state.system.q .= q
    st.state.system.q̇ .= q̇
    move_rigids!(st)
end

function move_rigids!(st)
    (;bodies,state) = st
    foreach(bodies) do rb
        rbid = rb.prop.id
        (;q,q̇) = state.parts[rbid]
        move_rigid!(rb,q,q̇)
    end
end

"""
Update Rigid Body State 。
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
    foreach(bodies) do rb
        rbid = rb.prop.id
        (;q, q̇) = state.parts[rbid]
        update_rigid!(rb,q,q̇)
        update_transformations!(rb,q)
        move_rigid!(rb,q,q̇)
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
    foreach(bodies) do rb
        (;F) = parts[rb.prop.id]
        generalize_force!(F,rb.state)
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
    foreach(bodies) do rb
        rb.state.f .+= gravity_acceleration*rb.prop.mass
    end
end