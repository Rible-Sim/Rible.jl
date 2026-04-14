
"""
Robot Type.
$(TYPEDEF)
"""
struct Robot{stType,hubType,trajType,contacts_trajType,contact_caches_trajType,control_trajType}
    structure::stType
    hub::hubType
    traj::trajType
    contacts_traj::contacts_trajType
    contact_caches_traj::contact_caches_trajType
    control_traj::control_trajType
end

"""
Robot Type Constructor.
$(TYPEDSIGNATURES)
"""
function Robot(structure,hub=ControlHub(
            structure,
            Int[], #captum gauges
            Int[], #error gauges
            Int[], #actuators
            Coalition(structure,Int[],Int[],Int[])
        )
    )
    (;bodyid2sys_locus_id) = structure.connectivity
    
    update!(structure)
    
    traj = StructArray([deepcopy(structure.state.system)])
    contacts_traj = [
        reduce(
            vcat,
            map(get_bodies(structure)) do body
                (;prop,) = body
                (;loci) = prop
                bodyid = prop.id
                [
                    Contact(id,lo.friction_coefficient,lo.restitution_coefficient)
                    for (id,lo) in  zip(bodyid2sys_locus_id[bodyid],loci)
                ]
            end
        )
    ]
    contact_caches_traj = [
        activate_frictional_contacts!(
            structure,
            EmptyEnv(),
            traj.q[begin] #set q̇ to zero
        )
    ]
    #note restore q̇ 
    structure.state.system.q̇ .= traj.q̇[begin]

    control_traj = StructArray([deepcopy(hub.state)])
    
    bot = Robot(structure,hub,traj,contacts_traj,contact_caches_traj,control_traj)

    bot 
end


make_cstr_jacobian(bot::Robot) = (q) -> cstr_jacobian(bot.structure,q)
make_cstr_function(bot::Robot) = (q) -> cstr_function(bot.structure,q)

build_mass_matrices(bot::Robot) = build_mass_matrices(bot.structure)



"""
Reset System State.
$(TYPEDSIGNATURES)
"""
function reset!(bot::Robot)
    (;structure, hub, traj) = bot
    (;q, q̇, p, s, c) = traj
    inst_state = structure.state.system
    inst_state.q .= q[begin]
    inst_state.q̇ .= q̇[begin]
    inst_state.p .= p[begin]
    inst_state.s .= s[begin]
    inst_state.c .= c[begin]
    clear_forces!(structure)
    stretch!(structure,c[begin])
    update_bodies!(structure, inst_state)
    update_apparatuses!(structure, inst_state)
    (;M) = build_mass_matrices(structure)
    p[begin] .= M*q̇[begin]
    hub.state.e .= Inf
    hub.state.u .= 0.0
end

"""
Set new nitial conditions.
$(TYPEDSIGNATURES)
"""
function set_initial!(bot::Robot,q,q̇,s=bot.traj.s[begin];)
    (;structure, traj) = bot
    traj.q[begin] .= q
    traj.q̇[begin] .= q̇
    traj.s[begin] .= s
    reset!(bot)
end

function set_initial!(bot::Robot,partial_initial; idx=:)
    (;structure, traj) = bot
    # nq = get_num_of_coords(structure)
    # ns = get_num_of_aux_var(structure)
    initial = ComponentArray(;
        q = traj.q[begin],
        q̇ = traj.q̇[begin],
        s = traj.s[begin]
    )
    initial[idx] .= partial_initial
    traj.q[begin] .= initial.q
    traj.q̇[begin] .= initial.q̇
    traj.s[begin] .= initial.s
    reset!(bot)
end

function set_structure_params!(bot::Robot,c; idx=:)
    (;traj) = bot
    traj.c[begin][idx] .= c
    reset!(bot)
end

"""
Set robot state from a dictionary of CartesianFrames mapped to body IDs.
$(TYPEDSIGNATURES)
"""
function set_robot_state!(bot::Robot, frames::AbstractVector{<:CartesianFrame})
    (;structure, traj) = bot
    bodyid2sys_full_coords = structure.connectivity.bodyid2sys_full_coords
    
    foreach(structure.bodies) do body
        id = body.prop.id
        frame = frames[id]
        q, q̇ = cartesian_frame2coords(body.coords, frame)
        sys_idx = bodyid2sys_full_coords[id]
        traj.q[begin][sys_idx] .= q
        traj.q̇[begin][sys_idx] .= q̇
    end
    structure.state.system.q .= traj.q[begin]
    structure.state.system.q̇ .= traj.q̇[begin]
    update!(structure)
    reset!(bot)
end


"""
Update System State to specific step.
$(TYPEDSIGNATURES)
"""
function goto_step!(bot::Robot,that_step;policy=NoPolicy(),)
    (;structure, traj) = bot
    structure.state.system.c .= traj.c[that_step]
    structure.state.system.q .= traj.q[that_step]
    structure.state.system.q̇ .= traj.q̇[that_step]
    structure.state.system.t = traj.t[that_step]
    actuate!(bot,policy,structure.state.system)
    update!(structure)
    bot
end

get_num_of_dims(bot::Robot) = get_num_of_dims(bot.structure)

get_numbertype(bot::Robot) = get_numbertype(bot.structure)

has_constant_mass_matrix(bot::Robot) = has_constant_mass_matrix(bot.structure)

get_gravity(bot::Robot) = get_gravity(bot.structure)

get_bodies(bot::Robot) = get_bodies(bot.structure)

get_apparatuses(bot::Robot) = get_apparatuses(bot.structure)

get_params(bot::Robot) = get_params(bot.structure)

get_auxiliary(bot::Robot) = get_auxiliary(bot.structure)

get_num_of_capta(bot::Robot) = bot.hub.coalition.num_of_capta