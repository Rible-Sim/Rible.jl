
"""
Robot Type.
$(TYPEDEF)
"""
struct Robot{stType,cacheType,hubType,trajType,contacts_trajType}
    structure::stType
    structure_cache::cacheType
    control_hub::hubType
    traj::trajType
    contacts_traj::contacts_trajType
end

"""
Robot Type Constructor.
$(TYPEDSIGNATURES)
"""
function Robot(structure,hub=nothing)
    (;numbered) = structure.connectivity
    (;bodyid2sys_loci_idx) = numbered
    structure_cache = StructureCache(structure)
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
                    for (id,lo) in  zip(bodyid2sys_loci_idx[bodyid],loci)
                ]
            end
        )
    ]
    Robot(structure,structure_cache,hub,traj,contacts_traj)
end


make_cstr_jacobian(bot::Robot) = make_cstr_jacobian(bot.structure)
make_cstr_function(bot::Robot) = make_cstr_function(bot.structure)

"""
Return System 质量矩阵。
$(TYPEDSIGNATURES)
"""
function build_mass_matrices(bot::Robot)
    (;structure) = bot
    (;num_of_free_coords,num_of_pres_coords,sys_free_idx,sys_pres_idx) = structure.connectivity.indexed
    M = assemble_M(structure)
    Ḿ = M[sys_free_idx,:]
    M̌ = Symmetric(M[sys_free_idx,sys_free_idx])
    M̄ =           M[sys_free_idx,sys_pres_idx]
    invM̌_raw = inv(Matrix(M̌))
    invM̌ = Symmetric(sparse(invM̌_raw))
    @eponymtuple(Ḿ,M̌,M̄,invM̌)
end


function build_Y(bot)
    (;structure, hub) = bot
    (;actuators) = hub
    (;cables) = structure.tensiles
    ncables = length(cables)
    nact = length(actuators)
    ret = spzeros(Int,ncables,nact)
    foreach(actuators) do act
        (;id,coupler,reg) = act
        if coupler isa Serial
            ret[act.reg.ids,id] .= 1
        elseif coupler isa Ganged
            is1,is2 = act.reg.ids
            ret[is1,id] =  1
            ret[is2,id] = -1
        else
            error("Unknown actuator type")
        end
    end
    ret
end

"""
重置System State.
$(TYPEDSIGNATURES)
"""
function reset!(bot::Robot)
    (;structure, traj) = bot
    (;q, q̇) = traj
    clear_forces!(structure)
    update_bodies!(structure,q[begin],q̇[begin])
    update_tensiles!(structure)
    resize!(traj,1)
end

"""
更改Initial 条件。
$(TYPEDSIGNATURES)
"""
function set_new_initial!(bot::Robot,q,q̇=zero(q))
    (;structure, traj) = bot
    traj.q[begin] .= q
    traj.q̇[begin] .= q̇
    reset!(bot)
end

"""
Update System 到指定时间步State.
$(TYPEDSIGNATURES)
"""
function goto_step!(bot::Robot,that_step;actuate=false)
    (;structure, traj) = bot
    structure.state.system.c .= traj.c[that_step]
    structure.state.system.q .= traj.q[that_step]
    structure.state.system.q̇ .= traj.q̇[that_step]
    if actuate
        actuate!(bot,[traj.t[that_step]])
    end
    update!(structure)
    bot
end


function mechanical_energy!(bot::Robot;actuate=false,gravity=true)
    (;structure,traj) = bot
    StructArray([
        begin
            structure.state.system.q .= trajstate.q
            structure.state.system.q̇ .= trajstate.q̇
            if actuate
                actuate!(bot,[trajstate.t])
            end
            update!(structure)
            mechanical_energy(structure;gravity)
        end
        for trajstate in traj
    ])
end
