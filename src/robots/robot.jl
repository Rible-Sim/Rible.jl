
"""
Robot Type.
$(TYPEDEF)
"""
struct Robot{tgT,hubT,trajT,contacts_trajT}
    st::tgT
    hub::hubT
    traj::trajT
    contacts_traj::contacts_trajT
end


"""
Robot Type Constructor.
$(TYPEDSIGNATURES)
"""
function Robot(st,hub=nothing)
    (;numbered) = st.connectivity
    (;mem2num) = numbered
    update!(st)
    traj = StructArray([deepcopy(st.state.system)])
    
    contacts_traj = [
        reduce(
            vcat,
            map(get_bodies(st)) do body
                (;prop,) = body
                (;loci) = prop
                bodyid = prop.id
                [
                    Contact(id,lo.friction_coefficient,lo.restitution_coefficient)
                    for (id,lo) in  zip(mem2num[bodyid],loci)
                ]
            end
        )
    ]
    Robot(st,hub,traj,contacts_traj)
end


make_constraints_jacobian(bot::Robot) = make_constraints_jacobian(bot.st)
make_constraints_function(bot::Robot) = make_constraints_function(bot.st)

"""
Return System 质量矩阵。
$(TYPEDSIGNATURES)
"""
function build_MassMatrices(bot::Robot)
    (;st) = bot
    (;nfree,npres,sysfree,syspres) = st.connectivity.indexed
    M = build_M(st)
    Ḿ = M[sysfree,:]
    M̌ = Symmetric(M[sysfree,sysfree])
    M̄ =           M[sysfree,syspres]
    invM̌_raw = inv(Matrix(M̌))
    invM̌ = Symmetric(sparse(invM̌_raw))
    @eponymtuple(Ḿ,M̌,M̄,invM̌)
end


function build_Y(bot)
    (;st, hub) = bot
    (;actuators) = hub
    (;cables) = st.tensiles
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
    (;st, traj) = bot
    (;q, q̇) = traj
    clear_forces!(st)
    update_rigids!(st,q[begin],q̇[begin])
    update_tensiles!(st)
    resize!(traj,1)
end

"""
更改Initial 条件。
$(TYPEDSIGNATURES)
"""
function set_new_initial!(bot::Robot,q,q̇=zero(q))
    (;st, traj) = bot
    traj.q[begin] .= q
    traj.q̇[begin] .= q̇
    reset!(bot)
end

"""
Update System 到指定时间步State.
$(TYPEDSIGNATURES)
"""
function goto_step!(bot::Robot,that_step;actuate=false)
    (;st, traj) = bot
    st.state.system.c .= traj.c[that_step]
    st.state.system.q .= traj.q[that_step]
    st.state.system.q̇ .= traj.q̇[that_step]
    if actuate
        actuate!(bot,[traj.t[that_step]])
    end
    update!(st)
    bot
end


function mechanical_energy!(bot::Robot;actuate=false,gravity=true)
    (;st,traj) = bot
    StructArray([
        begin
            st.state.system.q .= trajstate.q
            st.state.system.q̇ .= trajstate.q̇
            if actuate
                actuate!(bot,[trajstate.t])
            end
            update!(st)
            mechanical_energy(st;gravity)
        end
        for trajstate in traj
    ])
end
