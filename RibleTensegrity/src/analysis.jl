

function string_potential(tgstruct,sol)
    [
    [begin
        reset_forces!(tgstruct)
        distribute_q_to_rbs!(tgstruct,q)
        update_cables_apply_forces!(tgstruct)
        potential_energy(s)
    end
    for q in sol.qs] for s in tgstruct.cables]
end

function analyse_slackness(tgstruct,q)
    distribute_q_to_rbs!(tgstruct,q)
    for (i,s) in enumerate(tgstruct.cables)
        len = s.state.length
        restlen = s.state.restlen
        Δl = len - restlen
        if Δl < 0
            @info "String $i is slack."
        else
            @info "Slackness check pass. Δl_$i = $Δl."

        end
    end
end

function test_slackness(tr)
    (;structure, traj) = tr
     "Test Slackness" 
     begin
         for (j,q) in enumerate(traj.qs)
            "Solution $j" 
            distribute_q_to_rbs!(structure,q)
             for (i,s) in enumerate(structure.cables)
                "String $i" 
                len = s.state.length
                restlen = s.state.restlen
                @assert  len > restlen
            end
        end
    end
end


function analyse_slack(st::AbstractStructure,verbose=false)
    (;cables) = st.apparatuses
    slackcases = [cable.id for cable in cables if cable.slack && (cable.state.length <= cable.state.restlen)]
    if verbose && !isempty(slackcases)
        @show slackcases
    end
    slackcases
end


get_cables_len(bot::Robot) = get_cables_len(bot.structure)
get_cables_deform(bot::Robot) = get_cables_deform(bot.structure)
get_cables_restlen(bot::Robot) = get_cables_restlen(bot.structure)
get_cables_len_dot(bot::Robot) = get_cables_len_dot(bot.structure)
get_cables_tension(bot::Robot) = get_cables_tension(bot.structure)
get_cables_stiffness(bot::Robot) = get_cables_stiffness(bot.structure)
get_cables_force_density(bot::Robot) = get_cables_force_density(bot.structure)
get_cables(bot::Robot) = get_cables(bot.structure)
get_cables(st::AbstractStructure) = get_cables(st.apparatuses)

function get_cables_len!(st::TensegrityStructure,q)
    update_bodies!(st::AbstractStructure,q,zero(q))
    update_cables_apply_forces!(st::TensegrityStructure)
    get_cables_len(st::TensegrityStructure)
end

"""
Return System DistanceSpringDamper stiffness coefficients
$(TYPEDSIGNATURES)
"""
function get_cables_stiffness(st::TensegrityStructure)
    cables = get_cables(st)
    [s.force.k for s in cables]
end

"""
Return System DistanceSpringDamper current Length.
$(TYPEDSIGNATURES)
"""
function get_cables_len(st::TensegrityStructure)
    cables = get_cables(st)
    [s.force.state.length for s in cables]
end

function get_cables_len_dot(st::TensegrityStructure)
    cables = get_cables(st)
    [s.force.state.lengthdot for s in cables]
end

"""
Return System DistanceSpringDamper deformation.
$(TYPEDSIGNATURES)
"""
function get_cables_deform(st::TensegrityStructure)
    cables = get_cables(st)
    [s.force.state.length - s.force.state.restlen for s in cables]
end

"""
Return System DistanceSpringDamper Restlength.
$(TYPEDSIGNATURES)
"""
function get_cables_restlen(st::TensegrityStructure)
    cables = get_cables(st)
    [s.force.state.restlen for s in cables]
end

"""
Return System DistanceSpringDamper Tension.
$(TYPEDSIGNATURES)
"""
function get_cables_tension(st::TensegrityStructure)
    cables = get_cables(st)
    [s.force.state.tension for s in cables]
end

"""
Set cables' tension
$(TYPEDSIGNATURES)
"""
function set_cables_tension!(st::TensegrityStructure,fs)
    cables = get_cables(st)
    for (s,f) in zip(cables,fs)
        s.force.state.tension = f
    end
end


"""
Return System DistanceSpringDamper force density.
$(TYPEDSIGNATURES)
"""
function get_cables_force_density(st::TensegrityStructure)
    cables = get_cables(st)
    [s.force.state.tension/s.force.state.length for s in cables]
end


"""
Return System DistanceSpringDamper Initial Length.
$(TYPEDSIGNATURES)
"""
function get_original_restlen(botinput::Robot)
    bot = deepcopy(botinput)
    T = get_numbertype(bot)
    actuate!(bot,zeros(T,length(bot.hub.actuators)))
    u0 = get_cables_restlen(bot.structure)
end

function force_densities_to_restlen(st::TensegrityStructure,γs)
    cables = get_cables(st)
    [
    begin
        l = s.force.state.length
        l̇ = s.force.state.lengthdot
        k = s.k
        c = s.c
        u = l-(γ*l-c*l̇)/k
    end
        for (γ,s) in zip(γs,cables)]
end


function get_tension!(bot::Robot,cid::Int,step_range=:)
    (; structure, traj)= bot
    (; cables) = structure.apparatuses
    T = get_numbertype(structure)
    f = Vector{T}()
    h = traj.t[2] - traj.t[1]
    q_mids = [(traj.q[k] .+ traj.q[k+1])./2 for k = 1:length(traj)-1]
    q̇_mids = [(traj.q[k] .- traj.q[k+1])./h for k = 1:length(traj)-1]
    for (q,q̇) in zip(q_mids, q̇_mids)
        update!(structure,q)
        push!(f,cables[cid].state.tension)
    end
    f
end

function test_fvector(st,q0)
    function L(q)
        reset_forces!(st)
        distribute_q_to_rbs!(st,q,zero(q))
        update_cables_apply_forces!(st)
        fvector(st)
        [st.cables[i].state.length for i = 1:2]
    end
    FiniteDiff.finite_difference_jacobian(L,q0)
end

get_cables(apparatuses::AbstractVector) = filter!((x)->x.joint isa CableJoint, sort(apparatuses))

get_cables(apparatuses::TypeSortedCollection) = filter!((x)->x.joint isa CableJoint, sort(apparatuses))
