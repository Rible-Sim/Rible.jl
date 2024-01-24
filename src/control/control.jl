
abstract type Controller end

include("gauges.jl")
include("actuators.jl")
include("PIDController.jl")

struct Coalition{ntType}
    nt::ntType
end

function Coalition(structure::AbstractStructure,gauges,actuators)
    gauges_ids,num_of_gauges = check_id_sanity(gauges) 
    actuators_ids,num_of_actuators = check_id_sanity(actuators)
    num_of_error_gauges = 0
    num_of_capta_gauges = 0
    sys_errors_idx = Int[]
    gaugeid2error_idx = Vector{Vector{Int}}(undef,num_of_gauges)
    if num_of_gauges > 0
        foreach(gauges) do gauge
            gid = gauge.id
            nerr = get_num_of_errors(gauge)
            gaugeid2error_idx[gid] = collect(length(sys_errors_idx)+1:length(sys_errors_idx)+nerr)
            append!(sys_errors_idx,gaugeid2error_idx[gid])
            if !(gauge isa ErrorGauge)
                num_of_error_gauges += 1
            end
            if !(gauge isa CaptaGauge)
                num_of_capta_gauges += 1
            end
        end
    end

    num_of_errors = length(sys_errors_idx)

    sys_actions_idx = Int[]
    actid2sys_actions = Vector{Int}[]
    ntotal_by_act = zeros(Int,num_of_actuators)
    pres_idx_by_act = Vector{Vector{Int}}(undef,num_of_actuators)
    free_idx_by_act = Vector{Vector{Int}}(undef,num_of_actuators)
    foreach(actuators) do actuator
        actid = actuator.id
        ntotal_by_act[actid] = get_num_of_actions(actuator)
        pres_idx_by_act[actid] = Int[]
        free_idx_by_act[actid] = collect(1:ntotal_by_act[actid])
    end
    for actid = 1:num_of_actuators
        ntotal = ntotal_by_act[actid]
        push!(actid2sys_actions,fill(-1,ntotal))
        unshareds = collect(1:ntotal)
        nusi = length(unshareds)
        actid2sys_actions[actid][unshareds] = collect(length(sys_actions_idx)+1:length(sys_actions_idx)+nusi)
        append!(sys_actions_idx,actid2sys_actions[actid][unshareds])
    end

    num_of_actions = length(sys_actions_idx)

    Coalition(
        @eponymtuple(
            num_of_actuators, num_of_gauges,
            num_of_error_gauges, num_of_capta_gauges,
            num_of_errors,gaugeid2error_idx,
            num_of_actions, actid2sys_actions, 
        )
    )
end

struct ControlHub{gaugesType,actuatorsType,coalitionType,stateType}
    gauges::gaugesType
    actuators::actuatorsType
    coalition::coalitionType
    state::stateType
    function ControlHub(structure::AbstractStructure,gauges,actuators,coalition::Coalition)
        e = get_errors(structure,gauges,coalition)
        u = get_actions(structure,actuators,coalition)
        state = @eponymtuple(
                e,
                u,
            )
        new{typeof(gauges),typeof(actuators),typeof(coalition),typeof(state)}(
            gauges,actuators,coalition,state
        )
    end
end

function get_actions(structure,actuators,coalition)
    (;num_of_actions,actid2sys_actions) = coalition.nt
    T = get_numbertype(structure)
    u = zeros(T,num_of_actions)
    foreach(actuators) do actuator
        act_idx = actid2sys_actions[actuator.id]
        u[act_idx] .= get_actions(structure,actuator)
    end
    u
end

function get_actions!(bot::Robot,t::Number)
    (;structure,hub) = bot
    (;actuators,coalition) = hub
    (;num_of_actions,actid2sys_actions) = coalition.nt
    structure.state.system.t = t
    T = get_numbertype(structure)
    u = zeros(T,num_of_actions)
    foreach(actuators) do actuator
        act_idx = actid2sys_actions[actuator.id]
        u[act_idx] .= get_actions(structure,actuator)
    end
    u
end

function actuate!(bot::Robot,q::AbstractVector,q̇::AbstractVector,u::AbstractVector,t)
    (;structure,hub) = bot
    (;actuators) = hub
    foreach(actuators) do actuator
        execute!(structure,actuator)
    end
    update!(structure,q,q̇,t)
end

function actuate!(bot::Robot,t::Number=bot.structure.state.system.t)
    (;structure,hub) = bot
    bot.structure.state.system.t = t
    (;actuators) = hub
    foreach(actuators) do actuator
        execute!(structure,actuator)
    end
end

# For now, the cost is simply the sum of errors
#todo weighted sum
#todo user-defined cost
function cost!(bot::Robot,q::AbstractVector,q̇::AbstractVector,u::AbstractVector,t)
    ## (;q,p,λ) = x
    (;structure,hub) = bot
    (;gauges) = hub
    T = get_numbertype(bot)
    actuate!(bot,q,q̇,u,t)
    ϕ = zero(T)
    foreach(gauges) do gauge
        ϕ += measure(structure,gauge)
    end
    ϕ
end

function cost!(bot::Robot,)
    ## (;q,p,λ) = x
    (;q,q̇,t) = bot.structure.state.system
    (;u) = bot.hub.state
    cost!(bot,q,q̇,u,t)
end

function cost_jacobian!(∂ϕ∂qᵀ,∂ϕ∂q̇ᵀ,bot::Robot,q::AbstractVector,q̇::AbstractVector,t)
    (;structure,hub) = bot
    (;coalition) = hub
    update!(structure,q,q̇,t)
    ∂e∂q, ∂e∂q̇ = errors_jacobian(bot,)
    ∂ϕ∂qᵀ .= transpose(sum(∂e∂q,dims=1))
    ∂ϕ∂q̇ᵀ .= transpose(sum(∂e∂q̇, dims=1))
end

function cost_jacobian!(∂ϕ∂uᵀ,bot::Robot,q::AbstractVector,q̇::AbstractVector,u::AbstractVector,t)
    (;structure,hub) = bot
    (;coalition) = hub
    actuate!(bot,q,q̇,u,t)
    (;num_of_errors, gaugeid2error_idx, num_of_actions, actid2sys_actions) = coalition.nt
    ∂e∂q, ∂e∂q̇ = errors_jacobian(bot,)
    ∂ϕ∂uᵀ .= transpose(sum(∂e∂q,dims=1))
end

function cost_jacobian!(bot::Robot)
    (;q,q̇,t) = bot.structure.state.system
    (;u) = bot.hub.state
    T = get_numbertype(bot)
    ∂ϕ∂qᵀ = zeros(T,length(q))
    ∂ϕ∂q̇ᵀ = zeros(T,length(q̇))
    ∂ϕ∂uᵀ = zeros(T,length(u))
    cost_jacobian!(∂ϕ∂qᵀ,∂ϕ∂q̇ᵀ,∂ϕ∂uᵀ,bot,q,q̇,u,t)
end

function set_restlen!(structure,u)
    for (i,s) in enumerate(structure.apparatuses.apparatuses)
        s.state.restlen = u[i]
    end
end

function make_pres_actor(μ0,μ1,start,stop)
    nμ = length(μ0)

    function itp(t)
        scaled_itps = extrapolate(
            Interpolations.scale(
                interpolate(
                    hcat(μ0,μ1),
                    (NoInterp(),BSpline(Linear()))
                    # (NoInterp(),BSpline(Quadratic(Flat(OnGrid()))))
                ),
                1:nμ, start:stop-start:stop
            ),
            (Throw(),Flat())
        )
        [scaled_itps(j,t) for j in 1:nμ]
    end

    PrescribedActuator(
        1,
        RegisterActuator(1,collect(1:nμ),zeros(nμ),Uncoupled()),
        itp
    )
end
