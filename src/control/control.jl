
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
        state = ComponentArray(
            @eponymtuple(
                t = 0.0,
                e,
                u,
            )
        )
        new{typeof(gauges),typeof(actuators),typeof(coalition),typeof(state)}(
            gauges,actuators,coalition,state
        )
    end
end

function get_errors(structure,gauges,coalition)
    (;num_of_errors, gaugeid2error_idx) = coalition.nt
    T = get_numbertype(structure)
    e = zeros(T,num_of_errors)
    foreach(gauges) do gauge
        err_idx = gaugeid2error_idx[gauge.id]
        e[err_idx] .= measure(structure,gauge)
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
end

function cost!(bot::Robot,q::AbstractVector,q̇::AbstractArray,t)
    ## (;q,p,λ) = x
    (;hub) = bot
    (;gauges) = hub
    T = get_numbertype(bot)
    update!(bot.structure,q,q̇)
    ϕ = zero(T)
    foreach(gauges) do gauge
        ϕ += measure(bot,gauge)
    end
    ϕ
end

function cost_jacobian!(∂ϕ∂qᵀ,∂ϕ∂q̇ᵀ,bot::Robot,q::AbstractArray,q̇::AbstractArray,t)
    (;structure,hub) = bot
    (;gauges) = hub
    update!(bot.structure,q)
    T = get_numbertype(bot)
    nq = get_num_of_free_coords(structure)
    nλ = get_num_of_cstr(structure)
    nx = 2nq+nλ
    ## nu = get_num_of_actions(bot)
    ## ∂ϕ∂uᵀ = zeros(T,nu)
    foreach(gauges) do gauge
        measure_jacobian!(∂ϕ∂qᵀ,∂ϕ∂q̇ᵀ,bot,gauge)
    end
    vcat(
        ∂ϕ∂qᵀ, 
        ## ∂ϕ∂uᵀ
    )
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
        ManualActuator(1,collect(1:nμ),zeros(nμ),Uncoupled()),
        itp
    )
end

get_actions(bot::Robot) = get_actions(bot.structure,bot.hub.actuators,bot.hub.coalition)

function get_actions(bot::Robot,q,q̇,t)
    u
end