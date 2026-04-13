
abstract type Controller end

include("actuators.jl")
include("reg_actuator.jl")

include("capta.jl")
include("gauges.jl")

struct Coalition
    num_of_actuators::Int
    num_of_capta_gauges::Int
    num_of_error_gauges::Int
    num_of_capta::Int
    num_of_errors::Int
    captum_gauge_id2sys_capta_idx::Vector{Vector{Int}}
    error_gauge_id2sys_errors_idx::Vector{Vector{Int}}
    num_of_actions::Int 
    actid2sys_actions_idx::Vector{Vector{Int}} 
end

function Coalition(structure::AbstractStructure,capta_gauges,error_gauges,actuators)
    capta_gauges_ids,num_of_capta_gauges = check_id_sanity(capta_gauges) 
    error_gauges_ids,num_of_error_gauges = check_id_sanity(error_gauges) 
    actuators_ids,num_of_actuators = check_id_sanity(actuators)
    
    # Process captum gauges - for feedback policy (captum -> action)
    num_of_capta = 0
    sys_capta_idx = Int[]
    num_of_capta_per_gauge = Vector{Int}(undef,num_of_capta_gauges)
    captum_gauge_id2sys_capta_idx = Vector{Vector{Int}}(undef,num_of_capta_gauges)
    if num_of_capta_gauges > 0
        foreach(capta_gauges) do gauge
            gid = gauge.id
            # For captum gauges, we track each measurement component (not errors)
            num_of_capta_per_gauge[gid] = get_num_of_capta(gauge)
            captum_gauge_id2sys_capta_idx[gid] = collect(
                length(sys_capta_idx)+1:length(sys_capta_idx)+num_of_capta_per_gauge[gid]
            )
            append!(sys_capta_idx,captum_gauge_id2sys_capta_idx[gid])
        end
    end
    num_of_capta = length(sys_capta_idx)
    
    # Process error gauges - for cost function (error -> cost)
    sys_errors_idx = Int[]
    num_of_errors_ref = Ref(0)
    error_gauge_id2sys_errors_idx = Vector{Vector{Int}}(undef,num_of_error_gauges)
    if num_of_error_gauges > 0
        foreach(error_gauges) do gauge
            gid = gauge.id
            # ErrorGauge produces 1 scalar error
            error_gauge_id2sys_errors_idx[gid] = collect(gid:gid)
            num_of_errors_ref[] += 1
        end
    end
    num_of_errors = num_of_errors_ref[]

    sys_actions_idx = Int[]
    actid2sys_actions_idx = Vector{Int}[]
    ntotal_by_act = zeros(Int,num_of_actuators)
    foreach(actuators) do actuator
        actid = actuator.id
        ntotal_by_act[actid] = get_num_of_actions(actuator)
    end
    for actid = 1:num_of_actuators
        ntotal = ntotal_by_act[actid]
        push!(actid2sys_actions_idx,fill(-1,ntotal))
        unshareds = collect(1:ntotal)
        nusi = length(unshareds)
        actid2sys_actions_idx[actid][unshareds] = collect(length(sys_actions_idx)+1:length(sys_actions_idx)+nusi)
        append!(sys_actions_idx,actid2sys_actions_idx[actid][unshareds])
    end

    num_of_actions = length(sys_actions_idx)

    Coalition(
        num_of_actuators, 
        num_of_capta_gauges,
        num_of_error_gauges,
        num_of_capta,
        num_of_errors,
        captum_gauge_id2sys_capta_idx,
        error_gauge_id2sys_errors_idx,
        num_of_actions, 
        actid2sys_actions_idx, 
    )
end


struct ControlHub{captumGaugesType,errorGaugesType,actuatorsType,coalitionType,stateType}
    capta_gauges::captumGaugesType
    error_gauges::errorGaugesType
    actuators::actuatorsType
    coalition::coalitionType
    state::stateType
    function ControlHub(structure::AbstractStructure,capta_gauges,error_gauges,actuators,coalition::Coalition,)
        (;num_of_capta, num_of_errors, captum_gauge_id2sys_capta_idx, error_gauge_id2sys_errors_idx) = coalition
        T = get_numbertype(structure)
        
        # Initialize captum measurements
        c = zeros(T,num_of_capta)
        foreach(capta_gauges) do gauge
            capta_idx = captum_gauge_id2sys_capta_idx[gauge.id]
            c[capta_idx] .= measure(structure,gauge)
        end
        
        # Initialize error measurements  
        e = zeros(T,num_of_errors)
        foreach(error_gauges) do gauge
            err_idx = error_gauge_id2sys_errors_idx[gauge.id]
            e[err_idx] .= measure(structure,gauge)
        end
        
        u = get_initial_actions(structure,actuators,coalition)

        state = @eponymtuple(
                c,
                e,
                u,
            )
        new{typeof(capta_gauges),typeof(error_gauges),typeof(actuators),typeof(coalition),typeof(state)}(
            capta_gauges,error_gauges,actuators,coalition,state,
        )
    end
end


get_num_of_actions(bot) = bot.hub.coalition.num_of_actions

function get_initial_actions(structure::AbstractStructure,actuators,coalition::Coalition)
    (;num_of_actions,actid2sys_actions_idx) = coalition
    T = get_numbertype(structure)
    u = zeros(T,num_of_actions)
    foreach(actuators) do actuator
        act_idx = actid2sys_actions_idx[actuator.id]
        u[act_idx] .= get_initial_actions(structure,actuator)
    end
    u
end

function get_initial_actions!(bot::Robot,t::Real=bot.structure.state.system.t)
    (;structure,hub) = bot
    (;actuators,coalition) = hub
    (;num_of_actions,actid2sys_actions_idx) = coalition
    structure.state.system.t = t
    T = get_numbertype(structure)
    u = zeros(T,num_of_actions)
    foreach(actuators) do actuator
        act_idx = actid2sys_actions_idx[actuator.id]
        u[act_idx] .= get_initial_actions(structure,actuator)
    end
    u
end

function set_control_params!(bot::Robot,policy::AbstractPolicy,partial_params; idx=:)
    (;structure, hub, traj) = bot
    params = get_params(policy)
    params[idx] .= partial_params
    set_params!(policy,params)
    update!(policy)
end


function actuate!(bot::Robot,policy::NoPolicy,inst_state::AbstractCoordinatesState)
end


function actuate!(bot::Robot,policy::TimeFunctionPolicy,inst_state::AbstractCoordinatesState)
    (;structure,hub) = bot
    (;t) = inst_state
    structure.state.system.t = t
    (;actuators,state,coalition,) = hub
    (;actid2sys_actions_idx) = coalition
    control = policy.f
    @debug "Actuating TimeFunctionPolicy" length(control(t)) length(state.u)
    state.u .= control(t)
    foreach(actuators) do actuator
        idx = actid2sys_actions_idx[actuator.id]
        execute!(
            structure,
            actuator,
            (@view state.u[idx])
        )
    end
end




"""
    control!(bot::Robot,policy,state::ComponentArray)

For policy that acts directly on discrete solver's discretized states (e.g. iLQR)
"""
function control!(bot::Robot,policy,solver_cache,solver_state)
end

"""
    control_jacobian!(workspace,bot::Robot,policy,solver_cache,solver_state)
Jacobian of the control function w.r.t. the discrete solver's discretized states.
For policy that acts directly on discrete solver's discretized states (e.g. iLQR)
"""
function control_jacobian!(workspace,bot::Robot,policy,solver_cache,solver_state)
end

function vjp_wrt_state(v,policy::NoPolicy,bot::Robot,num_of_actions,solver,solver_state)
    (;structure) = bot
    (;tₖ, tₖ₊₁) = solver_state
    nq = get_num_of_full_coords(structure)
    nλ = get_num_of_cstr(structure)
    ns = get_num_of_aux_var(structure)
    ∂ϕ∂qₖᵀ = spzeros(typeof(tₖ),nq)
    ∂ϕ∂qₖ₊₁ᵀ = spzeros(typeof(tₖ),nq)
    ∂ϕ∂pₖᵀ = spzeros(typeof(tₖ),nq)
    ∂ϕ∂pₖ₊₁ᵀ = spzeros(typeof(tₖ),nq)
    ∂ϕ∂λᵀ = spzeros(typeof(tₖ),nλ)
    ∂ϕ∂sₖᵀ = spzeros(typeof(tₖ),ns)
    ∂ϕ∂sₖ₊₁ᵀ = spzeros(typeof(tₖ),ns)
    ∂ϕ∂qₖᵀ,  ∂ϕ∂qₖ₊₁ᵀ, ∂ϕ∂pₖᵀ, ∂ϕ∂pₖ₊₁ᵀ, ∂ϕ∂λᵀ, ∂ϕ∂sₖᵀ, ∂ϕ∂sₖ₊₁ᵀ
end

function accumulate_param_grad!(grad_storage, policy::NoPolicy, v_total, solver_state, bot) end

function get_params(policy::NoPolicy)
    nothing
end


function get_num_of_params(policy::NoPolicy)
    0
end
function get_num_of_params(policy::TimeFunctionPolicy)
    0
end
