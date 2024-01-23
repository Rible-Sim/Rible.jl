function actuate!(bot::Robot,u::Nothing,t::Number) end

function actuate!(bot::Robot,u::AbstractVector,t::Number)
    (;actuators) = bot.hub
    foreach(actuators) do actuator
        execute!(bot,actuator)
    end
end

function action_jacobian(bot::Robot,q,q̇,λ,u)
    (;num_of_actions) = bot.hub.scheme
    ∂f∂u = zeros(T,nu)
    foreach(actuators) do actuator
        ∂f∂u .+= action_jacobian(bot,actuator)
    end
end

abstract type AbstractOperator end
struct ManualOperator <: AbstractOperator end
struct TimeOperator <: AbstractOperator end

abstract type AbstractActuator end

get_id(actuator::AbstractActuator) = actuator.id
get_numbertype(actuator::AbstractActuator) = get_numbertype(actuator.signifier)

struct ExternalForceActuator{sigType,operType,forceType} <: AbstractActuator
    id::Int
    signifier::sigType
    operator::operType
    force::forceType
end

get_num_of_actions(::ExternalForceActuator) = 1

function generalized_force(structure::Structure,actuator::ExternalForceActuator{sigType,<:TimeOperator}) where {sigType}
    (;state) = structure
    (;t) = state.system
    (;signifier,operator,force) = actuator
    (;body,pid,) = signifier
    (;prop,coords) = body
    (;nmcs) = coords
    bid = body.prop.id
    (;q) = state.members[bid]
    c = to_local_coords(nmcs,prop.loci[pid].position)
    Tbody = build_T(structure,bid)
    C = to_position_jacobian(nmcs,q,c)*Tbody
    transpose(C)*force(t)
end

function generalized_force_jacobian(structure::Structure,actuator::ExternalForceActuator{sigType,<:TimeOperator}) where {sigType}
    (;state) = structure
    (;t) = state.system
    (;signifier,operator) = actuator
    (;body,pid,) = signifier
    (;prop,coords) = body
    (;nmcs) = coords
    bid = body.prop.id
    (;q) = state.members[bid]
    c = to_local_coords(nmcs,prop.loci[pid].position)
    Tbody = build_T(structure,bid)
    C = to_position_jacobian(nmcs,q,c)*Tbody
    transpose(C)*force(t)
end

function execute!(structure::Structure,actuator::ExternalForceActuator{sigType,<:TimeOperator}) where {sigType}
    (;id,signifier,operator,force) = actuator
    structure.state.system.F .+= generalized_force(structure,actuator)
end

function get_actions(structure::Structure,actuator::ExternalForceActuator{sigType,<:TimeOperator}) where {sigType}
    (;signifier,operator,force) = actuator
    (;state) = structure
    (;t) = state.system
    t
end

function generalized_force(actuator::ExternalForceActuator{sigType,<:ManualOperator}) where {sigType}
    (;state) = structure
    (;t) = state.system
    (;signifier,operator,force) = actuator
    force(signifier,t)
end

function generalized_force_jacobian(structure::Structure,actuator::ExternalForceActuator{sigType,<:ManualOperator}) where {sigType}
    (;state) = structure
    (;t) = state.system
    (;signifier,operator) = actuator
    (;body,pid,) = signifier
    (;prop,coords) = body
    (;nmcs) = coords
    bid = body.prop.id
    (;q) = state.members[bid]
    c = to_local_coords(nmcs,prop.loci[pid].position)
    Tbody = build_T(structure,bid)
    C = to_position_jacobian(nmcs,q,c)*Tbody
    transpose(C)*operator(t)
    force(signifier,t)
end

function execute!(structure::Structure,actuator::ExternalForceActuator{sigType,<:ManualOperator}) where {sigType}
    (;id,signifier,operator,force) = actuator
    structure.state.system.F .+= generalized_force(actuator)
end

function get_actions(structure::Structure,actuator::ExternalForceActuator{sigType,<:ManualOperator}) where {sigType}
    (;signifier,operator,force) = actuator
    (;state) = structure
    (;t) = state.system
    t
end

struct RegisterActuator{RT,CT,registerType} <: AbstractActuator
    id::Int
    signifier::RT
    operator::CT
    register::registerType
end

get_num_of_actions(actuator::RegisterActuator) = size(actuator.register.matrix,2)

function get_actions(structure::Structure,actuator::RegisterActuator)
    (;signifier,operator,register) = actuator
    zeros(eltype(register.values),get_num_of_actions(actuator))
end

function execute!(structure::Structure,actuator::RegisterActuator)
    (;id,signifier,register) = actuator
    (;values) = register
    foreach(signifier) do apparatus
        apparatus.force.state.restlen = original_restlen + u[id]
    end
end

# heater
struct SMAHeater{CT,operType,HT,TT} <:AbstractActuator
    id::Int
    signifier::CT
    operator::operType
    heating_law::HT
end

function SMAHeater(signifier::Signifier,heating_law)
    SMAHeater(signifier,heating_law,traj)
end

function execute!(actuator::SMAHeater,structure,u;inc=false,abs=true)
    (;SMA_cables) = structure.apparatuses
    (;id_string, original_value, heating_law) = actuator
    s = (SMA_cables,id_string)
    if abs
        s.state.temp = u
    else
        if inc
            s.state.temp += u
        else
            s.state.temp = original_value + u
        end
    end
    execute!(s,heating_law)
end

function execute!(actuator::SMAHeater,heating_law)
    s.law.F0, s.law.k = heating_law(s.state.temp)
end
