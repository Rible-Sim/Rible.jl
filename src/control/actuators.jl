function actuate!(bot::Robot,u::Nothing,t::Number) end

function actuate!(bot::Robot,u::AbstractVector,t::Number)
    (;actuators) = bot.hub
    foreach(actuators) do actuator
        actuate!(bot,actuator)
    end
end

abstract type AbstractOperator end
struct ManualOperator <: AbstractOperator end

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

function get_actions(structure::Structure,actuator::ExternalForceActuator)
    (;signifier,operator) = actuator
    T = get_numbertype(signifier.body)
    zero(T)
end

function generalized_force(bot::Robot,actuator::ExternalForceActuator)
    (;structure) = bot
    (;state) = structure
    (;t) = state.system
    (;signifier,operator) = actuator
    transpose(C)*operator(t)
end

function action_jacobian(bot::Robot,actuator::ExternalForceActuator)
    (;structure) = bot
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
end

function action_jacobian(bot::Robot,q,q̇,λ,u)
    (;num_of_actions) = bot.hub.scheme
    ∂f∂u = zeros(T,nu)
    foreach(actuators) do actuator
        ∂f∂u .+= action_jacobian(bot,actuator)
    end
end

function actuate!(bot::Robot,actuator::ExternalForceActuator)
    (;id,signifier,operator) = actuator
    actuate!(signifier,operator)
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

function actuate!(bot::Robot,actuator::RegisterActuator)
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

function actuate!(actuator::SMAHeater,structure,u;inc=false,abs=true)
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
    actuate!(s,heating_law)
end

function actuate!(actuator::SMAHeater,heating_law)
    s.law.F0, s.law.k = heating_law(s.state.temp)
end
