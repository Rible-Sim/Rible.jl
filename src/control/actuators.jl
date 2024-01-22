
## function actuator_jacobian(bot::Robot,act::Actuator{<:Signifier,<:AngleCapta})

## end


abstract type AbstractCoupler end
struct Uncoupled <: AbstractCoupler end
struct Ganged <: AbstractCoupler end
struct Serial <: AbstractCoupler end

abstract type AbstractActuator end

get_numbertype(act::AbstractActuator) = get_numbertype(act.sig)

function get_actions(structure::Structure,act::AbstractActuator)
    (;sig,force) = act
    T = get_numbertype(sig.body)
    zero(T)
end

struct ExternalForceActuator{RT,forceType} <: AbstractActuator
    id::Int
    sig::RT
    force::forceType
end

get_num_of_actions(::ExternalForceActuator) = 1

function generalized_force(bot::Robot,act::ExternalForceActuator)
    (;structure) = bot
    (;state) = structure
    (;t) = state.system
    (;sig,force) = act
    transpose(C)*force(t)
end

function action_jacobian(bot::Robot,act::ExternalForceActuator)
    (;structure) = bot
    (;state) = structure
    (;t) = state.system
    (;sig,force) = act
    (;body,pid,) = sig
    (;prop,coords) = body
    (;nmcs) = coords
    bid = body.prop.id
    (;q) = state.members[bid]
    c = to_local_coords(nmcs,prop.loci[pid].position)
    Tbody = build_T(structure,bid)
    C = to_position_jacobian(nmcs,q,c)*Tbody
    transpose(C)*force(t)
end

function action_jacobian(bot::Robot,q,q̇,λ,u)
    (;num_of_actions) = bot.hub.scheme
    ∂f∂u = zeros(T,nu)
    foreach(actuators) do act
        ∂f∂u .+= action_jacobian(bot,act)
    end
end

struct RestLengthActuator{CT<:AbstractCoupler,RT} <: AbstractActuator
    id::Int
    coupler::CT
    apparatuses::RT
end

get_num_of_actions(::RestLengthActuator) = 1

function RestLengthActuator(actid,id::Int,value::Number)
    RestLengthActuator(actid,Serial(),(ids=[id],values=[value]))
end

function actuate!(bot::Robot,act::ExternalForceActuator)
    (;id,sig,force) = act
    actuate!(sig,force)
end

function actuate!(bot::Robot,act::RestLengthActuator{<:Uncoupled})
    (;sig) = act
    (;ids, values) = sig
    foreach(apparatuses) do apparatus
        apparatus.state.restlen = original_restlen + u[id]
    end
end

function actuate!(bot::Robot,act::RestLengthActuator{<:Serial})
    (;sig) = act
    (;ids, values) = sig
    foreach(apparatuses) do apparatus
        apparatus
        apparatus.state.restlen = original_restlen + u
    end
end

function actuate!(bot::Robot,act::RestLengthActuator{<:Ganged})
    (;sig) = act
    (;ids, values) = sig
    cable1 = (apparatuses,ids[1])
    cable2 = (apparatuses,ids[2])
    cable1.state.restlen = values[1] + u
    cable2.state.restlen = values[2] - u
end

function actuate!(bot::Robot,u::Nothing,t::Number) end

function actuate!(bot::Robot,u::AbstractVector,t::Number)
    (;actuators) = bot.hub
    foreach(actuators) do actuator
        actuate!(bot,actuator)
    end
end

struct SMAHeater{CT,HT,TT} <:AbstractActuator
    sig::CT
    heating_law::HT
    traj::TT
end
#
#
function SMAHeater(sig::Signifier,heating_law)
    SMAHeater(sig,heating_law,traj)
end
#
#
function actuate!(act::SMAHeater,structure,u;inc=false,abs=true)
    (;SMA_cables) = structure.apparatuses
    (;id_string, original_value, heating_law) = act
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

function actuate!(act::SMAHeater,heating_law)
    s.law.F0, s.law.k = heating_law(s.state.temp)
end
