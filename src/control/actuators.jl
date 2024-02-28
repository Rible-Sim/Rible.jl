abstract type AbstractOperator end
struct NonOperator <: AbstractOperator end
struct NaiveOperator <: AbstractOperator 
    num_of_actions::Int
end

abstract type AbstractActuator end

get_id(actuator::AbstractActuator) = actuator.id
get_numbertype(actuator::AbstractActuator) = get_numbertype(actuator.signifier)
get_num_of_actions(::AbstractActuator) = 0
function execute!(structure::Structure,actuator)
    (;id,signifier,operator,force) = actuator
    structure.state.system.F .+= generalized_force(structure,actuator)
end
function generalized_force_jacobian!(∂F∂u,structure::Structure,actuator)
end

struct ExternalForceActuator{sigType,operType,forceType,T} <: AbstractActuator
    id::Int
    signifier::sigType
    operator::operType
    force::forceType
    action::Vector{T}
end

get_num_of_actions(actuator::ExternalForceActuator) = length(actuator.action)

function get_actions(structure::Structure,actuator::ExternalForceActuator{sigType,<:NaiveOperator}) where {sigType}
    (;signifier,operator,) = actuator
    T = get_numbertype(structure)
    zeros(T,operator.num_of_actions)
end


# ExternalForceActuator 
function generalized_force(structure::Structure,actuator::ExternalForceActuator) 
    (;state) = structure
    (;signifier,operator,force) = actuator
    (;body,pid,) = signifier
    (;prop,coords) = body
    (;nmcs) = coords
    bid = body.prop.id
    (;q) = state.members[bid]
    c = to_local_coords(nmcs,prop.loci[pid].position)
    Tbody = build_T(structure,bid)
    C = to_position_jacobian(nmcs,q,c)*Tbody
    transpose(C)*force
end

function generalized_force_jacobian!(∂F∂u, structure::Structure,actuator::ExternalForceActuator)
    (;state) = structure
    (;signifier,operator,force) = actuator
    (;body,pid,) = signifier
    (;prop,coords) = body
    (;nmcs) = coords
    bid = body.prop.id
    (;q) = state.members[bid]
    c = to_local_coords(nmcs,prop.loci[pid].position)
    Tbody = build_T(structure,bid)
    C = to_position_jacobian(nmcs,q,c)*Tbody
    ∂F∂u .= transpose(C)*force
end

function execute!(structure::Structure,actuator::ExternalForceActuator{<:AbstractBody,<:NaiveOperator},u) 
    (;id,signifier,operator,force) = actuator
    structure.state.system.F .+= u.*generalized_force(structure,actuator)
end

function generalized_force(structure::Structure,actuator::ExternalForceActuator{<:AbstractBody,<:NaiveOperator},u) 
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
    u.*transpose(C)*force
end

function generalized_force_jacobian!(∂F∂u, structure::Structure,actuator::ExternalForceActuator{<:AbstractBody,<:NaiveOperator},u) 
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
    ∂F∂u .= transpose(C)*force
end

#DONE? 重要的地方
function execute!(structure::Structure,actuator::ExternalForceActuator{<:Apparatus{<:ClusterJoint},<:NaiveOperator},u) 
    (;id,signifier,operator,force) = actuator
    (;state) = structure
    (;t) = state.system
    segs = signifier.force
    #用驱动量u和joint的相关变量和参数， 计算驱动力， 加到系统中
    if id == 1
        segs[1].force.state.restlen = segs[1].force.original_restlen - u[1]
    end
    # INFO 1 不知道怎么计算驱动力，加到系统中，我先仿照之前的做法，也就是直接改变滑动绳索最末端的静止长度
    # 但是这个函数在计算force的时候调用（我也不知道为什么施加驱动要在构建force的时候执行），
    # 所以他每个计算步实际上要执行多次，这样不能简单地 -=。
end

function GravityActuator(id,body::AbstractBody)
    T = get_numbertype(body)
    signifier = body
    operator = NonOperator()
    action = T[]
    force = get_gravity(body)
    ExternalForceActuator(
        id,
        signifier,
        operator,
        force,
        action
    )
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
