
"""
$(TYPEDEF)

Actuator that operates on a register of values, such as the original rest lengths of cables.

# Fields
- `id::Int`: Unique identifier for the actuator
- `signifier::sigType`: The signifier type, typically representing what the actuator acts upon
- `operator::opType`: The operator type, defining how the actuator operates
- `register::regType`: The register containing values to be manipulated

# Methods
- `get_num_of_actions(actuator::RegisterActuator)`: Returns the number of actions available
- `get_initial_actions(structure::AbstractStructure, actuator::RegisterActuator)`: Returns a vector of zeros representing possible actions
- `execute!(structure::AbstractStructure, actuator::RegisterActuator, u)`: Executes the actuator's action on the structure

This actuator type is designed to work with a register of values, allowing for manipulation of multiple elements simultaneously based on the provided actions.
"""
struct RegisterActuator{sigType,opType,regType} <: AbstractActuator
    id::Int
    signifier::sigType
    operator::opType
    register::regType
end

get_num_of_actions(actuator::RegisterActuator) = length(actuator.register.values)

function get_initial_actions(structure::AbstractStructure,actuator::RegisterActuator)
    (;signifier,operator,register) = actuator
    zeros(eltype(register.values),get_num_of_actions(actuator))
end

function measure(structure::AbstractStructure,actuator::Rible.RegisterActuator,u::AbstractVector{T}) where T
    zero(T)
end

function measure_gradient!(∂g∂u,structure::AbstractStructure,actuator::Rible.RegisterActuator,u::AbstractVector{T}) where T
end

function measure_hessians!(∂²g∂q∂u, ∂²g∂q̇∂u, ∂²g∂u²,structure::AbstractStructure,actuator::RegisterActuator,u)
    #todo implement this for register actuators
    ∂²g∂u²[:,:] .= 0.0
end

function execute!(structure::AbstractStructure,actuator::RegisterActuator{<:Apparatus{<:PrototypeJoint,<:RheonomicJointForce},<:FunctionOperator},u)
    (;signifier,operator,register) = actuator
    (;func!, func_vals) = operator
    func!(func_vals, register.values, u)
    signifier.joint.violations .= func_vals
end
