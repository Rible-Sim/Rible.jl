abstract type AbstractOperator end
struct NonOperator <: AbstractOperator end
struct NaiveOperator <: AbstractOperator 
    num_of_actions::Int
end
struct FunctionOperator{F,JF,FV,JV} <: AbstractOperator
    func!::F
    jac!::JF
    func_vals::FV
    Jac_vals::JV
end

abstract type AbstractActuator end

get_id(actuator::AbstractActuator) = actuator.id
get_numbertype(actuator::AbstractActuator) = get_numbertype(actuator.signifier)


struct ExternalForceActuator{sigType,operType,forceType,T} <: AbstractActuator
    id::Int
    signifier::sigType
    operator::operType
    force::forceType
    action::Vector{T}
end


function measure(structure::AbstractStructure,actuator::ExternalForceActuator,u)
    transpose(u)*u ./2
end

function measure_gradient!(∂g∂u,structure::AbstractStructure,actuator::ExternalForceActuator,u)
    for i in eachindex(u)
        ∂g∂u[i] = u[i]
    end
end

function measure_hessians(structure::AbstractStructure,actuator::ExternalForceActuator,u)
    nu = get_num_of_actions(actuator)
    [I(nu)]
end

function measure_hessians!(∂²g∂q∂u, ∂²g∂q̇∂u, ∂²g∂u²,structure::AbstractStructure,actuator::ExternalForceActuator,u)
    ∂²g∂u²[:,:] .= I(get_num_of_actions(actuator))
end


"""
    get_num_of_actions(::AbstractActuator)

Return the number of actions available for the actuator.
"""
get_num_of_actions(::AbstractActuator) = 0

get_num_of_actions(actuator::ExternalForceActuator) = length(actuator.action)

"""
    get_initial_actions(structure::AbstractStructure,actuator)

Return a vector of actions on the structure by the actuator.
"""
function get_initial_actions(structure::AbstractStructure,actuator::ExternalForceActuator{sigType,<:NaiveOperator}) where {sigType}
    (;signifier,operator,) = actuator
    T = get_numbertype(structure)
    zeros(T,operator.num_of_actions)
end

"""
    execute!(structure::AbstractStructure,actuator,u)

Execute the actuator on the structure with the provided actions `u`.
"""
function execute!(structure::AbstractStructure,actuator::ExternalForceActuator{<:Signifier,<:NaiveOperator,<:AbstractMatrix},u) 
    rows, Clocal = _actuator_rows_and_jac(structure, actuator)
    local_view = @view structure.state.system.F[rows]
    Ct = transpose(Clocal)
    for j in axes(actuator.force, 2)
        mul!(local_view, Ct, @view(actuator.force[:, j]), u[j], one(eltype(local_view)))
    end
    structure
end


"""
    gen_force_actu_jacobian!(∂F∂u, structure::AbstractStructure,actuator,u)
 
Compute and store in `∂F∂u` the Jacobian of the generalized force with respect to the actuator's actions, with the provided actions `u`.
"""
function gen_force_actu_jacobian!(∂F∂u,structure::AbstractStructure,actuator::AbstractActuator,u)
end

@inline function _actuator_rows_and_jac(structure::AbstractStructure, actuator::ExternalForceActuator)
    (;state, connectivity) = structure
    (;signifier) = actuator
    (;body, pid,) = signifier
    bid = body.prop.id
    rows = connectivity.bodyid2sys_full_coords[bid]
    rows, body.cache.Cps[pid]
end


@inline function gen_force_actu_jacobian!(
        ∂F∂u::AbstractMatrix,
        structure::AbstractStructure,
        actuator::ExternalForceActuator{sigType, operType, forceType, T},
        u,
    ) where {sigType, operType, forceType<:AbstractMatrix, T}
    rows, Clocal = _actuator_rows_and_jac(structure, actuator)
    fill!(∂F∂u, zero(eltype(∂F∂u)))
    local_view = @view ∂F∂u[rows, :]
    mul!(local_view, transpose(Clocal), actuator.force)
    ∂F∂u
end
