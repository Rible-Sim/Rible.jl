
abstract type AbstractCapta end
struct PositionCapta <: AbstractCapta end
struct VelocityCapta <: AbstractCapta end
struct PosVelCapta <: AbstractCapta end
struct AngleCapta <: AbstractCapta end

function measure(structure::Structure,signifier::Signifier,capta::PositionCapta)
    (;signifier,) = gauge
    body = signifier.body
    (;pid) = signifier
    (;prop,state,cache) = body
    state.loci_states[pid].frame.position
end

function measure(structure::Structure,signifier::Signifier,capta::VelocityCapta)
    body = signifier.body
    (;pid) = signifier
    (;prop,state,cache) = body
    state.loci_states[pid].frame.velocity
end

function measure(structure::Structure,signifier::Signifier,capta::PosVelCapta)
    body = signifier.body
    (;pid) = signifier
    (;prop,state,cache) = body
    vcat(
        state.loci_states[pid].frame.position,
        state.loci_states[pid].frame.velocity
    )
end

function measure_jacobian(structure::Structure,signifier::Signifier,capta::PositionCapta)
    body = signifier.body
    (;pid) = signifier
    (;prop,state,cache,coords) = body
    (;nmcs) = coords
    q = structure.state.members[body.prop.id].q
    c = to_local_coords(nmcs,prop.loci[pid].frame.position)
    Tbody = build_T(structure,body.prop.id)
    C = to_position_jacobian(nmcs,q,c)*Tbody
    C,zero(C)
end

function measure_jacobian(structure::Structure,signifier::Signifier,capta::VelocityCapta)
    body = signifier.body
    (;pid) = signifier
    (;prop,state,cache,coords) = body
    (;nmcs) = coords
    q = structure.state.members[body.prop.id].q
    q̇ = structure.state.members[body.prop.id].q̇
    c = to_local_coords(nmcs,q,prop.loci[pid].frame.position)
    Tbody = build_T(structure,body.prop.id)
    C = to_position_jacobian(nmcs,q,c)*Tbody
    ∂Cq̇∂q = to_velocity_jacobian(nmcs,q,q̇,c)*Tbody
    ∂Cq̇∂q, C
    
end

function measure_jacobian(structure::Structure,signifier::Signifier,capta::PosVelCapta)
    body = signifier.body
    (;pid) = signifier
    (;prop,state,cache,coords) = body
    (;nmcs) = coords
    q = structure.state.members[body.prop.id].q
    q̇ = structure.state.members[body.prop.id].q̇
    c = to_local_coords(nmcs,prop.loci[pid].position)
    Tbody = build_T(structure,body.prop.id)
    C = to_position_jacobian(nmcs,q,c)*Tbody
    ∂Cq̇∂q = to_velocity_jacobian(nmcs,q,q̇,c)*Tbody
    vcat(
        C,
        ∂Cq̇∂q
    ),
    vcat(
        zero(C),
        C
    )
end

function measure(structure::Structure,signifier::Signifier,capta::AngleCapta)
    appar = signifier.body
    (;
        hen2egg,
        cache,
    ) = appar.joint
    (;
        relative_core
    ) = cache
    (;hen,egg) = hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    q_hen = structure.state.members[id_hen].q
    q_egg = structure.state.members[id_egg].q
    q_jointed = vcat(
        q_hen,
        q_egg
    )
    jointed2angles = make_jointed2angles(hen2egg,relative_core)
    angles = jointed2angles(q_jointed)
    angles[signifier.pid]
end

function measure_jacobian(structure::Structure,signifier::Signifier,capta::AngleCapta)
    appar = signifier.body
    (;
        hen2egg,
        cache,
    ) = appar.joint
    (;
        relative_core
    ) = cache
    (;hen,egg) = hen2egg
    id_hen = hen.body.prop.id
    id_egg = egg.body.prop.id
    q_hen = structure.state.members[id_hen].q
    q_egg = structure.state.members[id_egg].q
    q_jointed = vcat(
        q_hen,
        q_egg
    )
    jointed2angles = make_jointed2angles(hen2egg,relative_core)
    angles_jacobian = ForwardDiff.jacobian(jointed2angles,q_jointed)
    angles_jacobian[signifier.pid,:]
end

abstract type AbstractGauge end

struct CaptaGauge{sigType,captaType} <: AbstractGauge
    id::Int
    signifier::sigType
    capta::captaType
end

measure(structure::Structure,gauge::CaptaGauge) = measure(structure,gauge.signifier,gauge.capta)

struct ErrorGauge{sigType,captaType,referenceType} <: AbstractGauge
    id::Int
    signifier::sigType
    capta::captaType
    reference::referenceType
end

get_numbertype(gauge::AbstractGauge) = get_numbertype(gauge.signifier.body)
get_num_of_errors(::CaptaGauge) = 0
get_num_of_errors(::ErrorGauge) = 1

function Base.isless(a::AbstractGauge,b::AbstractGauge)
    isless(a.id,b.id)
end

get_id(gauge::AbstractGauge) = gauge.id

function Base.isless(a::ErrorGauge,b::ErrorGauge)
    isless(a.id,b.id)
end

struct ErrorCost{ref_errorsType}
    ref_errors::ref_errorsType
end

function measure(structure::Structure,ref_err::ErrorGauge)
    (;signifier,capta,reference) = ref_err
    1/2*sum((measure(structure,signifier,capta).-reference).^2)
end

function measure_jacobian!(∂ϕ∂qᵀ,∂ϕ∂q̇ᵀ,structure::Structure,ref_err::ErrorGauge)
    (;gauge,reference) = ref_err
    m = measure(structure,gauge)
    Jq,Jq̇ = measure_jacobian(structure,gauge)
    ∂ϕ∂qᵀ .+= transpose(transpose(m)*Jq)
    ∂ϕ∂q̇ᵀ .+= transpose(transpose(m)*Jq̇)
end
