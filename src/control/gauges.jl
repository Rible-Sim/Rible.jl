
abstract type AbstractCapta end
struct PositionCapta <: AbstractCapta end
struct VelocityCapta <: AbstractCapta end
struct PosVelCapta <: AbstractCapta end
struct AngleCapta <: AbstractCapta end

function measure(structure::Structure,signifier::Signifier,capta::PositionCapta)
    (;pid,body) = signifier
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
    (;body,pid) = signifier
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
    c = to_local_coords(nmcs,prop.loci[pid].position)
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
    c = to_local_coords(nmcs,prop.loci[pid].position)
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

function measure_jacobian!(∂e∂q,∂e∂q̇,structure::Structure,ref_err::ErrorGauge)
    (;signifier,capta,reference) = ref_err
    m = measure(structure,signifier,capta)
    Jq,Jq̇ = measure_jacobian(structure,signifier,capta)
    ∂e∂q .+= transpose(m)*Jq
    ∂e∂q̇ .+= transpose(m)*Jq̇
end


function get_errors(structure,gauges,coalition)
    (;num_of_errors, gaugeid2sys_errors_idx) = coalition.nt
    T = get_numbertype(structure)
    e = zeros(T,num_of_errors)
    foreach(gauges) do gauge
        err_idx = gaugeid2sys_errors_idx[gauge.id]
        e[err_idx] .= measure(structure,gauge)
    end
    e
end

function errors_jacobian!(∂e∂q,∂e∂q̇,bot::Robot)
    (;structure,hub) = bot
    (;gauges,coalition) = hub
    (;num_of_errors, gaugeid2sys_errors_idx,) = coalition.nt
    foreach(gauges) do gauge
        err_idx = gaugeid2sys_errors_idx[gauge.id]
        measure_jacobian!(
            (@view ∂e∂q[err_idx,:]),
            (@view ∂e∂q̇[err_idx,:]),
            structure,gauge
        )
    end
end

function errors_jacobian(bot::Robot)
    (;structure,hub) = bot
    (;num_of_free_coords) = structure.connectivity.indexed
    (;num_of_errors, gaugeid2sys_errors_idx,) = hub.coalition.nt
    T = get_numbertype(structure)
    ∂e∂q = zeros(T,num_of_errors,num_of_free_coords)
    ∂e∂q̇ = zeros(T,num_of_errors,num_of_free_coords)
    errors_jacobian!(∂e∂q,∂e∂q̇,bot)
    ∂e∂q, ∂e∂q̇
end