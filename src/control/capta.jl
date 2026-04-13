abstract type AbstractCaptum end
struct PositionCaptum <: AbstractCaptum end
struct VelocityCaptum <: AbstractCaptum end
struct PosVelCaptum <: AbstractCaptum end
struct AngularPositionCaptum <: AbstractCaptum end
struct AngularVelocityCaptum <: AbstractCaptum end
struct CoMPositionCaptum <: AbstractCaptum end
struct CoMVelocityCaptum <: AbstractCaptum end
struct CoMPosVelCaptum <: AbstractCaptum end

"""
    FullStateCaptum

A captum that measures the full system state [q; v] (generalized coordinates and velocities).
Unlike other Captum types which measure from a specific signifier (body + locus),
this captum measures system-level quantities.
"""
struct FullStateCaptum <: AbstractCaptum end

function measure(structure::AbstractStructure,signifier,captum::FullStateCaptum)
    (; q, q̇) = structure.state.system
    return ComponentArray((q = q, q̇ = q̇))
end

"""
    measure(structure::AbstractStructure,signifier::Signifier,captum)

Return the measurement of the captum on the signifier.
"""
function measure(structure::AbstractStructure,signifier::Signifier,captum::PositionCaptum)
    (;pid,body) = signifier
    (;prop,state,cache) = body
    state.loci_states[pid].frame.position
end

function measure(structure::AbstractStructure,signifier::Signifier,captum::VelocityCaptum)
    body = signifier.body
    (;pid) = signifier
    (;prop,state,cache) = body
    state.loci_states[pid].frame.velocity
end

function measure(structure::AbstractStructure,signifier::Signifier,captum::PosVelCaptum)
    (;body,pid) = signifier
    (;prop,state,cache) = body
    vcat(
        state.loci_states[pid].frame.position,
        state.loci_states[pid].frame.velocity
    )
end

"""
    measure(structure::AbstractStructure, signifier, captum::CoMPositionCaptum)

Measure center of mass position for constant mass system: r = B*q
where B maps generalized coordinates to center of mass position.
"""
function measure(structure::AbstractStructure, signifier, captum::CoMPositionCaptum)
    (; q) = structure.state.system
    B = compute_com_position_jacobian(structure)
    return B * q
end

"""
    measure(structure::AbstractStructure, signifier, captum::CoMVelocityCaptum)

Measure center of mass velocity for constant mass system: v = B*q̇
where B maps generalized velocities to center of mass velocity.
"""
function measure(structure::AbstractStructure, signifier, captum::CoMVelocityCaptum)
    (; q̇) = structure.state.system
    B = compute_com_position_jacobian(structure)
    return B * q̇
end

"""
    measure(structure::AbstractStructure, signifier, captum::CoMPosVelCaptum)

Measure center of mass position and velocity for constant mass system: [r; v] = [B*q; B*q̇]
"""
function measure(structure::AbstractStructure, signifier, captum::CoMPosVelCaptum)
    (; q, q̇) = structure.state.system
    B = compute_com_position_jacobian(structure)
    return vcat(B * q, B * q̇)
end

"""
    measure_jacobian(structure::AbstractStructure,signifier::Signifier,captum)

Return the Jacobian of the captum's measurement w.r.t. the signifier's coordinates and velocities.
"""
function measure_jacobian!(Jq, Jq̇, Js, structure::AbstractStructure, signifier::Signifier, captum::PositionCaptum)
    body = signifier.body
    (;pid) = signifier
    (;prop,coords) = body
    q = structure.state.members[body.prop.id].q
    c = to_local_coords(coords,prop.loci[pid].position)
    idx = structure.connectivity.bodyid2sys_full_coords[body.prop.id]
    C = to_position_jacobian(coords,q,c)
    n = size(C,1)
    fill!(Jq, zero(eltype(Jq)))
    fill!(Jq̇, zero(eltype(Jq̇)))
    fill!(Js, zero(eltype(Js)))
    copyto!(@view(Jq[1:n, idx]), C)
    return Jq,Jq̇
end

function measure_jacobian!(Jq, Jq̇, Js, structure::AbstractStructure, signifier::Signifier, captum::VelocityCaptum)
    body = signifier.body
    (;pid) = signifier
    (;prop,coords) = body
    q = structure.state.members[body.prop.id].q
    q̇ = structure.state.members[body.prop.id].q̇
    c = to_local_coords(coords,prop.loci[pid].position)
    idx = structure.connectivity.bodyid2sys_full_coords[body.prop.id]
    C = to_position_jacobian(coords,q,c)
    ∂Cq̇∂q = to_velocity_jacobian(coords,q,q̇,c)
    fill!(Jq, zero(eltype(Jq)))
    fill!(Jq̇, zero(eltype(Jq̇)))
    fill!(Js, zero(eltype(Js)))
    copyto!(@view(Jq[1:size(∂Cq̇∂q,1), idx]), ∂Cq̇∂q)
    copyto!(@view(Jq̇[1:size(C,1), idx]), C)
    return Jq,Jq̇
end

function measure_jacobian!(Jq, Jq̇, Js, structure::AbstractStructure, signifier::Signifier, captum::PosVelCaptum)
    body = signifier.body
    (;pid) = signifier
    (;prop,coords) = body
    q = structure.state.members[body.prop.id].q
    q̇ = structure.state.members[body.prop.id].q̇
    c = to_local_coords(coords,prop.loci[pid].position)
    idx = structure.connectivity.bodyid2sys_full_coords[body.prop.id]
    C = to_position_jacobian(coords,q,c)
    ∂Cq̇∂q = to_velocity_jacobian(coords,q,q̇,c)
    n = size(C,1)
    fill!(Jq, zero(eltype(Jq)))
    fill!(Jq̇, zero(eltype(Jq̇)))
    fill!(Js, zero(eltype(Js)))
    copyto!(@view(Jq[1:n, idx]), C)
    copyto!(@view(Jq[n+1:2n, idx]), ∂Cq̇∂q)
    copyto!(@view(Jq̇[n+1:2n, idx]), C)
    return Jq,Jq̇
end

function measure(structure::AbstractStructure, signifier::Hen2Egg, captum::AbstractCaptum)
    (;hen, egg) = signifier
    return measure(structure, egg, captum) .- measure(structure, hen, captum)
end

function measure_jacobian!(Jq, Jq̇, Js, structure::AbstractStructure, signifier::Hen2Egg, captum::AbstractCaptum)
    (;hen, egg) = signifier
    measure_jacobian!(Jq, Jq̇, Js, structure, egg, captum)
    Jq_hen = similar(Jq)
    Jq̇_hen = similar(Jq̇)
    Js_hen = similar(Js)
    measure_jacobian!(Jq_hen, Jq̇_hen, Js_hen, structure, hen, captum)
    Jq .-= Jq_hen
    Jq̇ .-= Jq̇_hen
    Js .-= Js_hen
    return Jq, Jq̇, Js
end

function measure_hessians(structure::AbstractStructure, signifier::Hen2Egg, captum::AbstractCaptum)
    (;hen, egg) = signifier
    Hqq_egg, Hq̇q̇_egg, Hq̇q_egg = measure_hessians(structure, egg, captum)
    Hqq_hen, Hq̇q̇_hen, Hq̇q_hen = measure_hessians(structure, hen, captum)
    Hqq = Hqq_egg .- Hqq_hen
    Hq̇q̇ = Hq̇q̇_egg .- Hq̇q̇_hen
    Hq̇q = Hq̇q_egg .- Hq̇q_hen
    return Hqq, Hq̇q̇, Hq̇q
end

function measure_jacobian!(Jq, Jq̇, Js, structure::AbstractStructure, signifier, captum::CoMPositionCaptum)
    B = compute_com_position_jacobian(structure)
    fill!(Jq, zero(eltype(Jq)))
    fill!(Jq̇, zero(eltype(Jq̇)))
    fill!(Js, zero(eltype(Js)))
    copyto!(@view(Jq[1:size(B,1), 1:size(B,2)]), B)
    return Jq,Jq̇
end

function measure_jacobian!(Jq, Jq̇, Js, structure::AbstractStructure, signifier, captum::CoMVelocityCaptum)
    B = compute_com_position_jacobian(structure)
    fill!(Jq, zero(eltype(Jq)))
    fill!(Jq̇, zero(eltype(Jq̇)))
    fill!(Js, zero(eltype(Js)))
    copyto!(@view(Jq̇[1:size(B,1), 1:size(B,2)]), B)
    return Jq,Jq̇
end

function measure_jacobian!(Jq, Jq̇, Js, structure::AbstractStructure, signifier, captum::CoMPosVelCaptum)
    B = compute_com_position_jacobian(structure)
    m = size(B,1)
    fill!(Jq, zero(eltype(Jq)))
    fill!(Jq̇, zero(eltype(Jq̇)))
    fill!(Js, zero(eltype(Js)))
    copyto!(@view(Jq[1:m, 1:size(B,2)]), B)
    copyto!(@view(Jq̇[m+1:2m, 1:size(B,2)]), B)
    return Jq,Jq̇
end

function measure_jacobian!(Jq, Jq̇, Js, structure::AbstractStructure, signifier, captum::FullStateCaptum)
    nq = structure.connectivity.num_of_full_coords
    fill!(Jq, zero(eltype(Jq)))
    fill!(Jq̇, zero(eltype(Jq̇)))
    fill!(Js, zero(eltype(Js)))
    for i in 1:nq
        Jq[i,i] = one(eltype(Jq))
        Jq̇[nq+i,i] = one(eltype(Jq̇))
    end
    return Jq,Jq̇
end

function measure_jacobian!(Jq, Jq̇, Js, structure::AbstractStructure, signifier, captum)
    throw(MethodError(measure_jacobian!, (Jq, Jq̇, Js, structure, signifier, captum)))
end

"""
    measure_hessians(structure::AbstractStructure,signifier::Signifier,captum)

Return the Hessian of the captum's measurement w.r.t. the signifier's coordinates and velocities.
"""
function measure_hessians(structure::AbstractStructure,signifier::Signifier,captum::PositionCaptum)
    body = signifier.body
    (;pid) = signifier
    (;prop,state,cache,coords) = body
    q = structure.state.members[body.prop.id].q
    npos = prop.loci[pid].position |> length
    nq = structure.connectivity.num_of_full_coords
    T = get_numbertype(structure)
    Hqq = [spzeros(T,nq,nq) for _ in 1:npos]
    Hq̇q̇ = [spzeros(T,nq,nq) for _ in 1:npos]
    Hq̇q = [spzeros(T,nq,nq) for _ in 1:npos]
    Hqq, Hq̇q̇, Hq̇q
end

function measure_hessians(structure::AbstractStructure,signifier::Signifier,captum::VelocityCaptum)
    body = signifier.body
    (;pid) = signifier
    (;prop,state,cache,coords) = body
    q = structure.state.members[body.prop.id].q
    q̇ = structure.state.members[body.prop.id].q̇
    c = to_local_coords(coords,prop.loci[pid].position)
    npos = prop.loci[pid].position |> length
    nq = structure.connectivity.num_of_full_coords
    # ∂Cq̇∂q∂q = to_velocity_hessian(coords,q,q̇,c)*Tbody
    T = get_numbertype(structure)
    Hqq = [spzeros(T,nq,nq) for _ in 1:npos]
    Hq̇q̇ = [spzeros(T,nq,nq) for _ in 1:npos]
    Hq̇q = [spzeros(T,nq,nq) for _ in 1:npos]
    Hqq, Hq̇q̇, Hq̇q
end

function measure_hessians(structure::AbstractStructure,signifier::Signifier,captum::PosVelCaptum)
    body = signifier.body
    (;pid) = signifier
    (;prop,state,cache,coords) = body
    q = structure.state.members[body.prop.id].q
    q̇ = structure.state.members[body.prop.id].q̇
    c = to_local_coords(coords,prop.loci[pid].position)
    Tbody = build_T(structure,body.prop.id)
    C = to_position_jacobian(coords,q,c)*Tbody
    ∂Cq̇∂q = to_velocity_jacobian(coords,q,q̇,c)*Tbody
    T = get_numbertype(structure)
    npos = prop.loci[pid].position |> length
    nq = structure.connectivity.num_of_full_coords
    # ∂Cq̇∂q∂q = to_velocity_hessian(coords,q,q̇,c)*Tbody
    Hqq = [spzeros(T,nq,nq) for _ in 1:2npos]
    Hq̇q̇ = [spzeros(T,nq,nq) for _ in 1:2npos]
    Hq̇q = [spzeros(T,nq,nq) for _ in 1:2npos]
    Hqq, Hq̇q̇, Hq̇q
end

"""
    measure_hessians(structure::AbstractStructure, signifier, captum::CoMPositionCaptum)

Return Hessians for center of mass position (all zeros since r = B*q is linear).
"""
function measure_hessians(structure::AbstractStructure, signifier, captum::CoMPositionCaptum)
    nq = structure.connectivity.num_of_full_coords
    T = get_numbertype(structure)
    ndim = 3  # 3D center of mass position
    Hqq = [spzeros(T, nq, nq) for _ in 1:ndim]
    Hq̇q̇ = [spzeros(T, nq, nq) for _ in 1:ndim]
    Hq̇q = [spzeros(T, nq, nq) for _ in 1:ndim]
    return Hqq, Hq̇q̇, Hq̇q
end

"""
    measure_hessians(structure::AbstractStructure, signifier, captum::CoMVelocityCaptum)

Return Hessians for center of mass velocity (all zeros since v = B*q̇ is linear).
"""
function measure_hessians(structure::AbstractStructure, signifier, captum::CoMVelocityCaptum)
    nq = structure.connectivity.num_of_full_coords
    T = get_numbertype(structure)
    ndim = 3  # 3D center of mass velocity
    Hqq = [spzeros(T, nq, nq) for _ in 1:ndim]
    Hq̇q̇ = [spzeros(T, nq, nq) for _ in 1:ndim]
    Hq̇q = [spzeros(T, nq, nq) for _ in 1:ndim]
    return Hqq, Hq̇q̇, Hq̇q
end

"""
    measure_hessians(structure::AbstractStructure, signifier, captum::CoMPosVelCaptum)

Return Hessians for center of mass position and velocity (all zeros since both are linear).
"""
function measure_hessians(structure::AbstractStructure, signifier, captum::CoMPosVelCaptum)
    nq = structure.connectivity.num_of_full_coords
    T = get_numbertype(structure)
    ndim = 6  # 3D position + 3D velocity
    Hqq = [spzeros(T, nq, nq) for _ in 1:ndim]
    Hq̇q̇ = [spzeros(T, nq, nq) for _ in 1:ndim]
    Hq̇q = [spzeros(T, nq, nq) for _ in 1:ndim]
    return Hqq, Hq̇q̇, Hq̇q
end


"""
    compute_com_position_jacobian(structure::AbstractStructure)

Compute the Jacobian matrix B that maps generalized coordinates to center of mass position.
For constant mass system: r_com = B*q, where B = (1/M_total) * Σ m_i * J_i
J_i is the position Jacobian of body i's center of mass.
"""
function compute_com_position_jacobian(structure::AbstractStructure)
    (;bodies) = structure
    (;bodyid2sys_full_coords,num_of_full_coords) = structure.connectivity
    
    T = get_numbertype(structure)
    N = get_num_of_dims(structure)
    
    # Initialize B matrix (N x num_of_full_coords)
    B = spzeros(T, N, num_of_full_coords)
    total_mass = Ref(zero(T))
    
    # Sum over all bodies: mass-weighted position Jacobians
    foreach(bodies) do body
        (; coords) = body
        (; id, mass, mass_locus) = body.prop
        mass_center = mass_locus.position
        total_mass[] += mass
        # Convert mass center to local coordinates
        cₘ = to_local_coords(coords, mass_center)
        # Compute position Jacobian of body's center of mass
        Cₘ = to_position_jacobian(coords, cₘ)
        # Add mass-weighted contribution to system coordinates
        B[:, bodyid2sys_full_coords[id]] .+= mass .* Cₘ
    end
    
    # Normalize by total mass
    if total_mass[] > zero(T)
        B ./= total_mass[]
    end
    
    return B
end
