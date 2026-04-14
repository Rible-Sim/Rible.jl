
abstract type AbstractGauge end

struct CaptumGauge{sigType,captumType} <: AbstractGauge
    id::Int
    signifier::sigType
    captum::captumType
end

struct ErrorGauge{sigType,captumType,referenceType} <: AbstractGauge
    id::Int
    signifier::sigType
    captum::captumType
    reference::referenceType
end

function get_numbertype(gauge::AbstractGauge) 
    error("get_numbertype not implemented for $(typeof(gauge))")
    get_numbertype(gauge.signifier.body)
end

get_numbertype(gauge::CaptumGauge{<:Hen2Egg}) = get_numbertype(gauge.signifier.hen.body)

get_num_of_capta(gauge::CaptumGauge{<:Signifier,<:PositionCaptum}) = get_num_of_dims(gauge.signifier.body)
get_num_of_capta(gauge::CaptumGauge{<:Signifier,<:VelocityCaptum}) = get_num_of_dims(gauge.signifier.body)
get_num_of_capta(gauge::CaptumGauge{<:Signifier,<:PosVelCaptum}) = 2*get_num_of_dims(gauge.signifier.body)

get_num_of_capta(gauge::CaptumGauge{<:Hen2Egg,<:PositionCaptum}) = get_num_of_dims(gauge.signifier.hen.body)
get_num_of_capta(gauge::CaptumGauge{<:Hen2Egg,<:VelocityCaptum}) = get_num_of_dims(gauge.signifier.hen.body)
get_num_of_capta(gauge::CaptumGauge{<:Hen2Egg,<:PosVelCaptum}) = 2*get_num_of_dims(gauge.signifier.hen.body)

get_num_of_capta(gauge::CaptumGauge{<:Apparatus,<:AngularPositionCaptum}) = 1
get_num_of_capta(gauge::CaptumGauge{<:Hen2Egg,<:AngularPositionCaptum}) = 1
get_num_of_capta(gauge::CaptumGauge{<:Hen2Egg,<:AngularVelocityCaptum}) = 1

get_num_of_capta(gauge::CaptumGauge{<:AbstractStructure,<:FullStateCaptum}) = 2*get_num_of_full_coords(gauge.signifier)


get_numbertype(gauge::ErrorGauge) = get_numbertype(gauge.signifier.body)
get_numbertype(gauge::ErrorGauge{<:Hen2Egg}) = get_numbertype(gauge.signifier.hen.body)

get_num_of_capta(gauge::ErrorGauge) = get_num_of_capta(CaptumGauge(gauge.id, gauge.signifier, gauge.captum))


function Base.isless(a::AbstractGauge,b::AbstractGauge)
    isless(a.id,b.id)
end

get_id(gauge::AbstractGauge) = gauge.id

struct ErrorCost{ref_errorsType}
    ref_errors::ref_errorsType
end

"""
    measure(structure::AbstractStructure,gauge::CaptumGauge)

Return the measurement of the captum gauge.
"""
measure(structure::AbstractStructure,gauge::CaptumGauge) = measure(structure,gauge.signifier,gauge.captum)


function measure!(out, structure::AbstractStructure, signifier::Signifier, captum::PositionCaptum)
    (;body, pid) = signifier
    (;state) = body
    copyto!(out, state.loci_states[pid].frame.position)
end

function measure!(out, structure::AbstractStructure, signifier::Signifier, captum::VelocityCaptum)
    (;body, pid) = signifier
    (;state) = body
    copyto!(out, state.loci_states[pid].frame.velocity)
end

function measure!(out, structure::AbstractStructure, signifier::Signifier, captum::PosVelCaptum)
    (;body, pid) = signifier
    (;state) = body
    n = length(state.loci_states[pid].frame.position)
    copyto!(@view(out[1:n]), state.loci_states[pid].frame.position)
    copyto!(@view(out[n+1:2n]), state.loci_states[pid].frame.velocity)
end

function measure!(out, structure::AbstractStructure, signifier, captum::CoMPositionCaptum)
    (; q) = structure.state.system
    B = compute_com_position_jacobian(structure)
    mul!(@view(out[1:size(B, 1)]), B, q)
end

function measure!(out, structure::AbstractStructure, signifier, captum::CoMVelocityCaptum)
    (; q̇) = structure.state.system
    B = compute_com_position_jacobian(structure)
    mul!(@view(out[1:size(B, 1)]), B, q̇)
end

function measure!(out, structure::AbstractStructure, signifier, captum::CoMPosVelCaptum)
    (; q, q̇) = structure.state.system
    B = compute_com_position_jacobian(structure)
    m = size(B, 1)
    mul!(@view(out[1:m]), B, q)
    mul!(@view(out[m+1:2m]), B, q̇)
end

function measure!(out, structure::AbstractStructure, signifier::Nothing, captum::FullStateCaptum)
    (; q, q̇) = structure.state.system
    nq = length(q)
    copyto!(@view(out[1:nq]), q)
    copyto!(@view(out[nq+1:2nq]), q̇)
end

function measure!(out, structure::AbstractStructure, signifier::Hen2Egg, captum)
    (;hen, egg) = signifier
    measure!(out, structure, egg, captum)
    out .-= measure(structure, hen, captum)
end

"""
    measure(structure::AbstractStructure,err_gau::ErrorGauge{sigType,captumType,<:AbstractVector}) where {sigType,captumType}

Return the measurement of the error as deviation of the captum from the reference.
"""
function measure(structure::AbstractStructure,err_gau::ErrorGauge{sigType,captumType,<:AbstractVector}) where {sigType,captumType}
    (;signifier,captum,reference) = err_gau
    @debug "measure(structure::AbstractStructure,err_gau::ErrorGauge)" measure(structure,signifier,captum) reference
    d = measure(structure,signifier,captum).-reference
    e = 1/2*sum(d.^2)
end

function measure(structure::AbstractStructure,err_gau::ErrorGauge{sigType,captumType,<:Function}) where {sigType,captumType}
    (;signifier,captum,reference) = err_gau
    t = structure.state.system.t
    @debug "measure(structure::AbstractStructure,err_gau::ErrorGauge)" measure(structure,signifier,captum) reference
    d = measure(structure,signifier,captum).-reference(t)
    e = 1/2*sum(d.^2)
end

"""
    measure_jacobian!(∂g∂q,∂g∂q̇,∂g∂s,structure::AbstractStructure,err_gau::ErrorGauge)

Compute the gradient of the error function w.r.t. the coordinates and velocities of the structure.
"""
function measure_gradient!(∂g∂q,∂g∂q̇,∂g∂s,structure::AbstractStructure,err_gau::ErrorGauge)
    (;signifier,captum,reference) = err_gau
    len = length(reference)
    d = zeros(eltype(reference), len)
    measure!(d, structure, signifier, captum)
    d .-= reference
    nq = get_num_of_full_coords(structure)
    ns = get_num_of_aux_var(structure)
    Jq = zeros(eltype(d), len, nq)
    Jq̇ = similar(Jq)
    Js = zeros(eltype(d), len, ns)
    measure_jacobian!(Jq, Jq̇, Js, structure, signifier, captum)
    @debug "measure_gradient!(∂g∂q,∂g∂q̇,structure::AbstractStructure,err_gau::ErrorGauge)" d Jq Jq̇
    ∂g∂q .+= transpose(d)*Jq
    ∂g∂q̇ .+= transpose(d)*Jq̇
    ∂g∂s .+= transpose(d)*Js
end

function measure_gradient!(∂g∂q,∂g∂q̇,∂g∂s,structure::AbstractStructure,err_gau::ErrorGauge{sigType,captumType,<:Function}) where {sigType,captumType}
    (;signifier,captum,reference) = err_gau
    t = structure.state.system.t
    ref = reference(t)
    len = length(ref)
    d = zeros(eltype(ref), len)
    measure!(d, structure, signifier, captum)
    d .-= ref
    nq = get_num_of_full_coords(structure)
    Jq = zeros(eltype(d), len, nq)
    Jq̇ = similar(Jq)
    Js = zeros(eltype(d), len, ns)
    measure_jacobian!(Jq, Jq̇, Js, structure, signifier, captum)
    @debug "measure_gradient!(∂g∂q,∂g∂q̇,structure::AbstractStructure,err_gau::ErrorGauge)" d Jq Jq̇
    ∂g∂q .+= transpose(d)*Jq
    ∂g∂q̇ .+= transpose(d)*Jq̇
    ∂g∂s .+= transpose(d)*Js
end

function measure_gradient!(∂g∂q,∂g∂q̇,∂g∂s,structure::AbstractStructure,gauge::CaptumGauge)
    (;signifier,captum,) = gauge
    len = get_num_of_capta(gauge)
    d = zeros(eltype(structure.state.system.q), len)
    measure!(d, structure, signifier, captum)
    nq = get_num_of_full_coords(structure)
    ns = get_num_of_aux_var(structure)
    Jq = zeros(eltype(d), len, nq)
    Jq̇ = similar(Jq)
    Js = zeros(eltype(d), len, ns)
    measure_jacobian!(Jq, Jq̇, Js, structure, signifier, captum)
    @debug "measure_gradient!(∂g∂q,∂g∂q̇,structure::AbstractStructure,gauge::CaptumGauge)" d Jq Jq̇
    ∂g∂q .+= transpose(d)*Jq
    ∂g∂q̇ .+= transpose(d)*Jq̇
    ∂g∂s .+= transpose(d)*Js
end

function measure_gradient!(∂g∂q,∂g∂q̇,∂g∂s,structure::AbstractStructure,err_gau::ErrorGauge, gw)
    (;signifier,captum,reference) = err_gau
    len = length(reference)
    d = gw.d
    Jq_view = gw.Jq
    Jq̇_view = gw.Jq̇
    Js_view = gw.Js
    tmp_grad = gw.tmp_grad
    tmp_grad_s = gw.tmp_grad_s
    measure!(d, structure, signifier, captum)
    measure_jacobian!(Jq_view, Jq̇_view, Js_view, structure, signifier, captum)
    # @show err_gau.id d reference Js_view
    d .-= reference
    T = eltype(d)
    mul!(tmp_grad, transpose(Jq_view), d, one(T), zero(T))
    for i in eachindex(tmp_grad)
        ∂g∂q[i] += tmp_grad[i]
    end
    mul!(tmp_grad, transpose(Jq̇_view), d, one(T), zero(T))
    for i in eachindex(tmp_grad)
        ∂g∂q̇[i] += tmp_grad[i]
    end
    
    mul!(tmp_grad_s, transpose(Js_view), d, one(T), zero(T))
    for i in eachindex(tmp_grad_s)
        ∂g∂s[i] += tmp_grad_s[i]
    end
end

function measure_gradient!(∂g∂q,∂g∂q̇,∂g∂s,structure::AbstractStructure,err_gau::ErrorGauge{sigType,captumType,<:Function}, gw) where {sigType,captumType}
    (;signifier,captum,reference) = err_gau
    t = structure.state.system.t
    reference_val = reference(t)
    len = length(reference_val)
    d = gw.d
    measure!(d, structure, signifier, captum)
    d .-= reference_val
    Jq_view = gw.Jq
    Jq̇_view = gw.Jq̇
    Js_view = gw.Js
    tmp_grad = gw.tmp_grad
    tmp_grad_s = gw.tmp_grad_s
    measure_jacobian!(Jq_view, Jq̇_view, Js_view, structure, signifier, captum)
    T = eltype(d)
    mul!(tmp_grad, transpose(Jq_view), d, one(T), zero(T))
    for i in eachindex(tmp_grad)
        ∂g∂q[i] += tmp_grad[i]
    end
    mul!(tmp_grad, transpose(Jq̇_view), d, one(T), zero(T))
    for i in eachindex(tmp_grad)
        ∂g∂q̇[i] += tmp_grad[i]
    end
    
    mul!(tmp_grad_s, transpose(Js_view), d, one(T), zero(T))
    for i in eachindex(tmp_grad_s)
        ∂g∂s[i] += tmp_grad_s[i]
    end
end

function measure_gradient!(∂g∂q,∂g∂q̇,∂g∂s,structure::AbstractStructure,gauge::CaptumGauge, gw)
    len = get_num_of_capta(gauge)
    d = gw.d
    measure!(d, structure, gauge.signifier, gauge.captum)
    Jq_view = gw.Jq
    Jq̇_view = gw.Jq̇
    Js_view = gw.Js
    tmp_grad = gw.tmp_grad
    measure_jacobian!(Jq_view, Jq̇_view, Js_view, structure, gauge.signifier, gauge.captum)
    T = eltype(d)
    mul!(tmp_grad, transpose(Jq_view), d, one(T), zero(T))
    for i in eachindex(tmp_grad)
        ∂g∂q[i] += tmp_grad[i]
    end
    mul!(tmp_grad, transpose(Jq̇_view), d, one(T), zero(T))
    for i in eachindex(tmp_grad)
        ∂g∂q̇[i] += tmp_grad[i]
    end
    mul!(tmp_grad_s, transpose(Js_view), d, one(T), zero(T))
    for i in eachindex(tmp_grad_s)
        ∂g∂s[i] += tmp_grad_s[i]
    end
end

"""
    measure_hessians!(∂²g∂q²,∂²g∂q̇²,∂²g∂q̇∂q,structure::AbstractStructure,err_gau::ErrorGauge)

Compute the Hessian of the error function w.r.t. the coordinates and velocities of the structure.
"""
function measure_hessians!(∂²g∂q²,∂²g∂q̇²,∂²g∂q̇∂q,structure::AbstractStructure,err_gau::ErrorGauge)
    (;signifier,captum,reference) = err_gau
    len = length(reference)
    d = zeros(eltype(reference), len)
    measure!(d, structure, signifier, captum)
    d .-= reference
    nq = get_num_of_full_coords(structure)
    ns = get_num_of_aux_var(structure)
    Jq = zeros(eltype(d), len, nq)
    Jq̇ = similar(Jq)
    Js = zeros(eltype(d), len, ns)
    measure_jacobian!(Jq, Jq̇, Js, structure, signifier, captum)
    Hqq, Hq̇q̇, Hq̇q = measure_hessians(structure,signifier,captum)
    ∂²g∂q²[:,:]  .+= sum(d[i] * Hqq[i] for i in eachindex(d)) .+ transpose(Jq)*Jq
    ∂²g∂q̇²[:,:]  .+= sum(d[i] * Hq̇q̇[i] for i in eachindex(d)) .+ transpose(Jq̇)*Jq̇
    ∂²g∂q̇∂q[:,:] .+= sum(d[i] * Hq̇q[i] for i in eachindex(d)) .+ transpose(Jq̇)*Jq
end

function measure_hessians!(∂²g∂q²,∂²g∂q̇²,∂²g∂q̇∂q,structure::AbstractStructure,err_gau::ErrorGauge{sigType,captumType,<:Function}) where {sigType,captumType}
    (;signifier,captum,reference) = err_gau
    t = structure.state.system.t
    ref = reference(t)
    len = length(ref)
    d = zeros(eltype(ref), len)
    measure!(d, structure, signifier, captum)
    d .-= ref
    nq = get_num_of_full_coords(structure)
    ns = get_num_of_aux_var(structure)
    Jq = zeros(eltype(d), len, nq)
    Jq̇ = similar(Jq)
    Js = zeros(eltype(d), len, ns)
    measure_jacobian!(Jq, Jq̇, Js, structure, signifier, captum)
    Hqq, Hq̇q̇, Hq̇q = measure_hessians(structure,signifier,captum)

    ∂²g∂q²[:,:]  .+= sum(d[i] * Hqq[i] for i in eachindex(d)) .+ transpose(Jq)*Jq
    ∂²g∂q̇²[:,:]  .+= sum(d[i] * Hq̇q̇[i] for i in eachindex(d)) .+ transpose(Jq̇)*Jq̇
    ∂²g∂q̇∂q[:,:] .+= sum(d[i] * Hq̇q[i] for i in eachindex(d)) .+ transpose(Jq̇)*Jq
end

function measure_hessians!(∂²g∂q²,∂²g∂q̇²,∂²g∂q̇∂q,structure::AbstractStructure,gauge::CaptumGauge)
    (;signifier,captum,) = gauge
    len = get_num_of_capta(gauge)
    d = zeros(eltype(structure.state.system.q), len)
    measure!(d, structure, signifier, captum)
    nq = get_num_of_full_coords(structure)
    ns = get_num_of_aux_var(structure)
    Jq = zeros(eltype(d), len, nq)
    Jq̇ = similar(Jq)
    Js = zeros(eltype(d), len, ns)
    measure_jacobian!(Jq, Jq̇, Js, structure, signifier, captum)
    Hqq, Hq̇q̇, Hq̇q = measure_hessians(structure,signifier,captum)
    ∂²g∂q²[:,:]  .+= sum(d[i] * Hqq[i] for i in eachindex(d)) .+ transpose(Jq)*Jq
    ∂²g∂q̇²[:,:]  .+= sum(d[i] * Hq̇q̇[i] for i in eachindex(d)) .+ transpose(Jq̇)*Jq̇
    ∂²g∂q̇∂q[:,:] .+= sum(d[i] * Hq̇q[i] for i in eachindex(d)) .+ transpose(Jq̇)*Jq
end

# --- AngularPositionCaptum Implementation ---

function measure(structure::AbstractStructure, signifier::Hen2Egg, captum::AngularPositionCaptum)
    # 1. Get Relative Position [dx, dy]
    T = get_numbertype(structure)
    d_pos = zeros(T, 2) 
    measure!(d_pos, structure, signifier, PositionCaptum()) 
    
    atan(d_pos[2], d_pos[1])
end

function measure!(out, structure::AbstractStructure, signifier::Hen2Egg, captum::AngularPositionCaptum)
    out[1] = measure(structure, signifier, captum)
end

function measure_jacobian!(∂g∂q, ∂g∂q̇, ∂g∂s, structure::AbstractStructure, signifier::Hen2Egg, ::AngularPositionCaptum)
   
    ∂g∂q .= 0
    ∂g∂q̇ .= 0
    ∂g∂s .= 0
     
    nq = get_num_of_full_coords(structure)
    d_pos = zeros(eltype(∂g∂q), 2)
    measure!(d_pos, structure, signifier, PositionCaptum())
    
    J_pos_q = zeros(eltype(∂g∂q), 2, nq)
    J_pos_q̇ = zeros(eltype(∂g∂q), 2, nq)
    ns = size(∂g∂s, 2)
    J_pos_s = zeros(eltype(∂g∂s), 2, ns)
    measure_jacobian!(J_pos_q, J_pos_q̇, J_pos_s, structure, signifier, PositionCaptum())
    
    x = d_pos[1]
    y = d_pos[2]
    r2 = x^2 + y^2

    dθ_dx = -y / r2
    dθ_dy = x / r2
    
    for j in 1:nq
        ∂g∂q[1, j] = dθ_dx * J_pos_q[1, j] + dθ_dy * J_pos_q[2, j]
        ∂g∂q̇[1, j] = dθ_dx * J_pos_q̇[1, j] + dθ_dy * J_pos_q̇[2, j]
    end
end

function measure_hessians!(∂²g∂q², ∂²g∂q̇², ∂²g∂q̇∂q, structure::AbstractStructure, signifier::Hen2Egg, ::AngularPositionCaptum)
    nq = get_num_of_full_coords(structure)
    
    # 1. Get Position and Jacobians
    d_pos = zeros(eltype(∂²g∂q²), 2)
    measure!(d_pos, structure, signifier, PositionCaptum())
    
    J_pos_q = zeros(eltype(∂²g∂q²), 2, nq)
    J_pos_q̇ = zeros(eltype(∂²g∂q²), 2, nq)
    measure_jacobian!(J_pos_q, J_pos_q̇, structure, signifier, PositionCaptum())
    
    # Position Hessian? (Usually zero for Hen2Egg if linear coords? NC is linear)
    # But let's assume it exists to be safe, though Hen2Egg is linear diff.
    # measure_hessians! for PositionCaptum might be zero. 
    # Let's ignore term 2 (d^2 pos / dq^2) * dtheta/dpos
    # We focus on term 1: (d^2 theta / dpos^2) * (dpos/dq)^T * (dpos/dq)
    
    x = d_pos[1]
    y = d_pos[2]
    r2 = x^2 + y^2
    r4 = r2^2
    
    # H_theta_pos = [ d/dx(y/r2)  d/dy(y/r2) ]
    #               [ d/dx(-x/r2) d/dy(-x/r2) ]
    
    # d/dx(y/r2) = y * (-1/r4) * 2x = -2xy/r4
    # d/dy(-x/r2) = -x * (-1/r4) * 2y = 2xy/r4
    # d/dy(y/r2) = (1*r2 - y*2y)/r4 = (x^2 + y^2 - 2y^2)/r4 = (x^2 - y^2)/r4
    # d/dx(-x/r2) = (-1*r2 - (-x)*2x)/r4 = (-x^2 - y^2 + 2x^2)/r4 = (x^2 - y^2)/r4
    
    # H_pos = [ -2xy/r4    (x^2-y^2)/r4 ]
    #         [ (x^2-y^2)/r4   2xy/r4   ]
    
    H11 = -2*x*y / r4
    H12 = (x^2 - y^2) / r4
    H22 = 2*x*y / r4
    
    # J_pos_q is (2 x nq) matrix [Jx; Jy]
    # Chain rule: H_q = J_pos^T * H_theta * J_pos
    # H_q[i,j] = sum( J_pos[k,i] * H_theta[k,l] * J_pos[l,j] )
    
    for j in 1:nq, i in 1:nq
        Jx_i = J_pos_q[1, i]; Jy_i = J_pos_q[2, i]
        Jx_j = J_pos_q[1, j]; Jy_j = J_pos_q[2, j]
        
        val = Jx_i * (H11 * Jx_j + H12 * Jy_j) +
              Jy_i * (H12 * Jx_j + H22 * Jy_j)
              
        ∂²g∂q²[i, j] += val
    end
    
end

# --- AngularVelocityCaptum Implementation ---
function measure(structure::AbstractStructure, signifier::Hen2Egg, ::AngularVelocityCaptum)
    
    T = get_numbertype(structure)
    
    d_pos = zeros(T, 2)
    measure!(d_pos, structure, signifier, PositionCaptum())
    
    d_vel = zeros(T, 2)
    measure!(d_vel, structure, signifier, VelocityCaptum())
    
    return cartesian_to_angular_velocity(d_pos[1], d_pos[2], d_vel[1], d_vel[2])
end

function measure!(out, structure::AbstractStructure, signifier::Hen2Egg, captum::AngularVelocityCaptum)
    out[1] = measure(structure, signifier, captum)
end

function measure_jacobian!(∂g∂q, ∂g∂q̇, ∂g∂s, structure::AbstractStructure, signifier::Hen2Egg, ::AngularVelocityCaptum)
    
    ∂g∂q .= 0
    ∂g∂q̇ .= 0
    ∂g∂s .= 0

    nq = get_num_of_full_coords(structure)
    ns = size(∂g∂s,2)
    
    d_pos = zeros(eltype(∂g∂q), 2)
    measure!(d_pos, structure, signifier, PositionCaptum())
    
    d_vel = zeros(eltype(∂g∂q), 2)
    measure!(d_vel, structure, signifier, VelocityCaptum())
    
    J_pos_q = zeros(eltype(∂g∂q), 2, nq)
    J_pos_q̇ = zeros(eltype(∂g∂q), 2, nq)
    J_pos_s = zeros(eltype(∂g∂s), 2, ns)
    measure_jacobian!(J_pos_q, J_pos_q̇, J_pos_s, structure, signifier, PositionCaptum())
    
    J_vel_q = zeros(eltype(∂g∂q), 2, nq)
    J_vel_q̇ = zeros(eltype(∂g∂q), 2, nq)
    J_vel_s = zeros(eltype(∂g∂s), 2, ns)
    measure_jacobian!(J_vel_q, J_vel_q̇, J_vel_s, structure, signifier, VelocityCaptum())
    
    # Use ForwardDiff to compute gradients
    # This avoids manual derivation errors and stays consistent with the measure function
    
    # Define a helper function that takes the full state vector [q; q̇; s]
    # and returns the measurement.
    # However, measure! works on structure. We need to temporarily Mutate structure
    # or use a functional approach. 
    # Since structure state is mutable, we can use a closure but need to be careful.
    # Better approach: Extract q, q̇ from input vector and call measure on a shadow structure?
    # Or, given we have J_pos and J_vel from sub-captums, we can just use Chain Rule on 
    # the function f(pos, vel) = angular_velocity(pos, vel).
    
    # f(p, v) = (p[1]*v[2] - p[2]*v[1]) / (p[1]^2 + p[2]^2)
    # We need df/dp and df/dv.
    
    f_ang_vel = (x) -> begin
        p = x[1:2]
        v = x[3:4]
        cartesian_to_angular_velocity(p[1], p[2], v[1], v[2])
    end
    
    # Current values
    x_val = [d_pos; d_vel] # [x, y, vx, vy]
    
    gradient = ForwardDiff.gradient(f_ang_vel, x_val)
    
    dω_dp = gradient[1:2]
    dω_dv = gradient[3:4]
    
    # Chain Rule:
    # dω/dq = dω/dp * dp/dq + dω/dv * dv/dq
    
    # Manual loop to avoid allocations and shape mismatch issues (1xN vs N)
    for j in 1:nq
        ∂g∂q[1, j] = dω_dp[1] * J_pos_q[1, j] + dω_dp[2] * J_pos_q[2, j] +
                     dω_dv[1] * J_vel_q[1, j] + dω_dv[2] * J_vel_q[2, j]
                     
        ∂g∂q̇[1, j] = dω_dp[1] * J_pos_q̇[1, j] + dω_dp[2] * J_pos_q̇[2, j] +
                     dω_dv[1] * J_vel_q̇[1, j] + dω_dv[2] * J_vel_q̇[2, j]
    end

    for j in 1:ns
        ∂g∂s[1, j] = dω_dp[1] * J_pos_s[1, j] + dω_dp[2] * J_pos_s[2, j] +
                     dω_dv[1] * J_vel_s[1, j] + dω_dv[2] * J_vel_s[2, j]
    end
end

# Angular Apparatus

function measure(structure::AbstractStructure, appar::Apparatus, captum::AngularPositionCaptum)
    # 1. Get Relative Position [dx, dy]
    (;apparid2sys_aux_var_idx) = structure.connectivity
    (;s) = structure.state.system
    id = appar.id
    # 2. Return Angle 
    return s[apparid2sys_aux_var_idx[id]]
end

function measure!(out, structure::AbstractStructure, appar::Apparatus, captum::AngularPositionCaptum)
    
    out[1] = measure(structure, appar, captum)[begin]
end

function measure_jacobian!(∂g∂q, ∂g∂q̇, ∂g∂s, structure::AbstractStructure, appar::Apparatus, ::AngularPositionCaptum)
    
    ∂g∂q .= 0
    ∂g∂q̇ .= 0
    ∂g∂s .= 0
    (;apparid2sys_aux_var_idx) = structure.connectivity
    id = appar.id
    idx = apparid2sys_aux_var_idx[id]
    for j in idx
        ∂g∂s[1, j] = 1
    end
end
