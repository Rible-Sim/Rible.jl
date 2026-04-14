"""
Shared Jacobian computations for Zhong06 CCP Constant Mass Monolithic solver.

This module provides unified Jacobian computation routines used by:
- Primal solver (Newton iterations)
- Adjoint sensitivity solver (discrete adjoint)
- Direct sensitivity solver (forward sensitivity)

Ensures consistency across all three methods and reduces code duplication.
"""

"""
Create line search functions (ϕ, dϕ, ϕdϕ) with shared setup for efficient evaluation.

Returns a tuple of three closures:
- ϕ(β): merit function ϕ(β) = 0.5 * ||Res(x + β*Δx)||²
- dϕ(β): derivative dϕ(β)/dβ
- ϕdϕ(β): both merit and derivative (for efficiency)

All three functions share the same temporary workspace setup to avoid redundant allocations.
"""
function create_line_search_functions(
        workspace::NewtonWorkspace,
        workspace_temp::NewtonWorkspace,
        solver_state, solver_cache, bot, policy, env
    )
    
    # Merit function ϕ(β)
    ϕ = function(β)
        workspace.xₖ .= workspace.x .+ β.*workspace.Δx
        compute_constant_mass_residual!(workspace_temp.Res, workspace_temp.x, solver_state, solver_cache, bot, policy, env)
        return 0.5 * transpose(workspace.Res)*workspace.Res
    end
    
    # Derivative dϕ(β)/dβ
    dϕ = function(β)
        workspace.xₖ .= workspace.x .+ β.*workspace.Δx
        compute_constant_mass_residual!(workspace_temp.Res, workspace_temp.x, solver_state, solver_cache, bot, policy, env)
        compute_constant_mass_jacobian!(workspace_temp.Jac, workspace_temp.x, solver_state, solver_cache, bot, policy, env)
        return dot(workspace.Res, workspace.Jac * workspace.Δx)
    end
    
    # Combined function ϕdϕ(β) - more efficient than calling both separately
    ϕdϕ = function(β)
        workspace.xₖ .= workspace.x .+ β.*workspace.Δx
        compute_constant_mass_residual!(workspace_temp.Res, workspace_temp.x, solver_state, solver_cache, bot, policy, env)
        compute_constant_mass_jacobian!(workspace_temp.Jac, workspace_temp.x, solver_state, solver_cache, bot, policy, env)
        f = 0.5 * transpose(workspace.Res)*workspace.Res
        df = dot(workspace.Res, workspace.Jac * workspace.Δx)
        return f, df
    end
    
    return ϕ, dϕ, ϕdϕ
end

function create_line_search_functions(
        workspace::NewtonWorkspace,
        workspace_temp::NewtonWorkspace,
        call_residual!::Function,
        call_jacobian!::Function
    )
    
    # Merit function ϕ(β)
    ϕ = function(β)
        workspace.xₖ .= workspace.x .+ β.*workspace.Δx
        call_residual!(workspace_temp)
        return 0.5 * transpose(workspace.Res)*workspace.Res
    end
    
    # Derivative dϕ(β)/dβ
    dϕ = function(β)
        workspace.xₖ .= workspace.x .+ β.*workspace.Δx
        call_residual!(workspace_temp)
        call_jacobian!(workspace_temp)
        return dot(workspace.Res, workspace.Jac * workspace.Δx)
    end
    
    # Combined function ϕdϕ(β) - more efficient than calling both separately
    ϕdϕ = function(β)
        workspace.xₖ .= workspace.x .+ β.*workspace.Δx
        call_residual!(workspace_temp)
        call_jacobian!(workspace_temp)
        f = 0.5 * transpose(workspace.Res)*workspace.Res
        df = dot(workspace.Res, workspace.Jac * workspace.Δx)
        return f, df
    end
    
    return ϕ, dϕ, ϕdϕ
end

function interpolate!(solver_state::Zhong06SolverState)
    h = solver_state.dt
    (;state_k,state_kp1,state_mid) = solver_state
    @. state_mid.q = (state_k.q + state_kp1.q)/2
    @. state_mid.q̇ = (state_kp1.q - state_k.q)/h
    @. state_mid.p = (state_k.p + state_kp1.p)/2 #mind
    state_mid.t = (state_k.t + state_kp1.t)/2
    @. state_mid.s = (state_k.s + state_kp1.s)/2
    @. state_mid.λ = (state_k.λ + state_kp1.λ)/2 #caution
    @. state_mid.F = (state_k.F + state_kp1.F)/2 #caution
end
# ============================================================================
# Core Jacobian Computation
# ============================================================================


function compute_zhong06_jacobian_blocks!(
        jac_blocks::Zhong06JacobianBlocks,
        jacobian_workspace::Zhong06JacobianWorkspace,
        solver_state::Zhong06SolverState,
        consts::Zhong06Constants,
        contact_cache,
        bot, policy, field, forward_cache
    )
    (;structure, hub) = bot
    # Unpack jao acobian_workspace
    (;Fₘ, ∂F∂q, ∂F∂q̇, ∂C∂qₖ, ∂C∂pₖ, ∂C∂qₖ₊₁, ∂C∂pₖ₊₁, ∂Fₘ∂u, ∂Fₘ∂c, 
      Mₘ, M⁻¹ₘ, ∂Mₘhq̇ₘ∂qₘ, Aₖ₊₁, Aₖ, ∂Aᵀλ∂q) = jacobian_workspace
    
    # Unpack constants
    (;h, mass_norm, nq, nλ, nu, nc) = consts
    
    # Unpack solver state
    (;qₖ₊₁, qₖ, pₖ₊₁, pₖ, λₘ, qₘ, q̇ₘ, vₖ₊₁, vₖ, tₘ, tₖ, tₖ₊₁, Λₖ₊₁, Γₖ₊₁) = solver_state
    (;state_k, state_kp1, state_mid) = solver_state
    # Unpack Jacobian blocks
    (;Jacᵏ⁺¹ₖ₊₁, Jacᵏ⁺¹ₖ, Jacᵏ⁺¹ₘu, Jacᵏ⁺¹ₘc) = jac_blocks
    
    T = eltype(qₖ)
    n1 = nq
    n2 = 2nq
    n3 = n2 + nλ
    
    # Make structure-specific functions
    
    # Compute mass matrix and its Jacobian

    # commented out for constant mass
    # update_bodies!(structure,state_mid)
    # assemble_M!(Mₘ,structure)
    # assemble_∂Mq̇∂q!(∂Mₘhq̇ₘ∂qₘ,structure)
    # ∂Mₘhq̇ₘ∂qₘ .*=h

    # Compute forces and their Jacobians # Compute control Jacobians
    control_jacobian!(jacobian_workspace, bot, policy, forward_cache, solver_state)
    
    gen_force!(state_mid, bot, field, policy)
    Fₘ .= state_mid.F
    gen_force_jacobian!(∂F∂q, ∂F∂q̇, ∂Fₘ∂u, ∂Fₘ∂c, bot, field, policy, state_mid)
    
   
    # Compute constraint Jacobians
    cstr_jacobian!(Aₖ₊₁,structure,state_kp1)
    cstr_jacobian!(Aₖ, structure,state_k)
    
    # ========================================================================
    # Core (non-contact) Jacobian blocks
    # ========================================================================
    
    # Jacᵏ⁺¹ₖ: Derivative w.r.t. previous state (qₖ, pₖ, λₖ₋₁)
    cstr_forces_jacobian!(∂Aᵀλ∂q, structure, qₖ, λₘ)
    Jacᵏ⁺¹ₖ .= 0.0
    Jacᵏ⁺¹ₖ[   1:n1,    1:n1] .=  1/2 .*∂Mₘhq̇ₘ∂qₘ .- Mₘ .- (h^2)/2 .*(1/2 .*∂F∂q .- 1/h.*∂F∂q̇ .+ ∂C∂qₖ) .- 
                                   mass_norm.*∂Aᵀλ∂q
    Jacᵏ⁺¹ₖ[n1+1:n2,    1:n1] .= -1/2 .*∂Mₘhq̇ₘ∂qₘ .+ Mₘ .- (h^2)/2 .*(1/2 .*∂F∂q .- 1/h.*∂F∂q̇ .+ ∂C∂qₖ)
    Jacᵏ⁺¹ₖ[   1:n1, n1+1:n2] .= -h*I(n1) .-(h^2)/2 .*∂C∂pₖ
    Jacᵏ⁺¹ₖ[n1+1:n2, n1+1:n2] .= .-(h^2)/2 .*∂C∂pₖ
    
    # Jacᵏ⁺¹ₖ₊₁: Derivative w.r.t. current state (qₖ₊₁, pₖ₊₁, λₘ)
    cstr_forces_jacobian!(∂Aᵀλ∂q, structure, qₖ₊₁, λₘ)
    Jacᵏ⁺¹ₖ₊₁ .= 0.0
    Jacᵏ⁺¹ₖ₊₁[   1:n1,    1:n1] .=  1/2 .*∂Mₘhq̇ₘ∂qₘ .+ Mₘ .- (h^2)/2 .*(1/2 .*∂F∂q .+ 1/h.*∂F∂q̇ .+ ∂C∂qₖ₊₁)
    Jacᵏ⁺¹ₖ₊₁[n1+1:n2,    1:n1] .= -1/2 .*∂Mₘhq̇ₘ∂qₘ .- Mₘ .- (h^2)/2 .*(1/2 .*∂F∂q .+ 1/h.*∂F∂q̇ .+ ∂C∂qₖ₊₁) .- 
                                   mass_norm.*∂Aᵀλ∂q
    Jacᵏ⁺¹ₖ₊₁[n2+1:n3,    1:n1] .= -mass_norm.*Aₖ₊₁
    Jacᵏ⁺¹ₖ₊₁[   1:n1, n1+1:n2] .= -(h^2)/2 .*∂C∂pₖ₊₁
    Jacᵏ⁺¹ₖ₊₁[n1+1:n2, n1+1:n2] .=  h*I(n1) .- (h^2)/2 .*∂C∂pₖ₊₁
    Jacᵏ⁺¹ₖ₊₁[   1:n1, n2+1:n3] .= -mass_norm.*transpose(Aₖ)
    Jacᵏ⁺¹ₖ₊₁[n1+1:n2, n2+1:n3] .= -mass_norm.*transpose(Aₖ₊₁)
    
    # Parameter Jacobians
    Jacᵏ⁺¹ₘu .= 0.0
    Jacᵏ⁺¹ₘc .= 0.0
    Jacᵏ⁺¹ₘu[   1:n1, 1:nu] .= -(h^2)/2 .*∂Fₘ∂u
    Jacᵏ⁺¹ₘu[n1+1:n2, 1:nu] .= -(h^2)/2 .*∂Fₘ∂u
    Jacᵏ⁺¹ₘu[n2+1:n3, 1:nu] .= 0.0
    Jacᵏ⁺¹ₘc[   1:n1, 1:nc] .= -(h^2)/2 .*∂Fₘ∂c
    Jacᵏ⁺¹ₘc[n1+1:n2, 1:nc] .= -(h^2)/2 .*∂Fₘ∂c
    Jacᵏ⁺¹ₘc[n2+1:n3, 1:nc] .= 0.0
    
    # ========================================================================
    # Contact Jacobian blocks (if contacts are active)
    # ========================================================================
    
    (;na) = contact_cache
    if na != 0
        _compute_contact_jacobian_blocks!(
            jac_blocks, jacobian_workspace, solver_state, consts, contact_cache,
            structure, n1, n2, n3
        )
    end
    
    return nothing
end

"""
Compute contact-related Jacobian blocks.

This function adds the contribution from frictional contacts to the Jacobian.
"""
function _compute_contact_jacobian_blocks!(
        jac_blocks, workspace, solver_state, consts, contact_cache,
        structure, n1, n2, n3
    )
    
    (;Jacᵏ⁺¹ₖ₊₁, Jacᵏ⁺¹ₖ) = jac_blocks
    (;M⁻¹ₘ,) = workspace
    (;h, mass_norm) = consts
    (;qₖ, pₖ, pₖ₊₁, λₘ, Λₖ₊₁, Γₖ₊₁) = solver_state
    
    (;
        na,
        Λ, Γ,
        H,
        activated_restitution_coefficients,
        D,
        L, Lv,
        Dper, Dimp
    ) = contact_cache
    
    T = eltype(qₖ)
    nΛ = 3na
    n4 = n3 + nΛ
    n5 = n4 + nΛ
    
    Λ_split = split_by_lengths(Λ, 3)
    Γ_split = split_by_lengths(Γ, 3)
    
    Dₘ   = Dper
    Dₖ₊₁ = Dimp
    
    vₖ   = M⁻¹ₘ*pₖ
    vₖ₊₁ = M⁻¹ₘ*pₖ₊₁
    
    # Compute velocity Jacobians
    ∂vₘ∂qₖ   = -1/h * I
    ∂vₘ∂qₖ₊₁ =  1/h * I
    ∂vₖ∂pₖ   = M⁻¹ₘ
    ∂vₖ₊₁∂pₖ₊₁ = M⁻¹ₘ
    
    # Compute contact velocity
    vₘ = (solver_state.qₖ₊₁ .- solver_state.qₖ)./h
    v́⁺ = Dₘ*vₘ .+ Dₖ₊₁*vₖ₊₁
    
    # Compute Γ Jacobians
    ∂v́⁺∂qₖ   = Dₘ*∂vₘ∂qₖ
    ∂v́⁺∂qₖ₊₁ = Dₘ*∂vₘ∂qₖ₊₁
    ∂v́⁺∂pₖ₊₁ = Dₖ₊₁*∂vₖ₊₁∂pₖ₊₁
    
    ∂Γ∂qₖ = zeros(T, nΛ, n1)
    ∂Γ∂pₖ = zeros(T, nΛ, n1)
    ∂Γ∂qₖ₊₁ = zeros(T, nΛ, n1)
    ∂Γ∂pₖ₊₁ = zeros(T, nΛ, n1)
    ∂Γ∂qₖ   .= ∂v́⁺∂qₖ
    ∂Γ∂pₖ   .= 0.0
    ∂Γ∂qₖ₊₁ .= ∂v́⁺∂qₖ₊₁
    ∂Γ∂pₖ₊₁ .= ∂v́⁺∂pₖ₊₁
    
    # Add restitution and tangential velocity contributions
    v́ₖ = Dₖ₊₁*vₖ
    ∂v́ₖ∂pₖ = Dₖ₊₁*∂vₖ∂pₖ
    
    for i = 1:na
        is = 3(i-1)
        vⁱₖ = @view v́ₖ[is+1:is+3]
        vₙⁱₖ = vⁱₖ[1]
        
        # Restitution contribution
        ∂Γ∂pₖ[is+1, 1:n1] .+= activated_restitution_coefficients[i] * (vₙⁱₖ < 0) * ∂v́ₖ∂pₖ[is+1, :]
        
        # Tangential velocity norm contribution
        v́ₜ_norm = norm(v́⁺[is+2:is+3]) + 1e-14
        ∂Γ∂qₖ[  is+1, 1:n1] .+= (v́⁺[is+2]*∂v́⁺∂qₖ[  is+2, :] .+ v́⁺[is+3]*∂v́⁺∂qₖ[  is+3, :]) / v́ₜ_norm
        ∂Γ∂qₖ₊₁[is+1, 1:n1] .+= (v́⁺[is+2]*∂v́⁺∂qₖ₊₁[is+2, :] .+ v́⁺[is+3]*∂v́⁺∂qₖ₊₁[is+3, :]) / v́ₜ_norm
        ∂Γ∂pₖ₊₁[is+1, 1:n1] .+= (v́⁺[is+2]*∂v́⁺∂pₖ₊₁[is+2, :] .+ v́⁺[is+3]*∂v́⁺∂pₖ₊₁[is+3, :]) / v́ₜ_norm
    end
    
    # Add contact blocks to Jacobians
    Jacᵏ⁺¹ₖ[  n3+1:n4,    1:n1] .=  h.*∂Γ∂qₖ
    Jacᵏ⁺¹ₖ[  n3+1:n4, n1+1:n2] .=  h.*∂Γ∂pₖ
    
    Jacᵏ⁺¹ₖ₊₁[n3+1:n4,    1:n1] .=  h.*∂Γ∂qₖ₊₁
    Jacᵏ⁺¹ₖ₊₁[n3+1:n4, n1+1:n2] .=  h.*∂Γ∂pₖ₊₁
    Jacᵏ⁺¹ₖ₊₁[   1:n1, n3+1:n4] .= -mass_norm*h.*transpose(D)*H*(I+L)
    Jacᵏ⁺¹ₖ₊₁[n1+1:n2, n3+1:n4] .= -mass_norm*h.*transpose(D)*H*(I+L)
    Jacᵏ⁺¹ₖ₊₁[n4+1:n5, n3+1:n4] .=  BlockDiagonal(mat.(Γ_split))
    Jacᵏ⁺¹ₖ₊₁[n3+1:n4, n4+1:n5] .= -h.*I(nΛ)
    Jacᵏ⁺¹ₖ₊₁[n4+1:n5, n4+1:n5] .=  BlockDiagonal(mat.(Λ_split))
    
    return nothing
end

# ============================================================================
# CCP Mono Primal Solver Residual and Jacobian Functions
# ============================================================================

"""
Compute residual for primal Newton solver at timestep k+1.

Note: This uses a reduced state vector [qₖ₊₁, λₘ, Λ, Γ] where pₖ₊₁ is computed
from qₖ₊₁ via the momentum relation, unlike the adjoint/direct solvers which
treat p as an independent variable.
"""
function compute_primal_residual!(workspace,  # NewtonWorkspace containing Res, x, 𝐰, pₖ₊₁ (buffer), q̇ₖ₊₁ (buffer)
        solver_state::Zhong06SolverState,  # Unified solver state
        solver_cache, contact_cache
    )
    # Unpack for convenience
    𝐫𝐞𝐬 = workspace.Res
    x = workspace.x
    (;qₖ, q̇ₖ, pₖ, qₖ₊₁,pₖ₊₁, vₖ₊₁, qₘ, q̇ₘ, λₘ, tₖ, tₖ₊₁, dt) = solver_state
    (;state_k, state_kp1, state_mid) = solver_state
    vₖ = q̇ₖ
    h = dt
    𝐰 = workspace.𝐰
    
    (;bot, policy, field, consts) = solver_cache
    (;Mₘ, M⁻¹ₘ, Fₘ, Aₖ, Aₖ₊₁) = solver_cache.jacobian_workspace
    (;mass_norm) = consts
    (;hub, structure) = bot
    
    cstr_jacobian!(Aₖ, structure, state_k)
    cstr_jacobian!(Aₖ₊₁, structure, state_kp1)
    nq = length(qₖ)
    nλ = size(Aₖ, 1)
    na = contact_cache.na
    nΛ = 3na
    n1 = nq
    n2 = nq + nλ
    
    # Extract state variables
    qₖ₊₁ .= x[1:n1]
    λₘ .= x[n1+1:n2]
    
    # Compute midpoint quantities
    tₘ = (tₖ + tₖ₊₁) / 2
    qₘ .= (qₖ₊₁ .+ qₖ) ./ 2
    q̇ₘ .= (qₖ₊₁ .- qₖ) ./ h
    vₘ = q̇ₘ
    
    # Compute momentum (key: pₖ₊₁ is derived, not independent)
    pₖ₊₁ .= -pₖ .+ 2/h.*Mₘ*(qₖ₊₁.-qₖ) .+ mass_norm/h.*(transpose(Aₖ₊₁)-transpose(Aₖ))*λₘ
    control!(bot, policy, solver_cache, solver_state)
    
    # Compute generalized forces
    gen_force!(state_mid, bot, field, policy)
    Fₘ .= state_mid.F
    # Core residual (momentum and constraints)

    𝐫𝐞𝐬 .= 0.0
    𝐫𝐞𝐬[1:n1] .= -h .* pₖ .+ Mₘ * (qₖ₊₁ .- qₖ) .-
                  mass_norm .* transpose(Aₖ) * λₘ .-
                  (h^2) / 2 .* Fₘ
    ϕ = @view 𝐫𝐞𝐬[n1+1:n2]
    cstr_function!(ϕ, structure,state_kp1)
    ϕ .*= -mass_norm
    # Contact residual
    if na != 0
        Λ = @view x[(n2+1):n2+nΛ]
        Γ = @view x[n2+nΛ+1:n2+2nΛ]
        Λ_split = split_by_lengths(Λ, 3)
        Γ_split = split_by_lengths(Γ, 3)
        
        get_distribution_law!(structure, contact_cache, state_kp1)
        (;H, activated_restitution_coefficients, D, L) = contact_cache
        Dₘ = contact_cache.Dper
        Dₖ₊₁ = contact_cache.Dimp
        
        vₖ₊₁ .= M⁻¹ₘ * pₖ₊₁
        v́⁺ = Dₘ * vₘ .+ Dₖ₊₁ * vₖ₊₁
        ν = zero(v́⁺)
        v́ₖ = Dₖ₊₁ * vₖ
        
        for i = 1:na
            is = 3(i - 1)
            vⁱₖ = @view v́ₖ[is+1:is+3]
            vⁱ⁺ = @view v́⁺[is+1:is+3]
            vₙⁱₖ = vⁱₖ[1]
            vₜⁱ⁺ = norm(vⁱ⁺[2:3])
            
            𝐰[is+1:is+3] .= [activated_restitution_coefficients[i] * min(vₙⁱₖ, 0), 0, 0]
            ν[is+1:is+3] .= vⁱ⁺ + [vₜⁱ⁺, 0, 0]
        end
        
        𝐫𝐞𝐬[1:n1] .-= h .* mass_norm .* transpose(D) * H * (I + L) * Λ
        𝐫𝐞𝐬[(n2+1):(n2+nΛ)] .= h .* (ν .+ 𝐰 .- Γ)
        𝐫𝐞𝐬[(n2+nΛ+1):(n2+2nΛ)] .= reduce(vcat, Λ_split ⊙ Γ_split)
    end
    
    return nothing
end

"""
Compute Jacobian for primal Newton solver at timestep k+1.

Key difference from adjoint/direct: pₖ₊₁ = pₖ₊₁(qₖ₊₁, λₘ) is not independent,
so we use chain rule: ∂vₖ₊₁/∂qₖ₊₁ = (∂vₖ₊₁/∂pₖ₊₁) * (∂pₖ₊₁/∂qₖ₊₁).
"""
function compute_primal_jacobian!(
        workspace,  # NewtonWorkspace containing Jac, x, ∂Γ∂x, pₖ₊₁ (buffer), q̇ₖ₊₁ (buffer)
        solver_state::Zhong06SolverState,  # Unified solver state
        solver_cache,
        contact_cache
    )
    # Unpack for convenience
    𝐉 = workspace.Jac
    x = workspace.x
    (;qₖ, q̇ₖ, pₖ, qₖ₊₁, pₖ₊₁, vₖ₊₁, qₘ, q̇ₘ, λₘ, tₖ, tₖ₊₁, dt) = solver_state
    (; state_k, state_kp1, state_mid) = solver_state
    vₖ = q̇ₖ
    h = dt
    ∂Γ∂x = workspace.∂Γ∂x
    
    (;bot, policy, field, consts) = solver_cache
    (;Mₘ, M⁻¹ₘ, ∂F∂q, ∂F∂q̇, ∂Fₘ∂u, ∂F∂s, ∂C∂qₖ₊₁, ∂C∂pₖ₊₁, Fₘ, Aₖ, Aₖ₊₁,∂Aᵀλ∂q) = solver_cache.jacobian_workspace
    (;mass_norm) = consts
    (;hub, structure) = bot

    nq = length(qₖ)
    nλ = size(Aₖ, 1)
    n1 = nq
    n2 = nq + nλ
    na = contact_cache.na
    nΛ = 3na
    
    # Extract state variables
    qₖ₊₁ .= x[1:n1]
    λₘ .= x[n1+1:n2]
    
    # Compute midpoint quantities
    qₘ .= (qₖ₊₁ .+ qₖ) ./ 2
    q̇ₘ .= (qₖ₊₁ .- qₖ) ./ h
    vₘ = q̇ₘ
    tₘ = (tₖ + tₖ₊₁) / 2
    
    # Compute momentum (pₖ₊₁ is a function of qₖ₊₁, λₘ)
    pₖ₊₁ .= -pₖ .+ 2/h.*Mₘ*(qₖ₊₁.-qₖ) .+ mass_norm/h.*(transpose(Aₖ₊₁)-transpose(Aₖ))*λₘ
    
    control_jacobian!(solver_cache.jacobian_workspace, bot, policy, solver_cache, solver_state)
    cstr_forces_jacobian!(∂Aᵀλ∂q, structure,qₖ₊₁,λₘ)
    ∂p∂q = 2/h.*Mₘ  .+ mass_norm/h.*∂Aᵀλ∂q
    ∂C∂qₖ₊₁ .+= ∂C∂pₖ₊₁*∂p∂q
    #note: also need ∂p∂λ probably

    # Compute force Jacobians
    
    gen_force!(solver_state.state_mid, bot, field, policy, )
    gen_force_state_jacobian!(∂F∂q, ∂F∂q̇, ∂Fₘ∂u, bot, field, policy, solver_state.state_mid, ∂F∂s)

    # Constraint Jacobians
    cstr_jacobian!(Aₖ, structure, state_k)
    cstr_jacobian!(Aₖ₊₁, structure, state_kp1)
    
    # ========================================================================
    # CORE JACOBIAN (No Contacts)
    # ========================================================================
    𝐉 .= 0.0
    
    # ∂R₁/∂qₖ₊₁: Momentum balance equation w.r.t. qₖ₊₁
    𝐉[1:n1, 1:n1] .= Mₘ .- h^2 / 2 .* (1/2 .* ∂F∂q .+ 1/h .* ∂F∂q̇ .+ ∂C∂qₖ₊₁)
    
    # ∂R₁/∂λₘ: Momentum balance equation w.r.t. λₘ
    𝐉[1:n1, n1+1:n2] .= -mass_norm .* transpose(Aₖ)
    
    # ∂R₂/∂qₖ₊₁: Constraint equation w.r.t. qₖ₊₁
    𝐉[n1+1:n2, 1:n1] .= -mass_norm .* Aₖ₊₁
    
    # ========================================================================
    # CONTACT JACOBIAN (if contacts active)
    # ========================================================================
    if na != 0
        Λ = @view x[(n2+1):n2+nΛ]
        Γ = @view x[n2+nΛ+1:n2+2nΛ]
        Λ_split = split_by_lengths(Λ, 3)
        Γ_split = split_by_lengths(Γ, 3)
        
        get_distribution_law!(structure, contact_cache, state_kp1)
        (;H, activated_restitution_coefficients, D, L) = contact_cache
        Dₘ = contact_cache.Dper
        Dₖ₊₁ = contact_cache.Dimp
        
        vₖ₊₁ .= M⁻¹ₘ * pₖ₊₁
        
        # Chain rule: ∂vₖ₊₁/∂qₖ₊₁ = (∂vₖ₊₁/∂pₖ₊₁) * (∂pₖ₊₁/∂qₖ₊₁)
        # where pₖ₊₁ = -pₖ + (2/h)*Mₘ*(qₖ₊₁-qₖ) + (mass_norm/h)*(Aₖ₊₁ᵀ-Aₖᵀ)*λₘ
        # So: ∂pₖ₊₁/∂qₖ₊₁ = (2/h)*Mₘ + (mass_norm/h)*∂Aₖ₊₁ᵀ/∂qₖ₊₁*λₘ
        # And: ∂vₖ₊₁/∂qₖ₊₁ = M⁻¹ₘ*∂pₖ₊₁/∂qₖ₊₁ = (2/h)*I + (mass_norm/h)*M⁻¹ₘ*∂Aₖ₊₁ᵀ/∂qₖ₊₁*λₘ
        ∂vₘ∂qₖ₊₁ = 1 / h * I
        cstr_forces_jacobian!(∂Aᵀλ∂q, structure, qₖ₊₁, λₘ)
        ∂vₖ₊₁∂qₖ₊₁ = 2 / h * I + mass_norm / h .* M⁻¹ₘ * ∂Aᵀλ∂q
        
        # Similarly: ∂vₖ₊₁/∂λₘ = M⁻¹ₘ * ∂pₖ₊₁/∂λₘ
        ∂vₖ₊₁∂λₘ = mass_norm .* M⁻¹ₘ * transpose(Aₖ₊₁ - Aₖ) / h
        
        v́⁺ = Dₘ * vₘ .+ Dₖ₊₁ * vₖ₊₁
        ∂v́⁺∂qₖ₊₁ = Dₘ * ∂vₘ∂qₖ₊₁ .+ Dₖ₊₁ * ∂vₖ₊₁∂qₖ₊₁
        ∂v́⁺∂λₘ = Dₖ₊₁ * ∂vₖ₊₁∂λₘ
        
        ∂Γ∂x .= 0
        ∂Γ∂x[:, 1:n1] .= ∂v́⁺∂qₖ₊₁
        ∂Γ∂x[:, n1+1:n2] .= ∂v́⁺∂λₘ
        
        # Add tangential velocity contribution
        for i = 1:na
            is = 3(i - 1)
            v́ₜ_norm = norm(v́⁺[is+2:is+3]) + 1e-14
            ∂Γ∂x[is+1, 1:n1] .+= (v́⁺[is+2] * ∂v́⁺∂qₖ₊₁[is+2, :] .+ v́⁺[is+3] * ∂v́⁺∂qₖ₊₁[is+3, :]) / v́ₜ_norm
            ∂Γ∂x[is+1, n1+1:n2] .+= (v́⁺[is+2] * ∂v́⁺∂λₘ[is+2, :] .+ v́⁺[is+3] * ∂v́⁺∂λₘ[is+3, :]) / v́ₜ_norm
        end
        
        # Assemble contact blocks into Jacobian
        𝐉[(n2+1):(n2+nΛ), 1:n2] .= h .* ∂Γ∂x
        𝐉[1:n1, (n2+1):(n2+nΛ)] .= -mass_norm * h .* transpose(D) * H * (I + L)
        𝐉[(n2+nΛ+1):(n2+2nΛ), (n2+1):(n2+nΛ)] .= BlockDiagonal(mat.(Γ_split))
        𝐉[(n2+1):(n2+nΛ), (n2+nΛ+1):(n2+2nΛ)] .= -h .* I(nΛ)
        𝐉[(n2+nΛ+1):(n2+2nΛ), (n2+nΛ+1):(n2+2nΛ)] .= BlockDiagonal(mat.(Λ_split))
    end
    
    return nothing
end

function Zhong06_constant_mass_Res!(
        Res,
        solver_state,
        bot,env,policy,forward_cache,
        consts,jacobian_workspace,
    )
    (; qₖ₊₁, qₖ, pₖ₊₁, pₖ, λₘ, dt) = solver_state
    h = dt
    (; state_k, state_kp1, state_mid) = solver_state
    (; Mₘ, Aₖ, Aₖ₊₁) = jacobian_workspace
    (; mass_norm, nq, nλ) = consts
    (;structure) = bot
    control!(bot, policy, forward_cache, solver_state)
    gen_force!(state_mid,bot,env.field,policy)
    Fₘ = state_mid.F
    cstr_jacobian!(Aₖ, structure, state_k)
    cstr_jacobian!(Aₖ₊₁, structure, state_kp1)
    hpₘ = Mₘ*(qₖ₊₁.-qₖ)
    Res[    1:    nq    ] .= (hpₘ     .- h.*pₖ) .- mass_norm.*transpose(Aₖ  )*λₘ .- (h^2)/2 .*Fₘ
    Res[ nq+1:   2nq    ] .= (h.*pₖ₊₁ .- hpₘ  ) .- mass_norm.*transpose(Aₖ₊₁)*λₘ .- (h^2)/2 .*Fₘ
    ϕ = @view Res[2nq+1:2nq+nλ]
    cstr_function!(ϕ, structure,state_kp1)
    ϕ .*= -mass_norm
    # Res[2nq+nλ+1:end] .= S(qₘ,sₘ)
end

function Zhong06_constant_mass_Jac!(
        jac_blocks::Zhong06JacobianBlocks,
        solver_state::Zhong06SolverState,
        bot,env,policy,forward_cache,
        consts::Zhong06Constants,jacobian_workspace::Zhong06JacobianWorkspace
    )
    (;structure) = bot
    # Unpack constants
    (;h,mass_norm,nq,nλ,nu,nc,ns) = consts
    # Unpack jacobian_workspace
    (;Mₘ,∂Mₘhq̇ₘ∂qₘ,∂F∂q,∂F∂q̇,∂F∂s,∂C∂qₖ, ∂C∂pₖ, ∂C∂sₖ, ∂C∂qₖ₊₁,∂C∂pₖ₊₁,∂Fₘ∂u,∂Fₘ∂c,Aₖ₊₁,Aₖ,∂Aᵀλ∂q,∂S∂q,∂S∂s) = jacobian_workspace
    # Unpack solver state
    (;qₖ₊₁,qₖ,pₖ₊₁,pₖ,λₘ,qₘ,q̇ₘ,tₘ,tₖ,tₖ₊₁) = solver_state
    # Unpack Jacobian blocks
    (;Jacᵏ⁺¹ₖ₊₁,Jacᵏ⁺¹ₖ,Jacᵏ⁺¹ₘu,Jacᵏ⁺¹ₘc) = jac_blocks

    control_jacobian!(jacobian_workspace, bot, policy, forward_cache, solver_state)
    
    gen_force_jacobian!(∂F∂q,∂F∂q̇,∂Fₘ∂u,∂Fₘ∂c,bot,env.field,policy,solver_state.state_mid,∂F∂s)
    # Compute constraint Jacobians

    cstr_jacobian!(Aₖ₊₁, structure, solver_state.state_kp1)
    cstr_forces_jacobian!(∂Aᵀλ∂q, structure, qₖ, λₘ)

    Jacᵏ⁺¹ₖ .= 0.0
    Jacᵏ⁺¹ₖ[   1:nq ,   1:nq ]     .=  1/2 .*∂Mₘhq̇ₘ∂qₘ .- Mₘ .-(h^2)/2 .*(1/2 .*∂F∂q .-1/h.*∂F∂q̇ .+ ∂C∂qₖ) .- mass_norm.*∂Aᵀλ∂q
    Jacᵏ⁺¹ₖ[nq+1:2nq,   1:nq ]     .= -1/2 .*∂Mₘhq̇ₘ∂qₘ .+ Mₘ .-(h^2)/2 .*(1/2 .*∂F∂q .-1/h.*∂F∂q̇ .+ ∂C∂qₖ)
    Jacᵏ⁺¹ₖ[   1:nq ,nq+1:2nq]     .= -h.*I(nq) - (h^2)/2 .* ∂C∂pₖ
    Jacᵏ⁺¹ₖ[nq+1:2nq,nq+1:2nq]     .= -(h^2)/2 .* ∂C∂pₖ
    Jacᵏ⁺¹ₖ[   1:nq ,2nq+nλ+1:2nq+nλ+ns] .= -(h^2)/2 .* (∂C∂sₖ + 1/2 .*∂F∂s)
    Jacᵏ⁺¹ₖ[nq+1:2nq,2nq+nλ+1:2nq+nλ+ns] .= -(h^2)/2 .* (∂C∂sₖ + 1/2 .*∂F∂s)
    
    
    cstr_jacobian!(Aₖ,   structure, solver_state.state_k)
    cstr_forces_jacobian!(∂Aᵀλ∂q, structure, qₖ₊₁,λₘ)

    auxi_jacobian!(∂S∂q,∂S∂s,structure,solver_state.state_kp1)
    
    Jacᵏ⁺¹ₖ₊₁ .= 0.0
    Jacᵏ⁺¹ₖ₊₁[    1:nq ,    1:nq ]    .=  1/2 .*∂Mₘhq̇ₘ∂qₘ .+ Mₘ .-(h^2)/2 .*(1/2 .*∂F∂q .+1/h.*∂F∂q̇ .+ ∂C∂qₖ₊₁)
    Jacᵏ⁺¹ₖ₊₁[ nq+1:2nq,    1:nq ]    .= -1/2 .*∂Mₘhq̇ₘ∂qₘ .- Mₘ .-(h^2)/2 .*(1/2 .*∂F∂q .+1/h.*∂F∂q̇ .+ ∂C∂qₖ₊₁) .- mass_norm.*∂Aᵀλ∂q
    Jacᵏ⁺¹ₖ₊₁[2nq+1:2nq+nλ, 1:nq ]    .= -mass_norm.*Aₖ₊₁
    Jacᵏ⁺¹ₖ₊₁[    1:nq,  nq+1:2nq]    .=  - (h^2)/2 .* ∂C∂pₖ₊₁
    Jacᵏ⁺¹ₖ₊₁[ nq+1:2nq, nq+1:2nq]    .=  h.*I(nq) .- (h^2)/2 .* ∂C∂pₖ₊₁
    Jacᵏ⁺¹ₖ₊₁[    1:nq ,2nq+1:2nq+nλ] .= -mass_norm.*transpose(Aₖ)
    Jacᵏ⁺¹ₖ₊₁[ nq+1:2nq,2nq+1:2nq+nλ] .= -mass_norm.*transpose(Aₖ₊₁)
    Jacᵏ⁺¹ₖ₊₁[   1:nq ,2nq+nλ+1:2nq+nλ+ns] .= -(h^2)/2 .*(1/2 .*∂F∂s)
    Jacᵏ⁺¹ₖ₊₁[nq+1:2nq,2nq+nλ+1:2nq+nλ+ns] .= -(h^2)/2 .*(1/2 .*∂F∂s)
    Jacᵏ⁺¹ₖ₊₁[2nq+nλ+1:2nq+nλ+ns,        1:nq ]       .= ∂S∂q
    Jacᵏ⁺¹ₖ₊₁[2nq+nλ+1:2nq+nλ+ns, 2nq+nλ+1:2nq+nλ+ns] .= ∂S∂s

    @. Jacᵏ⁺¹ₘu[    1: nq,1:nu] = -h^2/2*∂Fₘ∂u
    @. Jacᵏ⁺¹ₘu[ nq+1:2nq,1:nu] = -h^2/2*∂Fₘ∂u
    Jacᵏ⁺¹ₘu[2nq+1:end,1:nu] .= 0.0

    @. Jacᵏ⁺¹ₘc[    1: nq,1:nc] = -h^2/2*∂Fₘ∂c
    @. Jacᵏ⁺¹ₘc[ nq+1:2nq,1:nc] = -h^2/2*∂Fₘ∂c
    Jacᵏ⁺¹ₘc[2nq+1:end,1:nc] .= 0.0
end

function populate!(solver_state::Zhong06SolverState{<:MonoContactCoordinatesState},x,structure,jacobian_workspace,consts::Zhong06Constants)
    (;nq,nλ,mass_norm,ns) = consts
    solver_state.qₖ₊₁ .= x[1:nq]
    solver_state.λₘ   .= x[nq+1:nq+nλ]
    # solver_state.sₖ₊₁ .= x[nq+nλ+1:end]
    (;state_k,state_kp1,dt) = solver_state
    h = dt
    (;Mₘ,M⁻¹ₘ,Aₖ,Aₖ₊₁) = jacobian_workspace
    qₖ = state_k.q
    qₖ₊₁ = state_kp1.q
    pₖ = state_k.p
    pₖ₊₁ = state_kp1.p
    q̇ₖ₊₁ = state_kp1.q̇
    λₘ = state_kp1.λ
    cstr_jacobian!(Aₖ₊₁, structure, state_kp1)
    pₖ₊₁ .= -pₖ.+2/h.*Mₘ*(qₖ₊₁.-qₖ) .+ mass_norm/h.*(transpose(Aₖ₊₁)-transpose(Aₖ))*λₘ
    q̇ₖ₊₁ .= M⁻¹ₘ*pₖ₊₁
end

# ============================================================================
# Base primal solver functions
# ============================================================================

function guess_newton!(x,solver_state::Zhong06SolverState{<:PresFreeCoordinatesState},consts::Zhong06Constants, bot::Union{Nothing, Robot}=nothing)
    (;nq̌,nλ) = consts
    (;q̌ₖ,λₘ,sₖ,dt,q̌̇ₖ) = solver_state
    x[1:nq̌ ] .= q̌ₖ .+ dt.*q̌̇ₖ
    x[nq̌+1:nq̌+nλ] .= λₘ .= 0
    x[nq̌+nλ+1:end] .= sₖ
end

function populate!(solver_state::Zhong06SolverState{<:PresFreeCoordinatesState},x,structure,jacobian_workspace,consts::Zhong06Constants)
    (;nq̌,nλ,mass_norm) = consts
    solver_state.q̌ₖ₊₁ .= x[1:nq̌]
    solver_state.λₘ   .= x[nq̌+1:nq̌+nλ]
    solver_state.sₖ₊₁ .= x[nq̌+nλ+1:end]
    (;state_k,state_kp1,dt) = solver_state
    h = dt
    (;Ḿ,Aₖ,M,M̌⁻¹,M̄,Aₖ₊₁) = jacobian_workspace
    cstr_jacobian!(Aₖ₊₁, structure, state_kp1)
    qₖ = state_k.q
    p̌ₖ = state_k.p̌
    qₖ₊₁ = state_kp1.q
    p̌ₖ₊₁ = state_kp1.p̌
    pₖ₊₁ = state_kp1.p
    q̇ₖ₊₁ = state_kp1.q̇
    q̃̇ₖ₊₁ = state_kp1.q̃̇
    q̌̇ₖ₊₁ = state_kp1.q̌̇
    λₘ = state_kp1.λ
    p̌ₖ₊₁ .= -p̌ₖ.+2/h.*Ḿ*(qₖ₊₁.-qₖ) .+ mass_norm/h.*( transpose(Aₖ₊₁) - transpose(Aₖ) )* λₘ
    q̌̇ₖ₊₁ .= M̌⁻¹*( p̌ₖ₊₁.-M̄*q̃̇ₖ₊₁ )
    pₖ₊₁ .= M*q̇ₖ₊₁
end

function compute_constant_mass_residual!(Res, x, solver_state::Zhong06SolverState{<:PresFreeCoordinatesState},
        solver_cache,
        bot, policy, env
    )
    (;consts, jacobian_workspace) = solver_cache
    (;structure) = bot
    (;h, mass_norm, nq̌, nλ, ns, ) = consts
    (;qₖ₊₁, q̌ₖ₊₁, λₘ, sₖ₊₁, qₖ, q̌ₖ, qₘ, q̇ₘ, sₘ, tₘ, tₖ, sₖ) = solver_state
    (;Aₖ, Ḿ) = jacobian_workspace
    F̌ = solver_state.F̌
    p̌ₖ₊₁ = solver_state.p̌ₖ₊₁
    p̌ₖ = solver_state.p̌ₖ
    nx = nq̌ + nλ + ns
    q̌ₖ₊₁ .= x[   1:nq̌   ]
    λₘ   .= x[nq̌+1:nq̌+nλ]
    sₖ₊₁ .= x[nq̌+nλ+1:nx]
    
    qₘ .= (qₖ₊₁.+qₖ)./2
    q̇ₘ .= (qₖ₊₁.-qₖ)./h
    sₘ .= (sₖ₊₁ .+ sₖ)./2
    
    Aᵀₖ = transpose(Aₖ)
    # TODO: update p̌ₖ₊₁ for iLQR to work
    control!(bot, policy, solver_cache, solver_state)
    gen_force!(solver_state.state_mid, bot, env.field, policy, )
    
    Res[   1:nq̌   ] .= Ḿ*(qₖ₊₁.-qₖ) .-
                       h.*p̌ₖ .-
                       (h^2)/2 .*F̌ .-
                       mass_norm.*Aᵀₖ*λₘ
    ϕ = @view Res[nq̌+1:nq̌+nλ]
    cstr_function!(ϕ, structure,solver_state.state_kp1)
    ϕ .*= -mass_norm
    S = @view Res[nq̌+nλ+1:nx]
    auxi_function!(S, structure,solver_state.state_mid)
end

function compute_constant_mass_jacobian!(Jac, x, solver_state::Zhong06SolverState{<:PresFreeCoordinatesState},
        solver_cache,
        bot, policy, env,
    )
    (;consts, jacobian_workspace) = solver_cache
    (;structure) = bot
    (;field) = env
    (;h, mass_norm, nq̌, nλ, ns) = consts
    nx = nq̌ + nλ + ns
    (;qₖ₊₁, q̌ₖ₊₁, λₘ, sₖ₊₁, qₖ, q̌ₖ, qₘ, q̇ₘ, sₘ, tₘ, tₖ, sₖ) = solver_state
    
    # Unpack Jacobian workspace
    (;∂F∂s, ∂S∂s, M̌, Ḿ, Aₖ, Aₖ₊₁, ∂Aᵀλ∂q) = jacobian_workspace
    # Map to local names - use views as workspace matrices are now size nq x nq
    ∂F∂q̌ = view(jacobian_workspace.∂F∂q, 1:nq̌, 1:nq̌)
    ∂F∂q̌̇ = view(jacobian_workspace.∂F∂q̇, 1:nq̌, 1:nq̌)
    ∂F∂u = view(jacobian_workspace.∂Fₘ∂u, 1:nq̌, :)
    ∂S∂q̌ = view(jacobian_workspace.∂S∂q, :, 1:nq̌)
    ∂C∂q̌ = view(jacobian_workspace.∂C∂qₖ₊₁, 1:nq̌, 1:nq̌)
    ∂C∂p̌ = view(jacobian_workspace.∂C∂pₖ₊₁, 1:nq̌, 1:nq̌)
    q̌ₖ₊₁ .= x[1:nq̌]
    qₘ .= (qₖ₊₁.+qₖ)./2
    q̇ₘ .= (qₖ₊₁.-qₖ)./h
    sₘ .= (sₖ₊₁ .+ sₖ)./2
    
    Aᵀₖ = transpose(Aₖ)
    # TODO: update p̌ₖ₊₁ for iLQR to work
    # TODO, this is not correct, we need to use jacobian_workspace, not its views
    control_jacobian!(jacobian_workspace, bot, policy, solver_cache, solver_state)
    ∂p∂q = 2/h.*M̌  .+ mass_norm/h.*∂Aᵀλ∂q
    ∂C∂q̌ .+= ∂C∂p̌*∂p∂q
    #note: also need ∂p∂λ probably
    gen_force_state_jacobian!(∂F∂q̌, ∂F∂q̌̇, ∂F∂u, bot, field, policy,solver_state.state_mid,∂F∂s)
    auxi_jacobian!(∂S∂q̌,∂S∂s,structure,solver_state.state_mid)
    cstr_jacobian!(Aₖ₊₁, structure, solver_state.state_kp1)
    Jac .= 0.0
    Jac[      1:nq̌,         1:nq̌   ] .= M̌.-(h^2)/2 .*(1/2 .*∂F∂q̌.+1/h.*∂F∂q̌̇ .+ ∂C∂q̌) 
    Jac[   nq̌+1:nq̌+nλ,      1:nq̌   ] .= -mass_norm.*Aₖ₊₁
    Jac[nq̌+nλ+1:nx,         1:nq̌   ] .= 1/2 .* ∂S∂q̌
    Jac[      1:nq̌,      nq̌+1:nq̌+nλ] .= -mass_norm.*Aᵀₖ
    Jac[      1:nq̌,   nq̌+nλ+1:nx   ] .= -(h^2)/2 .*(1/2 .*∂F∂s)
    Jac[nq̌+nλ+1:nx,   nq̌+nλ+1:nx   ] .=             1/2 .*∂S∂s
end

function guess_newton!(x,solver_state::Zhong06SolverState{<:CoordinatesState},consts::Zhong06Constants, bot::Union{Nothing, Robot}=nothing)
    (;nq,nλ) = consts
    (;qₖ,q̇ₖ,dt,λₘ,sₖ) = solver_state
    
    if isnothing(bot)
        # Fallback to linear extrapolation
        x[1:nq] .= qₖ .+ dt.*q̇ₖ
    else
        # Use SE(3) enhanced guess
        se3_guess!(x, solver_state, bot, dt)
    end

    x[nq+1:nq+nλ] .= λₘ .= 0
    x[nq+nλ+1:end] .= sₖ
end

function se3_guess!(x, solver_state, bot, dt)
    (;structure) = bot
    (;bodies, connectivity) = structure
    
    update_bodies!(structure, solver_state.state_k)
    
    foreach(bodies) do body
        bodyid = body.prop.id
        sys_full_idx = connectivity.bodyid2sys_full_coords[bodyid]

        r = body.state.origin_frame.position
        v = body.state.origin_frame.velocity

        R = body.state.origin_frame.axes.X
        ω = body.state.origin_frame.angular_velocity
        Ω = body.state.origin_frame.local_angular_velocity

        # 3. Advance in SE(3)
        r_next = r + dt * v
        
        # Rotation advance: R_next = exp(ω * dt) * R (spatial rotation)
        # Using Rotations.jl's RotationVec for exponential map
        # rot_step = RotationVec(Tuple(ω * 0.0)...)
        # R_next = rot_step * R
        R_next = R
        # 4. Convert back to Natural Coordinates
        # Reconstruct axes from rotation matrix
        axes_next = Axes(SMatrix{size(R_next)...}(R_next))
        
        origin_frame = (
            position = SVector(r_next),
            velocity = SVector(v),
            axes = axes_next,
            angular_velocity = SVector(ω)
        )
        
        q_next, _ = cartesian_frame2coords(body.coords, origin_frame)
        
        # Update the guess vector x
        x[sys_full_idx] .= q_next
    end
end

function populate!(solver_state::Zhong06SolverState{<:CoordinatesState},x,structure,jacobian_workspace,consts::Zhong06Constants)
    (;nq,nλ,ns,mass_norm) = consts
    copyto!(solver_state.qₖ₊₁, 1, x, 1, nq)
    copyto!(solver_state.λₘ, 1, x, nq+1, nλ)
    copyto!(solver_state.sₖ₊₁, 1, x, nq+nλ+1, ns)
    (;state_k,state_kp1,dt) = solver_state
    h = dt
    (;Mₘ,Aₖ,M̌⁻¹,Aₖ₊₁) = jacobian_workspace

    cstr_jacobian!(Aₖ₊₁, structure, state_kp1)
    qₖ = state_k.q
    qₖ₊₁ = state_kp1.q
    pₖ = state_k.p
    pₖ₊₁ = state_kp1.p
    q̇ₖ₊₁ = state_kp1.q̇
    λₘ = state_kp1.λ
    
    @. q̇ₖ₊₁ = qₖ₊₁ - qₖ
    @. pₖ₊₁ = -pₖ
    mul!(pₖ₊₁, Mₘ, q̇ₖ₊₁, 2/h, 1)
    fill!(q̇ₖ₊₁, zero(eltype(q̇ₖ₊₁)))
    mul!(q̇ₖ₊₁, transpose(Aₖ₊₁), λₘ)
    mul!(q̇ₖ₊₁, transpose(Aₖ), λₘ, -1, 1)
    scale = mass_norm/h
    @. pₖ₊₁ += scale * q̇ₖ₊₁
    
    if has_constant_mass_matrix(structure) === Val(:false)
        update_bodies!(structure, solver_state.state_kp1)
        assemble_M⁻¹!(M̌⁻¹,structure)
    end
    mul!(q̇ₖ₊₁, M̌⁻¹, pₖ₊₁)
    # @show qₖ₊₁ M̌⁻¹ Mₘ pₖ pₖ₊₁ λₘ  2/h.*Mₘ*(qₖ₊₁.-qₖ) (transpose(Aₖ₊₁)-transpose(Aₖ))*λₘ
end

function precompute!(jacobian_workspace::Zhong06JacobianWorkspace,solver_state::Zhong06SolverState,structure)
    (;Aₖ) = jacobian_workspace
    cstr_jacobian!(Aₖ, structure, solver_state.state_k)
end

function compute_constant_mass_residual!(Res, x, solver_state::Zhong06SolverState{<:CoordinatesState},
        solver_cache,
        bot, policy, env
    )
    (;consts, jacobian_workspace) = solver_cache
    (;structure) = bot
    (;h, mass_norm, nq, nλ, ns, ) = consts
    (;qₖ₊₁, λₘ, sₖ₊₁, qₖ, qₘ, q̇ₘ, sₘ, tₘ, tₖ, sₖ) = solver_state
    (;Aₖ, Mₘ , ϕbuf) = jacobian_workspace
    F = solver_state.state_mid.F
    pₖ₊₁ = solver_state.pₖ₊₁
    pₖ = solver_state.pₖ
    nx = nq + nλ + ns
    qₖ₊₁ .= @view x[   1:nq   ]
    λₘ   .= @view x[nq+1:nq+nλ]
    sₖ₊₁ .= @view x[nq+nλ+1:nx]

    @. qₘ = (qₖ₊₁ + qₖ) /2
    @. q̇ₘ = (qₖ₊₁ - qₖ) /h
    @. sₘ = (sₖ₊₁ + sₖ) /2

    Aᵀₖ = transpose(Aₖ)
    # note: update pₖ₊₁ for iLQR to work
    populate!(solver_state,x,structure,jacobian_workspace,consts)
    control!(bot, policy, solver_cache, solver_state)
    gen_force!(solver_state.state_mid, bot, env.field, policy)

    Resq = @view Res[1:nq]
    tmp = q̇ₘ
    @. tmp = qₖ₊₁ - qₖ
    mul!(Resq, Mₘ, tmp)                 # Resq = Mₘ*(qₖ₊₁ - qₖ)
    mul!(tmp, Aᵀₖ, λₘ)                 # tmp = Aᵀₖ * λₘ
    @. Resq -= h * pₖ + (h^2)/2 * F +  mass_norm * tmp                  # Resq .-= h .* pₖ + (h^2)/2 .* F + mass_norm .* tmp
    ϕ = @view Res[nq+1:nq+nλ]
    cstr_function!(ϕbuf, structure, solver_state.state_kp1)
    @. ϕ = -mass_norm * ϕbuf
    S = @view Res[nq+nλ+1:nx] 
    auxi_function!(S, structure,solver_state.state_kp1)
end

function compute_constant_mass_jacobian!(Jac, x, solver_state::Zhong06SolverState{<:CoordinatesState},
        solver_cache,
        bot, policy, env,
    )
    (;consts, jacobian_workspace) = solver_cache
    (;structure) = bot
    (;field) = env
    (;h, mass_norm, nq, nλ, ns) = consts
    nx = nq + nλ + ns
    (;qₖ₊₁, λₘ, sₖ₊₁, qₖ, qₘ, q̇ₘ, sₘ, tₘ, tₖ, sₖ) = solver_state

    # Unpack Jacobian workspace
    (;∂F∂s, ∂S∂s, Aₖ, Aₖ₊₁, ∂Aᵀλ∂q) = jacobian_workspace
    Mₘ = jacobian_workspace.Mₘ
    
    # Map to local names (CoordinateState uses full coords)
    ∂F∂q = jacobian_workspace.∂F∂q
    ∂F∂q̇ = jacobian_workspace.∂F∂q̇
    ∂F∂u = jacobian_workspace.∂Fₘ∂u
    ∂S∂q = jacobian_workspace.∂S∂q
    ∂C∂qₖ₊₁ = jacobian_workspace.∂C∂qₖ₊₁
    ∂C∂pₖ₊₁ = jacobian_workspace.∂C∂pₖ₊₁

    qₖ₊₁ .= x[1:nq]
    qₘ .= (qₖ₊₁.+qₖ)./2
    q̇ₘ .= (qₖ₊₁.-qₖ)./h
    sₘ .= (sₖ₊₁ .+ sₖ)./2

    Aᵀₖ = transpose(Aₖ)
    # note: update pₖ₊₁ for iLQR to work
    populate!(solver_state,x,structure,jacobian_workspace,consts)

    control_jacobian!(jacobian_workspace, bot, policy, solver_cache, solver_state)
    cstr_forces_jacobian!(∂Aᵀλ∂q, structure,qₖ₊₁,λₘ)
    ∂p∂q = 2/h.*Mₘ  .+ mass_norm/h.*∂Aᵀλ∂q
    ∂C∂qₖ₊₁ .+= ∂C∂pₖ₊₁*∂p∂q
    #note: also need ∂p∂λ ∂p∂s probably
    gen_force_state_jacobian!(∂F∂q, ∂F∂q̇, ∂F∂u, bot, field, policy,solver_state.state_mid,∂F∂s)
    auxi_jacobian!(∂S∂q,∂S∂s,structure,solver_state.state_kp1)
    cstr_jacobian!(Aₖ₊₁, structure, solver_state.state_kp1)
    Jac .= 0.0
    Jac[      1:nq,         1:nq   ] .= Mₘ.-(h^2)/2 .*(1/2 .*∂F∂q.+1/h.*∂F∂q̇ .+ ∂C∂qₖ₊₁)
    Jac[   nq+1:nq+nλ,      1:nq   ] .= -mass_norm.*Aₖ₊₁
    Jac[nq+nλ+1:nx,         1:nq   ] .=  ∂S∂q
    Jac[      1:nq,      nq+1:nq+nλ] .= -mass_norm.*Aᵀₖ
    Jac[      1:nq,   nq+nλ+1:nx   ] .= -(h^2)/2 .*(1/2 .*∂F∂s)
    Jac[nq+nλ+1:nx,   nq+nλ+1:nx   ] .=  ∂S∂s
end

# ============================================================================
# Inner Solver Helper Functions for Inner CCP
# ============================================================================

function compute_inner_contact_linearization!(solver_cache::Zhong06_CCP_Constant_Mass_Inner_Cache,solver_state::Zhong06SolverState,
        contact_cache,contact_workspace
    )
    (;bot, jacobian_workspace, consts) = solver_cache
    (;structure) = bot
    (;Mₘ, M⁻¹ₘ, Aₖ, Aₖ₊₁, ∂Aᵀλ∂q) = jacobian_workspace
    (;h, mass_norm, nq, nλ) = consts
    (;qₖ, pₖ, qₖ₊₁, pₖ₊₁, λₘ, q̇ₖ, vₖ, vₖ₊₁) = solver_state
    (;na, H, activated_restitution_coefficients, D, L) = contact_cache
    Dₘ = contact_cache.Dper
    Dₖ = contact_cache.Dimp
    (; state_k, state_kp1, state_mid) = solver_state
    (;𝐁, 𝐛, 𝐜ᵀ, 𝐲) = contact_workspace
    
    nΛ = 3na
    
    # Update B
    𝐁 .= 0
    # Note: B is (nx, nΛ). Top block is (nq, nΛ).
    # 𝐁[1:nq, 1:nΛ] .= h .* mass_norm .* transpose(D) * H
    # But wait, CCP_inner original code: 𝐁[1:nq,1:nΛ] .= h.*mass_norm.*transpose(D)*H
    # D is (3na, nq). transpose(D) is (nq, 3na). H is (3na, 3na)? No, H is usually selection matrix.
    # Let's assume dimensions match.
    
    𝐁[1:nq, :] .= h .* mass_norm .* transpose(D) * H
    
    # Compute linearization terms
    # pₖ is already updated in solver_state before calling this?
    # In CCP_inner original, pₖ is updated inside ns_stepk! if doin=true.
    # Here we assume pₖ is up to date in solver_state.
    
    # vₖ = M⁻¹ₘ * pₖ
    # vₖ is already in solver_state (q̇ₖ).
    
    # Let's use local vars to match math
    qₖ = qₖ
    pₖ = pₖ
    qₖ₊₁ = qₖ₊₁
    λₘ = λₘ

    cstr_jacobian!(Aₖ, structure, solver_state.state_k)
    cstr_jacobian!(Aₖ₊₁, structure, solver_state.state_kp1)

    # Derivatives
    # ∂vₘ∂q_curr = 1/h * I
    # ∂v_curr∂q_curr = 2/h * I + mass_norm/h * M⁻¹ₘ * ∂Aᵀλ∂q(qₖ₊₁, λₘ)
    cstr_forces_jacobian!(∂Aᵀλ∂q, structure, qₖ₊₁, λₘ)
    
    ∂v_curr∂q_curr = 2/h * I + mass_norm/h .* M⁻¹ₘ * ∂Aᵀλ∂q
    ∂v_curr∂λ_curr = mass_norm .* M⁻¹ₘ * transpose(Aₖ₊₁ - Aₖ) / h
    
    # vₘ = (qₖ₊₁ - qₖ) / h
    vₘ = (qₖ₊₁ .- qₖ) ./ h
    # v_curr = M⁻¹ₘ * p_curr
    # p_curr is computed in compute_inner_residual! usually.
    # But here we need it for linearization.
    # Momentum_k(qₖ, pₖ, qₖ₊₁, λₘ, Mₘ, A, mass_norm, h)
    p_curr  = -pₖ.+2/h.*Mₘ*(qₖ₊₁.-qₖ) .+ 
        mass_norm/h.*(transpose(Aₖ₊₁ - Aₖ))*λₘ

    v_curr = M⁻¹ₘ * p_curr
    
    v́⁺ = Dₘ * vₘ .+ Dₖ * v_curr
    ∂v́⁺∂q_curr = Dₘ * (1/h * I) .+ Dₖ * ∂v_curr∂q_curr
    
    𝐜ᵀ .= 0
    
    v́_prev = Dₖ * (M⁻¹ₘ * pₖ) # v_prev = M⁻¹ₘ * pₖ
    
    for i = 1:na
        is = 3(i-1)
        vⁱ_prev = @view v́_prev[is+1:is+3]
        vⁱ⁺ = @view v́⁺[is+1:is+3]
        vₙⁱ_prev = vⁱ_prev[1]
        
        v́ₜⁱ = norm(vⁱ⁺[2:3]) + activated_restitution_coefficients[i] * min(vₙⁱ_prev, 0)
        
        𝐛[is+1:is+3] .= [v́ₜⁱ, 0, 0]
        
        # Derivatives for C
        # 𝐜ᵀ[is+1, 1:nq] .= ...
        term = (v́⁺[is+2]*∂v́⁺∂q_curr[is+2,:] .+ v́⁺[is+3]*∂v́⁺∂q_curr[is+3,:]) ./ (norm(v́⁺[is+2:is+3]) + 1e-14)
        𝐜ᵀ[is+1, 1:nq] .= term
        𝐜ᵀ[is+1:is+3, 1:nq] .+= ∂v́⁺∂q_curr[is+1:is+3, :]
        
        # w.r.t λ
        Dⁱₖ = Dₖ[is+1:is+3, :]
        𝐜ᵀ[is+1:is+3, nq+1:nq+nλ] .= Dⁱₖ * ∂v_curr∂λ_curr
    end
    
    # 𝐲 .= (v́⁺ + 𝐛)
    𝐲 .= v́⁺ .+ 𝐛
    
    return nothing
end

function compute_inner_residual!(workspace::NewtonWorkspace,solver_state::Zhong06SolverState,
        solver_cache::Zhong06_CCP_Constant_Mass_Inner_Cache,
        contact_cache, contact_workspace
    )
    (;bot, policy, field, consts, jacobian_workspace,) = solver_cache
    (;structure) = bot
    (;h, mass_norm, nq, nλ, ns) = consts
    (;Mₘ, Aₖ) = jacobian_workspace
    (;qₖ, pₖ, qₖ₊₁, pₖ₊₁, λₘ, sₖ, sₖ₊₁, qₘ, q̇ₘ, sₘ, tₘ) = solver_state
    (;na) = contact_cache
    (;𝐁, Λ) = contact_workspace
    
    # Unpack workspace
    𝐫𝐞𝐬 = workspace.Res
    x = workspace.x
    
    n1 = nq
    n2 = nq + nλ
    nx = nq + nλ + ns
    
    # Extract state variables from x
    qₖ₊₁ .= x[1:n1]
    λₘ .= x[n1+1:n2]
    sₖ₊₁ .= x[n2+1:nx]
    
    # Compute midpoint quantities
    qₘ .= (qₖ₊₁ .+ qₖ) ./ 2
    q̇ₘ .= (qₖ₊₁ .- qₖ) ./ h
    sₘ .= (sₖ₊₁ .+ sₖ) ./ 2
    
    # Compute forces
    gen_force!(solver_state.state_mid, bot, field, policy, )
    Fₘ = solver_state.state_mid.F
    
    cstr_jacobian!(Aₖ, structure, solver_state.state_k)
    
    # Residuals
    # 𝐫𝐞𝐬[1:nq] .= -h.*pₖ .+ Mₘ*(qₖ₊₁.-qₖ) .- mass_norm.*transpose(Aₖ)*λₘ .- (h^2)/2 .*Fₘ
    𝐫𝐞𝐬[1:nq] .= -h .* pₖ .+ Mₘ * (qₖ₊₁ .- qₖ) .-
                  mass_norm .* transpose(Aₖ) * λₘ .-
                  (h^2) / 2 .* Fₘ
                  
    # Note: CCP_inner uses qₖ (current guess) which is qₖ₊₁ in solver_state
    ϕ = @view 𝐫𝐞𝐬[nq+1:nq+nλ] 
    cstr_function!(ϕ,structure,solver_state.state_kp1)
    ϕ.*= -mass_norm 
    
    # Aux Constraints
    S = @view 𝐫𝐞𝐬[nq+nλ+1:nx]
    auxi_function!(S, structure,solver_state.state_mid)
    
    # Contact contribution
    if na != 0
        # 𝐁 is (nx, nΛ). Λ is (nΛ).
        𝐫𝐞𝐬 .-= 𝐁 * Λ
    end
    
    return nothing
end

function compute_inner_jacobian!(workspace::NewtonWorkspace,solver_state::Zhong06SolverState,
        solver_cache::Zhong06_CCP_Constant_Mass_Inner_Cache,contact_cache
    )
    (;bot, policy, field, consts, jacobian_workspace,) = solver_cache
    (;structure) = bot
    (;h, mass_norm, nq, nλ, ns) = consts
    (;Mₘ, ∂F∂q, ∂F∂q̇, ∂F∂s, ∂Fₘ∂u, ∂S∂q, ∂S∂s, Aₖ, Aₖ₊₁) = jacobian_workspace
    (;qₖ, qₖ₊₁, qₘ, q̇ₘ, sₘ, tₘ) = solver_state
    
    𝐉 = workspace.Jac
    
    n1 = nq
    n2 = nq + nλ
    nx = nq + nλ + ns
    
    # Compute Jacobians
    gen_force_state_jacobian!(∂F∂q, ∂F∂q̇, ∂Fₘ∂u, bot, field, policy, solver_state.state_mid, ∂F∂s)
    
    auxi_jacobian!(∂S∂q, ∂S∂s, structure, solver_state.state_mid)
    
    cstr_jacobian!(Aₖ, structure, solver_state.state_k)
    cstr_jacobian!(Aₖ₊₁, structure, solver_state.state_kp1)
    
    𝐉 .= 0.0
    
    𝐉[1:nq, 1:nq] .= Mₘ .- (h^2)/2 .* (1/2 .* ∂F∂q .+ (1/h) .* ∂F∂q̇)
    
    𝐉[1:nq, nq+1:nq+nλ] .= -mass_norm .* transpose(Aₖ)
    
    𝐉[1:nq, nq+nλ+1:nx] .= -(h^2)/2 .* (1/2 .* ∂F∂s)
    
    𝐉[nq+1:nq+nλ, 1:nq] .= -mass_norm .* Aₖ₊₁
    
    𝐉[nq+nλ+1:nx, 1:nq] .= ∂S∂q
    𝐉[nq+nλ+1:nx, nq+nλ+1:nx] .= ∂S∂s
    
    return nothing
end

function populate!(solver_state::Zhong06SolverState{<:InnerContactCoordinatesState},x,structure,jacobian_workspace,consts::Zhong06Constants)
    (;nq,nλ,mass_norm) = consts
    solver_state.qₖ₊₁ .= x[1:nq]
    solver_state.λₘ   .= x[nq+1:nq+nλ]
    solver_state.sₖ₊₁ .= x[nq+nλ+1:end]
    (;state_k,state_kp1,dt) = solver_state
    h = dt
    (;Mₘ,M⁻¹ₘ,Aₖ,Aₖ₊₁) = jacobian_workspace
    qₖ = state_k.q
    qₖ₊₁ = state_kp1.q
    pₖ = state_k.p
    pₖ₊₁ = state_kp1.p
    q̇ₖ₊₁ = state_kp1.q̇
    λₘ = state_kp1.λ
    cstr_jacobian!(Aₖ₊₁, structure, state_kp1)
    pₖ₊₁ .= -pₖ.+2/h.*Mₘ*(qₖ₊₁.-qₖ) .+ mass_norm/h.*(transpose(Aₖ₊₁)-transpose(Aₖ))*λₘ
    q̇ₖ₊₁ .= M⁻¹ₘ*pₖ₊₁
end


function compute_nonconstant_mass_residual!(Res, x, solver_state::Zhong06SolverState,
        solver_cache::Zhong06_Nonconstant_Mass_Cache
    )
    (;
        bot, policy, env,
        consts,
        jacobian_workspace,
    ) = solver_cache
    (; structure) = bot
    (; h, mass_norm, nq, nλ) = consts
    (; Mₘ, Aₖ, Fₘ, ϕbuf) = jacobian_workspace

    pack_nonconstant_state!(solver_state, x, consts)

    update_bodies!(structure, solver_state.state_mid)
    assemble_M!(Mₘ, structure)
    control!(bot, policy, solver_cache, solver_state)
    gen_force!(solver_state.state_mid, bot, env.field, policy)
    Fₘ .= solver_state.state_mid.F

    cstr_jacobian!(Aₖ, structure, solver_state.state_k)

    Resq = @view Res[1:nq]
    Resq .= Mₘ * (solver_state.qₖ₊₁ .- solver_state.qₖ) .-
            h .* solver_state.pₖ .-
            (h^2) / 2 .* Fₘ .-
            mass_norm .* transpose(Aₖ) * solver_state.λₘ
    ϕ = @view Res[nq+1:nq+nλ]
    cstr_function!(ϕbuf, structure, solver_state.state_kp1)
    @. ϕ = -mass_norm * ϕbuf
    nothing
end

function compute_nonconstant_mass_jacobian!(Jac, x, solver_state::Zhong06SolverState,
        solver_cache::Zhong06_Nonconstant_Mass_Cache
    )
    (;
        bot, policy, env,
        consts,
        jacobian_workspace,
    ) = solver_cache
    (; structure) = bot
    (; h, mass_norm, nq, nλ) = consts
    (;
        Mₘ,
        ∂Mₘhq̇ₘ∂qₘ,
        Aₖ,
        Aₖ₊₁,
        ∂F∂q,
        ∂F∂q̇,
        ∂Fₘ∂u,
        ∂F∂s,
    ) = jacobian_workspace

    pack_nonconstant_state!(solver_state, x, consts)

    update_bodies!(structure, solver_state.state_mid)
    assemble_M!(Mₘ, structure)
    assemble_∂Mq̇∂q!(∂Mₘhq̇ₘ∂qₘ, structure)
    ∂Mₘhq̇ₘ∂qₘ .*= h

    gen_force_state_jacobian!(
        ∂F∂q, ∂F∂q̇, ∂Fₘ∂u, bot, env.field, policy, solver_state.state_mid, ∂F∂s
    )
    cstr_jacobian!(Aₖ, structure, solver_state.state_k)
    cstr_jacobian!(Aₖ₊₁, structure, solver_state.state_kp1)

    Jac .= 0.0
    Jac[1:nq, 1:nq] .= Mₘ .+ 1 / 2 .* ∂Mₘhq̇ₘ∂qₘ .-
                       (h^2) / 2 .* (1 / 2 .* ∂F∂q .+ (1 / h) .* ∂F∂q̇)
    Jac[1:nq, nq+1:nq+nλ] .= -mass_norm .* transpose(Aₖ)
    Jac[nq+1:nq+nλ, 1:nq] .= -mass_norm .* Aₖ₊₁
    Jac[nq+1:nq+nλ, nq+1:nq+nλ] .= 0.0
    nothing
end