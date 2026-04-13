
"""
    gen_force!(inst_state::AbstractCoordinatesState, bot::Robot, field::AbstractField, policy::AbstractPolicy, ; )

Assembles the total generalized force vector `F` acting on the robot given the configuration
encoded by `inst_state` (`q`, `q̇`, `t`, possibly `s`), the robot model `bot`, policy, and external field.

This method performs the following steps, in order:
1. Clears all internal force accumulators on the structure.
2. Applies actuation/controls according to `policy` and `inst_state`.
3. Updates body kinematics from the current state.
4. Updates any auxiliary apparatuses (possibly parameterized by `s`).
5. Applies the external field to the structure (e.g., gravity, magnetics).
6. Assembles all forces and stores the result in `F`.

All operations are done in-place. The result in `inst_state.F` reflects the *total* forces
resulting from both internal (structural, actuation) and external (field) contributions.

# Arguments
- `inst_state`: Pre-allocated `CoordinatesState` to store the force output.
- `bot`: Robot model containing structure and actuation.
- `field`: External field object (e.g., gravity, electromagnetics).
- `policy`: Actuation or control policy.
- (optional keyword arguments to pass through.)

# Example
```julia
gen_force!(inst_state, bot, field, mypolicy, )
```
Afterwards, `inst_state.F` will contain the total generalized forces for the system at this state.

"""
function gen_force!(inst_state::AbstractCoordinatesState, bot::Robot, field::AbstractField, policy::AbstractPolicy, ; )
    (;s) = inst_state
    (;structure) = bot
    clear_forces!(structure)
    actuate!(bot,policy,inst_state)
    lazy_update_bodies!(structure,inst_state)
    update_apparatuses!(structure, s)
    apply_field!(structure, field; )
    assemble_forces!(inst_state, structure)
end

"""
    jacobian of generalized force with respect to state (q, q̇
    
This is for forward dynamics.
"""
function gen_force_state_jacobian!(∂F∂q, ∂F∂q̇, ∂F∂u, bot::Robot, field::AbstractField, policy::AbstractPolicy, inst_state::AbstractCoordinatesState, ∂F∂s=nothing)
    (;q,q̇,t,s) = inst_state
    (;structure,hub) = bot
    ∂F∂q .= 0
    ∂F∂q̇ .= 0
    ∂F∂u .= 0
    if s !== nothing
        ∂F∂s .= 0
    end
    clear_forces!(structure)
    actuate!(bot,policy,inst_state)
    lazy_update_bodies!(structure,inst_state)
    update_apparatuses!(structure, s)
    assemble_forces!(inst_state, structure)

    gen_force_state_jacobian!(∂F∂q, ∂F∂q̇, ∂F∂s, ∂F∂u, policy, bot, inst_state)
    gen_force_auxi_jacobian!(∂F∂s,structure,inst_state,)
    
    add_tangent_stiffness_matrix!(∂F∂q, structure, field)
    add_tangent_damping_matrix!(∂F∂q̇, structure, field)
end

"""
    jacobian of generalized force with respect to state (q, q̇) and control (u)
    
This is for adjoint dynamics.
"""
function gen_force_jacobian!(∂F∂q, ∂F∂q̇, ∂F∂u, ∂F∂c, bot::Robot, field::AbstractField, policy::AbstractPolicy, inst_state::AbstractCoordinatesState, ∂F∂s=nothing)
    (;q,q̇,t,s) = inst_state
    (;structure,hub) = bot
    (;bodies, apparatuses) = structure
    (;actuators,coalition) = hub
    cnt = structure.connectivity
    (; num_of_params,apparid2params_idx) = cnt
    (; num_of_actions,actid2sys_actions_idx) = coalition
    
    structure.state.system.t = t
    structure.state.system.q .= q
    structure.state.system.q̇ .= q̇
    structure.state.system.s .= s
    
    ∂F∂q .= 0
    ∂F∂q̇ .= 0
    ∂F∂u .= 0
    if ∂F∂s !== nothing
        ∂F∂s .= 0.0
    end
    clear_forces!(structure)
    # Order matters: stretch! sets ALL apparatus params from structural params c,
    # then actuate! OVERRIDES the controlled params (e.g., rest lengths) from control u.
    # This ensures gradients work correctly for both c (structural) and θ (control).
    stretch!(structure)
    actuate!(bot,policy,inst_state)
    lazy_update_bodies!(structure,inst_state)
    update_apparatuses!(structure,s)
    apply_field!(structure, field)
    assemble_forces!(inst_state, structure)
    gen_force_state_jacobian!(∂F∂q,∂F∂q̇,∂F∂s,∂F∂u,policy,bot,inst_state)

    add_tangent_stiffness_matrix!(∂F∂q, structure, field)
    add_tangent_damping_matrix!(∂F∂q̇, structure, field)

    if ∂F∂s !== nothing
        gen_force_auxi_jacobian!(∂F∂s,structure,inst_state)
    end

    if ∂F∂c isa AbstractArray
        ∂F∂c .= 0.0
        # todo bodies's structural parameters, such as shape, inertia and elasticity 
        # foreach(bodies) do body
        #     c_idx = bodyid2sys_st_paras_idx[body.id]
        #     gen_force_para_jacobian!(
        #         (@view ∂F∂c[:,c_idx]),
        #         structure,
        #         body
        #     )
        # end
        foreach(apparatuses) do appar
            c_idx = apparid2params_idx[appar.id]
            gen_force_para_jacobian!(
                (@view ∂F∂c[:,c_idx]),
                appar,
                structure,
                q,q̇,t,
                (@view structure.state.system.c[c_idx]),
            )
        end
    end
end


function gen_force_state_jacobian!(∂F∂q,∂F∂q̇,∂F∂s,∂F∂u,policy::NoPolicy,bot::Robot,inst_state::AbstractCoordinatesState;)
	# do nothing
end

function gen_force_auxi_jacobian!(∂F∂s,structure::AbstractStructure,inst_state)
    (;q,q̇,t,s) = inst_state
    (;apparatuses,connectivity) = structure
    (;
        apparid2sys_aux_var_idx,
        apparid2sys_full_coords_idx
    ) = connectivity
    foreach(apparatuses) do appar
        appar_aux_var_idx = apparid2sys_aux_var_idx[appar.id]
        appar_full_coords_idx = apparid2sys_full_coords_idx[appar.id]
        gen_force_auxi_jacobian!(
            #note changed from (@view ∂F∂s[appar_full_coords_idx,appar_aux_var_idx]),
            #why this works?
            (@view ∂F∂s[:,appar_aux_var_idx]),
            appar,
            structure,
            q,q̇,t,
            (@view s[appar_aux_var_idx]),
        )
    end
end
