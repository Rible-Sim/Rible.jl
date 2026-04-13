"""
Rigid Body Interface Implementations

This module implements the body interface methods for rigid bodies.
"""



"""
Clear Rigid Body Forces and torques.
$(TYPEDSIGNATURES)
"""
function clear_forces!(body::AbstractRigidBody)
    (;mass_locus_state,loci_states,) = body.state
    mass_locus_state.force .= 0
    mass_locus_state.torque .= 0
    foreach(loci_states) do locus_state
        locus_state.force .= 0
        locus_state.torque .= 0
    end
end

function reset!(body::AbstractRigidBody)
    (;mass_locus_state,loci_states,) = body.state
    mass_locus_state.force .= 0
    mass_locus_state.torque .= 0
    foreach(loci_states) do locus_state
        locus_state.force .= 0
        locus_state.torque .= 0
        reset!(locus_state.contact_state)
    end
end

function body_state2coords_state(state::RigidBodyState,coords::AbstractCoordinates)
    (;origin_frame) = state
    cartesian_frame2coords(coords,origin_frame)
end

# Energy methods
"""
Return Rigid Body kinetic energy.
$(TYPEDSIGNATURES)
"""
function kinetic_energy(body::AbstractRigidBody{3,T}) where {T}
    (; mass, inertia) = body.prop
    (; mass_locus_state) = body.state
    v = mass_locus_state.frame.velocity
    Ω = to_3D(mass_locus_state.frame).local_angular_velocity
    1/2*mass*transpose(v)*v + 1/2 * transpose(Ω)*inertia*Ω
end

"""
Return Rigid Body kinetic energy.
$(TYPEDSIGNATURES)
"""
function kinetic_energy(body::AbstractRigidBody{2,T}) where {T}
    (; prop, state) = body
    (; mass, mass_locus, inertia) = prop
    o = zero(T)
    i = one(T)
    M = Diagonal(mass .* ones(3))
    c = vcat(
        mass_locus.position,
        o
    )
    MS = M * skew(c)

    I3 = SMatrix{3,3}(
        [
            inertia[1, 1] inertia[1, 2] o;
            inertia[2, 1] inertia[2, 2] o;
            o o inertia[1, 1]+inertia[2, 2]
        ]
    )

    inertia_tensor = [
        M -MS;
        MS I3+M*skew(c)*skew(c)'
    ]

    origin_frame_3d = CartesianFrame{3}(state.origin_frame)

    R = origin_frame_3d.axes.X
    Ω = origin_frame_3d.local_angular_velocity
    v = origin_frame_3d.velocity
    V = transpose(R) * v

    η = vcat(
        V,
        Ω
    )

    0.5η' * inertia_tensor * η
end

"""
Return Rigid Body gravity potential energy.
$(TYPEDSIGNATURES)
"""
function potential_energy_field(body::RigidBody, field::Gravity, q)
    (; mass) = body.prop
    (; mass_locus_state) = body.state
    gravity_acceleration = get_gravity(body)
    -transpose(mass_locus_state.frame.position) * gravity_acceleration * mass
end

"""
Return Rigid Body Strain potential energy.
$(TYPEDSIGNATURES)
"""
function potential_energy_strain(body::AbstractRigidBody)
    zero(get_numbertype(body))
end

# Field application
function apply_field!(F, body::RigidBody, field::Gravity, q; )
    gravity_acceleration = get_gravity(body)
    (; Cps, Cg, inertia_cache) = body.cache
    f_gra = gravity_acceleration * body.prop.mass
    mul!(F, transpose(Cg), f_gra, 1, 1)
end

function field_jacobian!(∂F∂q, body::AbstractBody, field::Gravity, q;)
end

function update_transformations!(coords::AbstractCoordinates,cache::AbstractBodyCache,state::AbstractFlexibleBodyState,prop::AbstractFlexibleBodyProperty, inst_state)
    # default no-op implementation
end

function stretch_loci!(coords::AbstractCoordinates,cache,prop::AbstractFlexibleBodyProperty, c)
    # to be implemented
end

function centrifugal_force!(F,state::RigidBodyState,cache::RigidBodyCache)
    (;mass_locus_state,loci_states) = state
    (;Cps,Cg,inertia_cache) = cache
    F .+= inertia_cache.∂T∂qᵀ
end

function mass_center_force!(F,state::RigidBodyState,cache::RigidBodyCache)
    (;mass_locus_state,loci_states) = state
    (;Cps,Cg,inertia_cache) = cache
    mul!(F,transpose(Cg),mass_locus_state.force,1,1)
end

function concentrated_force!(F,state::RigidBodyState,cache::RigidBodyCache)
    (;mass_locus_state,loci_states) = state
    (;Cps,Cg,inertia_cache) = cache
    for (pid,locus_state) in enumerate(loci_states)
        mul!(F,transpose(Cps[pid]),locus_state.force,1,1)
    end
end
