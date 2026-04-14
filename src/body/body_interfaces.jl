
"""
Body Interface

This module defines the interface for bodies in Rible.
All body implementations should implement these methods.
"""

# Utility functions for abstract types
function Base.isless(a::AbstractBody,b::AbstractBody)
    isless(a.prop.id,b.prop.id)
end

# Body property interface

"""
    get_id(body)

Get unique identifier for the body.
"""

# Interface method implementations
get_id(body::AbstractBody) = body.prop.id


# Default implementations for abstract bodies
update_inertia_cache!(body::AbstractBody,inst_state) = update_inertia_cache!(body.cache,body.coords,body.prop,inst_state)
lazy_update_state!(body::AbstractBody,inst_state) = lazy_update_state!(body.state,body.coords,body.cache,body.prop,inst_state)
update_state!(body::AbstractBody,inst_state) = update_state!(body.state,body.coords,body.cache,body.prop,inst_state)
update_loci_states!(body::AbstractBody,inst_state) = update_loci_states!(body.state,body.coords,body.cache,body.prop,inst_state)
stretch_loci!(body::AbstractBody,inst_state) = stretch_loci!(body.coords,body.cache,body.prop,inst_state)
update_transformations!(body::AbstractBody,inst_state) = update_transformations!(body.coords,body.cache,body.state,body.prop,inst_state)
kinetic_energy_coords(body::AbstractBody,inst_state) = kinetic_energy(body.coords,body.cache,body.state,body.prop,inst_state)
cstr_function!(ret, body::AbstractBody, inst_state::AbstractCoordinatesState) = cstr_function!(ret, body.coords, body.cache, inst_state)
cstr_jacobian!(ret, body::AbstractBody, inst_state::AbstractCoordinatesState) = cstr_jacobian!(ret, body.coords, body.cache, inst_state)
cstr_velocity_jacobian!(ret, body::AbstractBody, inst_state::AbstractCoordinatesState) = cstr_velocity_jacobian!(ret, body.coords, body.cache, inst_state)
cstr_velocity_jacobian!(ret, body::AbstractBody, q̇::AbstractVector) = cstr_velocity_jacobian!(ret, body.coords, body.cache, q̇)
cstr_forces_jacobian!(ret, body::AbstractBody, inst_state::AbstractCoordinatesState) = cstr_forces_jacobian!(ret, body.coords, body.cache, inst_state)


# Default implementations for Coordinates
function cstr_function!(ret, coords::AbstractCoordinates, cache::AbstractBodyCache, q,d = get_deform(coords))
    #todo make use cache.jac_cache as ret_full
    ret_full = zeros(eltype(ret),get_num_of_intrinsic_cstr(coords))
    cstr_function!(ret_full, coords, q ,d)
    ret .= ret_full[get_cstr_idx(coords)]
end

function cstr_jacobian!(ret, coords::AbstractCoordinates, cache::AbstractBodyCache, inst_state::AbstractCoordinatesState)
    #todo make use cache.jac_cache as ret_full
    ret_full = zeros(eltype(ret),get_num_of_intrinsic_cstr(coords),get_num_of_coords(coords))
    cstr_jacobian!(ret_full, coords, cache.jac_cache, inst_state.q)
    ret .= ret_full[get_cstr_idx(coords),:]
end


# Default implementations for body interface methods
function potential_energy_field(body::AbstractBody, field::NoField, q)
    zero(eltype(q))
end

function apply_field!(F, body::AbstractBody, field::NoField, q;)
end

function field_jacobian!(∂F∂q, body::AbstractBody, field::NoField, q;)
end

function body_state2coords_state(body::AbstractBody)
    body_state2coords_state(body.state,body.coords)
end

function strain!(F,prop::AbstractBodyProperty,coords,state::AbstractBodyState,cache::AbstractBodyCache) end


"""
    get_mass(body)

Get mass of the body.
"""
function get_mass end

"""
    get_inertia(body)

Get inertia tensor of the body.
"""
function get_inertia end

# Body state management interface

"""
    update_inertia_cache!(cache, coords, prop, inst_state)

Update inertia-related cache.
"""
function update_inertia_cache! end

"""
    lazy_update_state!(state, coords, cache, prop, inst_state)

Update only position and velocity (lazy update).
"""
function lazy_update_state! end

"""
    update_state!(state, coords, cache, prop, inst_state)

Update full state including rotation and angular velocity.
"""
function update_state! end

"""
    update_loci_states!(state, coords, cache, prop, inst_state)

Update loci states using cached transformation matrices.
"""
function update_loci_states! end

"""
    stretch_loci!(coords, cache, prop, inst_state)

Update transformation matrices from deformation parameters.
"""
function stretch_loci! end

"""
    update_transformations!(coords, cache, state, prop, inst_state)

Update transformation matrices from coordinates.
"""
function update_transformations! end

# Energy methods

"""
    kinetic_energy(body, inst_state)

Calculate kinetic energy of the body.
"""
function kinetic_energy end

"""
    potential_energy_field(body, field, inst_state)

Calculate potential energy from field.
"""
function potential_energy_field end

"""
    potential_energy_strain(body)

Calculate strain potential energy.
"""
function potential_energy_strain end



# State conversion

"""
    body_state2coords_state(state, coords)

Convert body state to coordinate state.
"""
function body_state2coords_state end

# Strain methods

"""
    strain!(F, prop, coords, state, cache)

Calculate strain from body state.
"""
function strain! end


function get_params!(params,body::AbstractRigidBody,r̄p)
    params .= to_local_coords(body.coords,r̄p)
end

function get_params!(params,body::FlexibleBody,r̄p)
    params .= to_local_coords(body.coords,r̄p)
end

get_num_of_dims(body::BodyType) where {BodyType<:AbstractBody{N,T}} where {N,T} = N

get_numbertype(body::BodyType) where {BodyType<:AbstractBody{N,T}} where {N,T} = T

get_num_of_intrinsic_cstr(body::AbstractBody) = get_num_of_intrinsic_cstr(body.coords)
get_num_of_coords(body::AbstractBody) = get_num_of_coords(body.coords)
get_num_of_dof(body::AbstractBody) = get_num_of_dof(body.coords)
get_num_of_local_dims(body::AbstractBody) = get_num_of_local_dims(body.coords)

has_constant_mass_matrix(body::AbstractBody) = has_constant_mass_matrix(body.coords)

get_gravity(body::AbstractBody) = get_gravity(typeof(body))
get_gravity(::Type{<:AbstractBody{2,T}}) where {T} = SVector(zero(T),        -9.81*one(T))
get_gravity(::Type{<:AbstractBody{3,T}}) where {T} = SVector(zero(T),zero(T),-9.81*one(T))
get_gravity(rbs::AbstractVector{<:AbstractRigidBody}) = get_gravity(eltype(rbs))
