"""
Rigid Body Coordinate-Specific Implementations

This module contains coordinate-specific implementations for rigid bodies.
These methods are specific to natural coordinate systems (NCF.NC).
"""

# NCF-specific implementations
@inline _make_jac_cache(::NCF.NC2D2C{T}) where {T} = nothing
@inline _make_jac_cache(::NCF.NC3D3C{T}) where {T} = nothing
@inline _make_jac_cache(::NCF.NC2D4C{T}) where {T} = @MMatrix zeros(T,1,4)
@inline _make_jac_cache(::NCF.NC3D6C{T}) where {T} = @MMatrix zeros(T,1,6)
@inline _make_jac_cache(::NCF.NC2D6C{T}) where {T} = @MMatrix zeros(T,3,6)
@inline _make_jac_cache(::NCF.NC3D12C{T}) where {T} = @MMatrix zeros(T,6,12)
@inline _make_jac_cache(::NCF.NC) = nothing

function RigidBodyCache(prop::RigidBodyProperty{N,T},nmcs::NCF.NC) where {N,T}
    (;mass,inertia,mass_locus,loci) = prop
    mass_center = mass_locus.position
    num_of_loci = length(loci)
    M = NCF.make_M(nmcs,mass,inertia,mass_center)
    M⁻¹ = inv(M)
    ∂Mq̇∂q = zero(M)
    ∂M⁻¹p∂q = zero(M)
    q = @MVector zeros(T,size(M,2))
    Ṁq̇ = @MVector zeros(T,size(M,2))
    ∂T∂qᵀ = @MVector zeros(T,size(M,2))
    c(x) = to_local_coords(nmcs,x)
    C(c) = to_position_jacobian(nmcs,q,c)
    Co = C(c(zero(mass_center)))
    Cg = C(c(mass_center))
    Cps = [typeof(Cg)(C(c(loci[i].position))) for i in 1:num_of_loci]
    inertia_cache = InertiaCache(
        sparse(M),
        sparse(M⁻¹),
        sparse(∂Mq̇∂q),
        sparse(∂M⁻¹p∂q),
        Ṁq̇,
        ∂T∂qᵀ,
        false
    )
    RigidBodyCache(
        Co,Cg,Cps,
        inertia_cache,
        _make_jac_cache(nmcs)
    )
end

function update_inertia_cache!(cache::RigidBodyCache,
        coords::AbstractCoordinates,
        prop::RigidBodyProperty,
        inst_state
    )
        # constant mass matrix, no need to update
end

function lazy_update_state!(state::RigidBodyState,
    coords::AbstractCoordinates, cache,
    prop::RigidBodyProperty, inst_state)
    (;
        origin_frame,
        mass_locus_state,
    ) = state
    (; Co, Cg) = cache
    (;q, q̇) = inst_state
    mul!(origin_frame.position, Co, q)
    mul!(origin_frame.velocity, Co, q̇)
    mul!(mass_locus_state.frame.position, Cg, q)
    mul!(mass_locus_state.frame.velocity, Cg, q̇)
end

function update_state!(state::RigidBodyState,
    coords::AbstractCoordinates, cache,
    prop::RigidBodyProperty, inst_state)
    (;
        origin_frame,
        mass_locus_state,
    ) = state
    (;q, q̇) = inst_state
    lazy_update_state!(state, coords, cache, prop, inst_state)
    axes = Axes(find_rotation(coords, q))
    origin_frame.axes = axes
    origin_frame.angular_velocity = find_angular_velocity(coords, q, q̇)
end

function stretch_loci!(coords::AbstractCoordinates,
        cache,
        prop::RigidBodyProperty, inst_state
    )
    (; loci) = prop
    (; Cps,) = cache
    # to be reimplemented
    # nlocaldim = get_num_of_local_dims(nmcs)
    # for pid in eachindex(loci)
    #     Cps[pid] = NCF.to_position_jacobian(nmcs,q,c[nlocaldim*(pid-1)+1:nlocaldim*pid])
    # end
end

function update_transformations!( coords::AbstractCoordinates, cache, state::RigidBodyState,
    prop::RigidBodyProperty, inst_state)
end

function update_loci_states!(state::RigidBodyState,
    coords::AbstractCoordinates, cache,
    prop::RigidBodyProperty, inst_state)
    (;q, q̇) = inst_state
    (; loci_states) = state
    (; Cps) = cache
    for (i, locus_state) in enumerate(loci_states)
        mul!(locus_state.frame.position, Cps[i], q)
        mul!(locus_state.frame.velocity, Cps[i], q̇)
    end
end

function cstr_function!(ret, coords::NC, cache, inst_state::AbstractCoordinatesState, )
    cstr_function!(ret, coords, inst_state.q)
end

function kinetic_energy(coords::AbstractCoordinates, 
        cache, state::RigidBodyState, prop::RigidBodyProperty,
        inst_state
    )
    (; M) = cache.inertia_cache
    (;q̇) = inst_state
    T = 1 / 2 * transpose(q̇) * M * q̇
end

function cstr_function!(ret, coords::NC, cache::AbstractBodyCache, inst_state::AbstractCoordinatesState)
     cstr_function!(ret, coords, inst_state.q)
end

function cstr_jacobian!(ret, coords::NC, cache::AbstractBodyCache, inst_state::AbstractCoordinatesState)
    cstr_jacobian!(ret, coords, cache.jac_cache, inst_state.q)
end

function cstr_velocity_jacobian!(ret, coords::NC, cache::AbstractBodyCache, inst_state::AbstractCoordinatesState)
    cstr_velocity_jacobian!(ret, coords, cache, inst_state.q̇)
end

function cstr_velocity_jacobian!(ret, coords::NC, cache::AbstractBodyCache, q̇::AbstractVector)
    cstr_velocity_jacobian!(ret, coords, cache.jac_cache, q̇)
end

# Explicit overload for NC without going through NCF module dispatch
function cstr_velocity_jacobian!(ret, coords::NC, jac_cache, q̇::AbstractVector)
    cstr_velocity_jacobian!(ret, coords, q̇)
end

function cstr_forces_jacobian!(ret, coords::NC, cache::AbstractBodyCache, inst_state::AbstractCoordinatesState)
    cstr_forces_jacobian!(ret, coords, cache.hes_cache, inst_state.q)
end
