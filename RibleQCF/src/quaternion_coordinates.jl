
function RigidBodyCache(prop::RigidBodyProperty{N,T},qcs::QC,
        pres_idx=Int[],
        cstr_idx=get_cstr_idx(qcs)
    ) where {N,T}
    (;mass,inertia,mass_locus,loci) = prop
    mass_center = mass_locus.position
    num_of_loci = length(loci)
    M = MMatrix{7,7}(Matrix(one(T)*I,7,7))
    M⁻¹ = MMatrix{7,7}(Matrix(one(T)*I,7,7))
    ∂Mq̇∂q = @MMatrix zeros(T,7,7)
    ∂M⁻¹p∂q = @MMatrix zeros(T,7,7)
    Ṁq̇ = @MVector zeros(T,7)
    ∂T∂qᵀ = @MVector zeros(T,7)
    Co = MMatrix{3,7}(zeros(T,3,7))
    for i = 1:3
        Co[i,i] = 1
    end
    Cg = deepcopy(Co)
    Cps = [deepcopy(Co) for i in 1:num_of_loci]
    inertia_cache = InertiaCache(
        sparse(M),
        sparse(M⁻¹),
        sparse(∂Mq̇∂q),
        sparse(∂M⁻¹p∂q),
        Ṁq̇,∂T∂qᵀ,
        true
    )
    RigidBodyCache(
        Co,Cg,Cps,
        inertia_cache,
        nothing
    )
end

function update_inertia_cache!(cache::RigidBodyCache,coords::QC,
        prop::RigidBodyProperty,
        inst_state)
    (; M, M⁻¹, ∂Mq̇∂q, ∂M⁻¹p∂q, ∂T∂qᵀ) = cache.inertia_cache
    # update inertia
    # @show x[4:7],ẋ[4:7]
    qcs = coords
    x = inst_state.q
    ẋ = inst_state.q̇
    M .= build_M(qcs, x)
    M⁻¹ .= build_M⁻¹(qcs, x)
    ∂Mq̇∂q .= build_∂Mẋ∂x(qcs, x, ẋ)
    ∂M⁻¹p∂q .= build_∂M⁻¹y∂x(qcs, x, M * ẋ)
    ∂T∂qᵀ .= build_∂T∂xᵀ(qcs, x, ẋ)
end

function lazy_update_state!(state::RigidBodyState,coords::QC, 
    cache, prop::RigidBodyProperty, inst_state)
    (;origin_frame, mass_locus_state) = state
    r = @view inst_state.q[1:3]
    v = @view inst_state.q̇[1:3]
    R = find_rotation(coords, inst_state.q)
    ω_local = find_local_angular_velocity(coords, inst_state.q, inst_state.q̇)
    ω_world = R * ω_local

    origin_frame.position = mass_locus_state.frame.position = r
    origin_frame.velocity = mass_locus_state.frame.velocity = v
    origin_frame.axes = mass_locus_state.frame.axes = Axes(R)
    origin_frame.local_angular_velocity .= ω_local
    origin_frame.angular_velocity .= ω_world
    mass_locus_state.frame.local_angular_velocity .= ω_local
    mass_locus_state.frame.angular_velocity .= ω_world
end

function update_state!(state::RigidBodyState, coords::QC, 
    cache, prop::RigidBodyProperty, inst_state)
    lazy_update_state!(state, coords, cache, prop, inst_state)
end

function stretch_loci!(coords::QC,
        cache,
        prop::RigidBodyProperty, c
    )
    # to be implemented
end

function update_transformations!(coords::QC, 
    cache, state::RigidBodyState,
    prop::RigidBodyProperty, inst_state)
    (; mass_locus, loci) = prop
    (; Cg, Cps) = cache
    for (Cp, lo) in zip(Cps, loci)
        Cp .= to_position_jacobian(coords, inst_state.q, lo.position)
    end
    Cg .= to_position_jacobian(coords, inst_state.q, mass_locus.position)
end

function update_loci_states!(state::RigidBodyState,coords::QC, 
    cache,prop::RigidBodyProperty, inst_state)
    (; loci, mass_locus) = prop
    (; loci_states, origin_frame, mass_locus_state) = state
    (; Cg, Cps) = cache
    x, ẋ = inst_state.q, inst_state.q̇
    r = @view x[1:3]
    q = @view x[4:7]
    R = Rmat(q)

    # Update mass locus to keep axes/ω in sync with origin frame
    mass_locus_state.frame.position = r .+ R * mass_locus.position
    mass_locus_state.frame.velocity = Cg * ẋ
    mass_locus_state.frame.axes = origin_frame.axes
    mass_locus_state.frame.angular_velocity .= origin_frame.angular_velocity
    mass_locus_state.frame.local_angular_velocity .= origin_frame.local_angular_velocity

    for (i, (locus, locus_state)) in enumerate(zip(loci, loci_states))
        locus_state.frame.position = r .+ R * locus.position
        locus_state.frame.velocity = Cps[i] * ẋ
        locus_state.frame.axes = origin_frame.axes * locus.axes
        locus_state.frame.angular_velocity .= origin_frame.angular_velocity
        locus_state.frame.local_angular_velocity .= origin_frame.local_angular_velocity
    end
end


function cstr_function!(ret, coords::QC, cache::AbstractBodyCache, inst_state::AbstractCoordinatesState)
     cstr_function!(ret, coords, inst_state.q)
end

function cstr_function!(ret, coords::QC, cache::AbstractBodyCache, q::AbstractVector)
    cstr_function!(ret, coords, q)
end

function cstr_jacobian!(ret, coords::QC, cache::AbstractBodyCache, inst_state::AbstractCoordinatesState)
    cstr_jacobian!(ret, coords, cache.jac_cache, inst_state.q)
end

function cstr_velocity_jacobian!(ret, coords::QC, cache::AbstractBodyCache, inst_state::AbstractCoordinatesState)
    cstr_velocity_jacobian!(ret, coords, cache.jac_cache, inst_state.q̇)
end

function cstr_forces_jacobian!(ret, coords::QC, cache::AbstractBodyCache, inst_state::AbstractCoordinatesState)
    cstr_forces_jacobian!(ret, coords, inst_state.q, inst_state.λ)
end

get_cstr_idx(::QC) = collect(1:1)

has_constant_mass_matrix(::QC) = Val{false}()
