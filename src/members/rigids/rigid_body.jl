abstract type AbstractRigidBody{N,T} <: AbstractBody{N,T}end
abstract type AbstractRigidBodyProperty{N,T} <: AbstractBodyProperty{N,T} end
abstract type AbstractRigidBodyState{N,T} <: AbstractBodyState{N,T} end

"""
Rigid Body Property Type 
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct RigidBodyProperty{N,T} <: AbstractRigidBodyProperty{N,T}
    "Is able to make contact with?"
    contactable::Bool
    "Is visible?"
    visible::Bool
    "id. Unique in a system"
    id::Int
    "Type or name."
    type::Symbol
    "Mass"
    mass::T
    "Inertia"
    inertia::SMatrix{N,N,T}
    "Centroid in local frame"
    mass_locus::Locus{N,T}
    "Anchor points in local frame"
    loci::Vector{Locus{N,T}}
end

"""
Rigid Body Property Constructor 
$(TYPEDSIGNATURES)
"""
function RigidBodyProperty(
        id::Integer,contactable::Bool,
        mass::T,inertia_input,
        mass_center_position::AbstractVector,
        positions=SVector{size(inertia_input,1),T}[],
        axes_normals=[
            ones(T,size(inertia_input,1)) |> SVector{size(inertia_input,1)}
            for i = 1:length(positions)
        ],
        friction_coefficients=zeros(T,length(positions)),
        restitution_coefficients=zeros(T,length(positions));
        visible = true,
        type = :generic
    ) where T
    mtype = StaticArrays.similar_type(inertia_input)
    VecType = SVector{size(inertia_input,1)}
    mass_locus = Locus(
        VecType(mass_center_position)
    )
    loci = [
        Locus(
            VecType(position),
            VecType(normal),
            friction_coefficient,
            restitution_coefficient
        )
        for (position,normal,friction_coefficient,restitution_coefficient) in zip(
            positions,
            axes_normals,
            friction_coefficients,
            restitution_coefficients
        )
    ]
    # @show loci
    return RigidBodyProperty(
        contactable,
        visible,
        id,type,
        mass,
        mtype(inertia_input),
        mass_locus,
        loci
    )
end

"""
Rigid Body State mutableType.所有坐标在同一个惯性系中表达。
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
mutable struct RigidBodyState{N,M,T} <: AbstractRigidBodyState{N,T}
    "Origin of local frame"
    origin_frame::CartesianFrame{N,M,T}
    "Position of mass center in global frame"
    mass_locus_state::LocusState{N,M,T}
    "Positions of anchor points in global frame"
    loci_states::Vector{LocusState{N,M,T}}
end

"""
Rigid Body State Constructor 
$(TYPEDSIGNATURES)
---
"""
function RigidBodyState(prop::RigidBodyProperty{N,T},
        origin_position_input,
        rotation_input,
        origin_velocity_input,
        angular_velocity_input
    ) where {N,T}
    (;mass_locus,loci) = prop
    num_of_loci = length(loci)
    origin_position = SVector{N}(origin_position_input)
    origin_velocity = SVector{N}(origin_velocity_input)
    if rotation_input isa Number
        rotation = rotation_matrix(rotation_input)
    else
        rotation = Matrix(rotation_input)
    end
    axes = Axes(SMatrix{N,N,T}(rotation))
    M = 2N-3
    ω = SVector{M}(angular_velocity_input)
    origin_frame = CartesianFrame(
        origin_position,
        origin_velocity,
        axes,
        ω
    )
    mass_locus_state = LocusState(
        mass_locus,
        origin_frame,
    )
    loci_states = [
        LocusState(
            lo,
            origin_frame,
        )
        for lo in loci
    ]
    RigidBodyState(
        origin_frame,
        mass_locus_state,
        loci_states,
    )
end

struct RigidBodyCache{CType,cacheType} <: AbstractBodyCache
    Co::CType
    Cg::CType
    Cps::Vector{CType}
    inertia_cache::cacheType
end


function BodyCache(prop::RigidBodyProperty{N,T},
        nmcs::NCF.NC,
        pres_idx=Int[],
        cstr_idx=get_cstr_idx(nmcs)) where {N,T}
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
    c(x) = NCF.to_local_coords(nmcs,x)
    C(c) = NCF.to_transformation(nmcs,q,c)
    Co = C(c(zero(mass_center)))
    Cg = C(c(mass_center))
    Cps = [typeof(Cg)(C(c(loci[i].position))) for i in 1:num_of_loci]
    inertia_cache = InertiaCache(
        M,M⁻¹,
        ∂Mq̇∂q,∂M⁻¹p∂q,
        Ṁq̇,∂T∂qᵀ
    )
    RigidBodyCache(
        Co,Cg,Cps,
        inertia_cache
    )
end

function BodyCache(
        prop::RigidBodyProperty{N,T},
        qcs::QCF.QC,
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
        M,M⁻¹,
        ∂Mq̇∂q,∂M⁻¹p∂q,
        Ṁq̇,∂T∂qᵀ
    )
    RigidBodyCache(
        Co,Cg,Cps,
        inertia_cache
    )
end

"""
通用Rigid Body Type 
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct RigidBody{N,M,T,coordsType,cacheType,meshType} <: AbstractRigidBody{N,T}
    "Rigid Body Property "
    prop::RigidBodyProperty{N,T}
    "Rigid Body State "
    state::RigidBodyState{N,M,T}
    "Coordinates State"
    coords::NonminimalCoordinates{coordsType}
    "Cache"
    cache::cacheType
    "Rigid Body Mesh"
    mesh::meshType
end

function RigidBody(prop,state,coords,mesh=nothing)
    (;nmcs,pres_idx,cstr_idx) = coords
    cache = BodyCache(prop,nmcs,pres_idx,cstr_idx)
    RigidBody(prop,state,coords,cache,mesh)
end

function body_state2coords_state(state::RigidBodyState,coords::NonminimalCoordinates)
    (;origin_frame) = state
    cartesian_frame2coords(coords.nmcs,origin_frame)
end

"""
    update_inertia_cache!(
        cache::RigidBodyCache,
        coords::NonminimalCoordinates{<:NCF.NC},
        prop::RigidBodyProperty,
        q,q̇
    )

update mass matrices/inertia related caches.
"""
function update_inertia_cache!(
        cache::RigidBodyCache,
        coords::NonminimalCoordinates{<:NCF.NC},
        prop::RigidBodyProperty,
        q,q̇
    )
end

function update_inertia_cache!(
        cache::RigidBodyCache,
        coords::NonminimalCoordinates{<:QCF.QC},
        prop::RigidBodyProperty,
        x,ẋ
    )
    (;M,M⁻¹,∂Mq̇∂q,∂M⁻¹p∂q,∂T∂qᵀ) = cache.inertia_cache
    # update inertia
    # @show x[4:7],ẋ[4:7]
    qcs = coords.nmcs
    M .= QCF.build_M(qcs,x)
    M⁻¹ .= QCF.build_M⁻¹(qcs,x)	
    ∂Mq̇∂q .= QCF.build_∂Mẋ∂x(qcs,x,ẋ)
    ∂M⁻¹p∂q .= QCF.build_∂M⁻¹y∂x(qcs,x,M*ẋ)
    ∂T∂qᵀ .= QCF.build_∂T∂xᵀ(qcs,x,ẋ)
end


"""
    lazy_update_state!(state::RigidBodyState,
        coords::NonminimalCoordinates{<:NCF.NC},cache,
        prop::RigidBodyProperty,q,q̇)

update only position and velocity of mass center and local frame's origin
"""
function lazy_update_state!(state::RigidBodyState,
        coords::NonminimalCoordinates{<:NCF.NC},cache,
        prop::RigidBodyProperty,q,q̇)
    (;
        origin_frame,
        mass_locus_state,
    ) = state
    (;Co,Cg) = cache
    origin_frame.position = Co*q
    origin_frame.velocity = Co*q̇
    mass_locus_state.frame.position = Cg*q
    mass_locus_state.frame.velocity = Cg*q̇
end

"""
    update_state!(state::RigidBodyState,
            coords::NonminimalCoordinates{<:NCF.NC},cache,
            prop::RigidBodyProperty,q,q̇)

update position and velocity of mass center and local frame's origin,
and rotation matrix and angular velocity of local frame
"""
function update_state!(state::RigidBodyState,
            coords::NonminimalCoordinates{<:NCF.NC},cache,
            prop::RigidBodyProperty,q,q̇)
    (;
        origin_frame,
        mass_locus_state,
    ) = state
    (;nmcs) = coords
    lazy_update_state!(state,coords,cache,prop,q,q̇)
    axes = Axes(find_rotation(nmcs,q))
    origin_frame.axes = axes
    origin_frame.angular_velocity = find_angular_velocity(nmcs,q,q̇)
end

function lazy_update_state!(state::RigidBodyState,
        coords::NonminimalCoordinates{<:QCF.QC},cache,
        prop::RigidBodyProperty,x,ẋ)
    (;
        origin_frame,
        mass_locus_state,
    ) = state
    origin_frame.position = mass_locus_state.frame.position = x[1:3]
    origin_frame.velocity = mass_locus_state.frame.velocity = ẋ[1:3]
end

function update_state!(state::RigidBodyState,
        coords::NonminimalCoordinates{<:QCF.QC},cache,
        prop::RigidBodyProperty,q,q̇)
    (;
        origin_frame,
        mass_locus_state,
    ) = state
    (;axes) = origin_frame
    (;nmcs) = coords
    lazy_update_state!(state,coords,cache,prop,q,q̇)
    origin_frame.angular_velocity = axes*find_local_angular_velocity(nmcs,q,q̇)
end

"""
    stretch_loci!(
        coords::NonminimalCoordinates{<:NCF.NC},
        cache,
        prop::RigidBodyProperty,c
    )

update transformation matrix of locu from new *c* parameters.
"""
function stretch_loci!(
        coords::NonminimalCoordinates{<:NCF.NC},
        cache,
        prop::RigidBodyProperty,c
    )
    (;loci) = prop
    (;Cps,) = cache
    (;nmcs) = coords
    # to be reimplemented
    # nlocaldim = get_num_of_local_dims(nmcs)
    # for pid in eachindex(loci)
    #     Cps[pid] = NCF.to_transformation(nmcs,q,c[nlocaldim*(pid-1)+1:nlocaldim*pid])
    # end
end

function stretch_loci!(
        coords::NonminimalCoordinates{<:QCF.QC},
        cache,
        prop::RigidBodyProperty,c
    )
    # to be implemented
end

"""
    update_transformations!(
        coords::NonminimalCoordinates{<:NCF.NC},cache,
        state::RigidBodyState,
        prop::RigidBodyProperty,q)

update transformation matrices of loci from new *coordinates*.
"""
function update_transformations!(
        coords::NonminimalCoordinates{<:NCF.NC},cache,
        state::RigidBodyState,
        prop::RigidBodyProperty,q)
end

function update_transformations!(
        coords::NonminimalCoordinates{<:QCF.QC},cache,
        state::RigidBodyState,
        prop::RigidBodyProperty,q)
    (;mass_locus,loci) = prop
    (;nmcs) = coords
    (;Cg,Cps) = cache
    for (Cp,lo) in zip(Cps,loci)
        Cp .= to_transformation(nmcs,q,lo.position)
    end 
    Cg .= to_transformation(nmcs,q,mass_locus.position)
end


"""
    update_loci_states!(state::RigidBodyState,
        coords::NonminimalCoordinates{<:NCF.NC},cache,
        prop::RigidBodyProperty,q,q̇)

update loci's postions and velocities using cached transformation matrices.

"""
function update_loci_states!(state::RigidBodyState,
        coords::NonminimalCoordinates{<:NCF.NC},cache,
        prop::RigidBodyProperty,q,q̇)
    (;loci_states) = state
    (;Cps) = cache
    for (i,locus_state) in enumerate(loci_states)
        locus_state.frame.position = Cps[i]*q
        locus_state.frame.velocity = Cps[i]*q̇
    end
end

function update_loci_states!(state::RigidBodyState,
        coords::NonminimalCoordinates{<:QCF.QC},cache,
        prop::RigidBodyProperty,x,ẋ)
    (;loci) = prop
    (;loci_states,) = state
    (;Cps) = cache
    r = @view x[1:3]
    q = @view x[4:7]
    for (i,(locus,locus_state)) in enumerate(zip(loci,loci_states))
        locus_state.frame.position = r .+ QCF.Rmat(q) * locus.position
        locus_state.frame.velocity = Cps[i]*ẋ
    end
end

function centrifugal_force!(F,state::AbstractBodyState,cache::AbstractBodyCache)
    (;mass_locus_state,loci_states) = state
    (;Cps,Cg,inertia_cache) = cache
    F .+= inertia_cache.∂T∂qᵀ
end

function mass_center_force!(F,state::AbstractBodyState,cache::AbstractBodyCache)
    (;mass_locus_state,loci_states) = state
    (;Cps,Cg,inertia_cache) = cache
    mul!(F,transpose(Cg),mass_locus_state.force,1,1)
end

function concentrated_force!(F,state::AbstractBodyState,cache::AbstractBodyCache)
    (;mass_locus_state,loci_states) = state
    (;Cps,Cg,inertia_cache) = cache
    for (pid,locus_state) in enumerate(loci_states)
        mul!(F,transpose(Cps[pid]),locus_state.force,1,1)
    end
end

get_id(body::AbstractBody) = body.prop.id

"""
Return 约束方程编号。
$(TYPEDSIGNATURES)
"""
function get_cstr_idx(nmcs::NCF.NC)
    num_of_cstr = NCF.get_num_of_cstr(nmcs)
    collect(1:num_of_cstr)
end

get_cstr_idx(::QCF.QC) = collect(1:1)
##

# operations on rigid body
"""
Return Rigid Body translational kinetic energy.
$(TYPEDSIGNATURES)
"""
function kinetic_energy_translation(body::AbstractRigidBody)
    (;mass) = body.prop
    (;mass_center_velocity) = body.state
    T = 1/2*transpose(mass_center_velocity)*mass*mass_center_velocity
end

"""
Return Rigid Body rotational kinetic energy.
$(TYPEDSIGNATURES)
"""
function kinetic_energy_rotation(body::AbstractRigidBody)
    (;inertia) = body.prop
    (;R,ω) = body.state
    Ω = inv(R)*ω
    T = 1/2*transpose(Ω)*inertia*Ω
end

"""
Return Rigid Body gravity potential energy.
$(TYPEDSIGNATURES)
"""
function potential_energy_gravity(body::AbstractRigidBody)
    (;mass) = body.prop
    (;mass_locus_state) = body.state
    gravity_acceleration = get_gravity(body)
    -transpose(mass_locus_state.frame.position)*gravity_acceleration*mass
end

"""
Return Rigid Body Strain potential energy.
$(TYPEDSIGNATURES)
"""
function potential_energy_strain(body::AbstractRigidBody)
    zero(get_numbertype(body))
end


"""
Clear Rigid Body Forces and torques。
$(TYPEDSIGNATURES)
"""
function clear_forces!(body::AbstractRigidBody)
    (;mass_locus_state,loci_states,) = body.state
    mass_locus_state.force .= 0
    mass_locus_state.torque .= 0
    foreach(loci_states) do locus_state
        locus_state.force .= 0
        locus_state.torque .= 0
        reset!(locus_state.contact_state)
    end
end