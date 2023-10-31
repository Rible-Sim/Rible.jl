abstract type AbstractRigidBody{N,T} <: AbstractBody{N,T}end
abstract type AbstractRigidBodyProperty{N,T} <: AbstractBodyProperty{N,T} end
abstract type AbstractRigidBodyState{N,T} <: AbstractBodyState{N,T} end
abstract type ExternalConstraints{T} end

"""
Rigid Body Property Type 
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct RigidBodyProperty{N,T} <: AbstractRigidBodyProperty{N,T}
    "Is movable?"
    movable::Bool
    "Is constrained?"
    constrained::Bool
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
        id::Integer,movable::Bool,
        mass::T,inertia_input,
        mass_center_position::AbstractVector,
        positions=SVector{size(inertia_input,1),T}[],
        axes_normals=SVector{size(inertia_input,1),T}[],
        friction_coefficients=zeros(T,length(positions)),
        restitution_coefficients=zeros(T,length(positions));
        constrained = false,
        type = :generic
    ) where T
    if !movable
        constrained = true
    end
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
    return RigidBodyProperty(
        movable,constrained,
        id,type,
        mass,
        mtype(inertia_input),
        mass_locus,
        loci
    )
end

struct NonminimalCoordinatesCache{FuncsType,MassMatrixType,JacobianType,HessianType,GenForceType,TransformMatrixType}
    pres_idx::Vector{Int}
    free_idx::Vector{Int}
    cstr_idx::Vector{Int}
    num_of_cstr::Int
    funcs::FuncsType
    cstr_hessians::Vector{HessianType}
    M::MassMatrixType
    M⁻¹::MassMatrixType
    ∂Mq̇∂q::JacobianType
    ∂M⁻¹p∂q::JacobianType
    Ṁq̇::GenForceType
    ∂T∂qᵀ::GenForceType
    Co::TransformMatrixType
    Cg::TransformMatrixType
    Cps::Vector{TransformMatrixType}
end

function get_CoordinatesCache(prop::RigidBodyProperty{N,T},
                                 nmcs::NCF.NC,
                                 pres_idx=Int[],
                                 cstr_idx=get_cstr_idx(nmcs)) where {N,T}
    (;mass,inertia,mass_locus,loci) = prop
    mass_center = mass_locus.position
    num_of_loci = length(loci)
    free_idx = NCF.get_free_idx(nmcs,pres_idx)
    num_of_cstr = length(cstr_idx)
    cf = NCF.CoordinateFunctions(nmcs,free_idx,cstr_idx)
    cstr_hessians = make_cstr_hessians(nmcs)
    mass_matrix = NCF.make_M(cf,mass,inertia,mass_center)
    M⁻¹ = inv(mass_matrix)
    ∂Mq̇∂q = zero(mass_matrix)
    ∂M⁻¹p∂q = zero(mass_matrix)
    Ṁq̇ = @MVector zeros(T,size(mass_matrix,2))
    ∂T∂qᵀ = @MVector zeros(T,size(mass_matrix,2))
    # (;C,c) = cf
    c(x) = NCF.to_local_coords(nmcs,x)
    C(c) = NCF.to_transformation(nmcs,c)
    Co = C(c(zero(mass_center)))
    Cg = C(c(mass_center))
    Cps = [typeof(Cg)(C(c(loci[i].position))) for i in 1:num_of_loci]
    if prop.movable
        if prop.constrained && pres_idx == Int[]
            @error "Rigid body constrained, but no index specified."
        elseif !prop.constrained && !(pres_idx == Int[])
            @error "Rigid body not constrained. No index should be specified."
        end
    end
    NonminimalCoordinatesCache(
        pres_idx,free_idx,
        cstr_idx,num_of_cstr,
        cf,cstr_hessians,
        mass_matrix,M⁻¹,
        ∂Mq̇∂q,∂M⁻¹p∂q,Ṁq̇,
        ∂T∂qᵀ,
        Co,Cg,Cps
    )
end

function get_CoordinatesCache(prop::RigidBodyProperty{N,T},
                                qcs::QBF.QC,
                                pres_idx=Int[],
                                cstr_idx=get_cstr_idx(qcs)) where {N,T}
    (;mass,inertia,mass_locus,loci) = prop
    mass_center = mass_locus.position
    num_of_loci = length(loci)
    free_idx = deleteat!(collect(1:7),pres_idx)
    num_of_cstr = length(cstr_idx)
    cf = QBF.CoordinateFunctions(qcs)
    mass_matrix = MMatrix{7,7}(Matrix(one(T)*I,7,7))
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
    if prop.movable
        if prop.constrained && pres_idx == Int[]
            @error "Rigid body constrained, but no index specified."
        elseif !prop.constrained && !(pres_idx == Int[])
            @error "Rigid body not constrained. No index should be specified."
        end
    end
    NonminimalCoordinatesCache(pres_idx,free_idx,cstr_idx,num_of_cstr,cf,mass_matrix,M⁻¹,∂Mq̇∂q,∂M⁻¹p∂q,Ṁq̇,∂T∂qᵀ,Co,Cg,Cps)
end

"""
Rigid Body State mutableType.所有坐标在同一个惯性系中表达。
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
mutable struct RigidBodyState{N,M,T,cacheType} <: AbstractRigidBodyState{N,T}
    "Origin of local frame"
    origin_position::MVector{N,T}
    "Rotation Matrix of local frame"
    R::MMatrix{N,N,T}
    "Translational velocity of local frame"
    origin_velocity::MVector{N,T}
    "Angular velocity of local frame (expressed in the global frame)"
    ω::MVector{M,T}
    "Position of mass center in global frame"
    mass_locus_state::LocusState{N,M,T}
    "Positions of anchor points in global frame"
    loci_states::Vector{LocusState{N,M,T}}
    "Other cache"
    cache::cacheType
end

"""
Rigid Body State Constructor 
$(TYPEDSIGNATURES)
---
`pres_idx`为约束坐标的索引。
`cstr_idx`为约束方程的索引。
"""
function RigidBodyState(prop::RigidBodyProperty{N,T},
                        nmcs,
                        origin_position_input,
                        rotation_input,
                        origin_velocity_input,
                        angular_velocity_input,
                        pres_idx=Int[],
                        cstr_idx=get_cstr_idx(nmcs)) where {N,T}
    (;mass_locus,loci) = prop
    num_of_loci = length(loci)
    origin_position = MVector{N}(origin_position_input)
    origin_velocity = MVector{N}(origin_velocity_input)
    if rotation_input isa Number
        rotation = rotation_matrix(rotation_input)
    else
        rotation = Matrix(rotation_input)
    end
    R = MMatrix{N,N,T}(rotation)
    M = 2N-3
    ω = MVector{M}(angular_velocity_input)
    mass_locus_state = LocusState(
        mass_locus,
        origin_position,R,
        origin_velocity,ω
    )
    loci_states = [
        LocusState(
            lo,
            origin_position,R,
            origin_velocity,ω
        )
        for lo in loci
    ]
    cache = get_CoordinatesCache(prop,nmcs,pres_idx,cstr_idx)
    RigidBodyState(
        origin_position,R,
        origin_velocity,ω,
        mass_locus_state,
        loci_states,
        cache,
    )
end

"""
通用Rigid Body Type 
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct RigidBody{N,M,T,cacheType,meshType} <: AbstractRigidBody{N,T}
    "Rigid Body Property "
    prop::RigidBodyProperty{N,T}
    "Rigid Body State "
    state::RigidBodyState{N,M,T,cacheType}
    "Rigid Body Mesh"
    mesh::meshType
end

RigidBody(prop,state) = RigidBody(prop,state,nothing)


function body2coords(body::RigidBody)
    (;origin_position,R,origin_velocity,ω,cache) = body.state
    cartesian_frame2coords(cache.funcs.nmcs,origin_position,R,origin_velocity,ω)
end

function update_body!(state::RigidBodyState,
            cache::NonminimalCoordinatesCache{<:NCF.CoordinateFunctions},
            prop::RigidBodyProperty,q,q̇)
    (;cache,
    origin_position,R,
    origin_velocity,ω,
    mass_locus_state,
    ) = state
    (;funcs,Co,Cg) = cache
    (;nmcs) = funcs
    mul!(origin_position, Co, q)
    mul!(origin_velocity, Co, q̇)
    mul!(mass_locus_state.position, Cg, q)
    mul!(mass_locus_state.velocity, Cg, q̇)
    R .= NCF.find_rotation(nmcs,q)
    ω .= NCF.find_angular_velocity(nmcs,q,q̇)
end

function update_body!(state::RigidBodyState,
            cache::NonminimalCoordinatesCache{<:QBF.CoordinateFunctions},
            prop::RigidBodyProperty,x,ẋ)
    (;cache,
    origin_position,R,
    origin_velocity,ω,
    mass_locus_state,
    ) = state
    (;M,M⁻¹,∂Mq̇∂q,∂M⁻¹p∂q,∂T∂qᵀ,funcs) = cache
    # update inertia
    # @show x[4:7],ẋ[4:7]
    M .= funcs.build_M(x)
    M⁻¹ .= funcs.build_M⁻¹(x)	
    ∂Mq̇∂q .= funcs.build_∂Mẋ∂x(x,ẋ)
    ∂M⁻¹p∂q .= funcs.build_∂M⁻¹y∂x(x,M*ẋ)
    ∂T∂qᵀ .= funcs.build_∂T∂xᵀ(x,ẋ)
    origin_position .= mass_locus_state.position .= x[1:3]
    origin_velocity .= mass_locus_state.velocity .= ẋ[1:3]
    R .= QBF.find_rotation(x)
    ω .= R*QBF.find_local_angular_velocity(x,ẋ)
end

function stretch_body!(cache::NonminimalCoordinatesCache{<:NCF.CoordinateFunctions},
                    prop::RigidBodyProperty,c)
    (;loci) = prop
    (;Cps,funcs) = cache
    nlocaldim = get_num_of_local_dims(cache)
    for pid in eachindex(loci)
        Cps[pid] = NCF.to_transformation(funcs.nmcs,c[nlocaldim*(pid-1)+1:nlocaldim*pid])
    end
end

function stretch_body!(cache::NonminimalCoordinatesCache{<:QBF.CoordinateFunctions},
                prop::RigidBodyProperty,c)
end

function move_body!(state::RigidBodyState,
                    cache::NonminimalCoordinatesCache{<:NCF.CoordinateFunctions},
                    prop::RigidBodyProperty,q,q̇)
    (;loci_states) = state
    (;Cps) = cache
    for (i,locus_state) in enumerate(loci_states)
        mul!(locus_state.position, Cps[i], q)
        mul!(locus_state.velocity, Cps[i], q̇)
    end
end

function move_body!(state::RigidBodyState,
                    cache::NonminimalCoordinatesCache{<:QBF.CoordinateFunctions},
                    prop::RigidBodyProperty,q,q̇)
    update_body!(state,cache,prop,q,q̇)
    (;loci) = prop
    (;loci_states,
    origin_position,R,
    origin_velocity,ω,
    ) = state
    for (locus,locus_state) in zip(loci,loci_states)
        relative_velocity = R * locus.position
        locus_state.position .= origin_position .+ relative_velocity
        locus_state.velocity .= origin_velocity .+ ω × relative_velocity
    end
end

function update_transformations!(
        cache::NonminimalCoordinatesCache{<:NCF.CoordinateFunctions},
        state::RigidBodyState,
        prop::RigidBodyProperty,q)
end

function update_transformations!(
        cache::NonminimalCoordinatesCache{<:QBF.CoordinateFunctions},
        state::RigidBodyState,
        prop::RigidBodyProperty,q)
    (;mass_locus,loci) = prop
    (;R,) = state
    (;Cg,Cps) = cache
    L = QBF.Lmat(q[4:7])
    for (Cp,lo) in zip(Cps,loci)
        for i = 1:3 
        	Cp[i,i] = 1
        end
        Cp[1:3,4:7] .= -2R*skew(lo.position)*L
    end 
    for i = 1:3 
    	Cg[i,i] = 1
    end
    Cg[1:3,4:7] .= -2R*skew(mass_locus.position)*L
end

function generalize_force!(F,state::AbstractRigidBodyState)
    (;cache,mass_locus_state,loci_states) = state
    (;Cps,Cg,∂T∂qᵀ) = cache
    for (pid,locus_state) in enumerate(loci_states)
        mul!(F,transpose(Cps[pid]),locus_state.force,1,1)
    end
    mul!(F,transpose(Cg),mass_locus_state.force,1,1)
    F .+= ∂T∂qᵀ
    F
end

# kinematic joint cstr

function get_bodies_ids(bodies)
    ids = mapreduce((body)->body.prop.id,vcat,bodies;init=Int[])
    nb = length(ids)
    ids,nb
end

"""
Return 约束方程编号。
$(TYPEDSIGNATURES)
"""
function get_cstr_idx(nmcs::NCF.NC)
    num_of_cstr = NCF.get_num_of_cstr(nmcs)
    collect(1:num_of_cstr)
end

get_cstr_idx(::QBF.QC) = collect(1:1)
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
    -transpose(mass_locus_state.position)*gravity_acceleration*mass
end

"""
Return Rigid Body Strain potential energy.
$(TYPEDSIGNATURES)
"""
function potential_energy_strain(body::AbstractRigidBody)
    zero(get_numbertype(body))
end


function clear_forces!(body::AbstractRigidBody)
    (;mass_locus_state,loci_states,) = body.state
    mass_locus_state.force .= 0
    mass_locus_state.torque .= 0
    foreach(loci_states) do locus_state
        locus_state.force .= 0
        locus_state.torque .= 0
    end
end