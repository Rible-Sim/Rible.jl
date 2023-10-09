abstract type AbstractBody{N,T} end
abstract type AbstractBodyProperty{N,T} end
abstract type AbstractBodyState{N,T} end

abstract type AbstractRigidBody{N,T} <: AbstractBody{N,T}end
abstract type AbstractRigidBodyProperty{N,T} <: AbstractBodyProperty{N,T} end
abstract type AbstractRigidBodyState{N,T} <: AbstractBodyState{N,T} end
abstract type ExternalConstraints{T} end

"""
刚体属性类
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
    r̄g::SVector{N,T}
    "Anchor points in local frame"
    r̄ps::Vector{SVector{N,T}}
    "Axes in local frame"
    ās::Vector{SVector{N,T}}
end

"""
刚体属性构造子
$(TYPEDSIGNATURES)
"""
function RigidBodyProperty(
        id::Integer,movable::Bool,
        mass::T,inertia_input,
        r̄g::AbstractVector,
        r̄ps=SVector{size(inertia_input,1),T}[],
        ās=SVector{size(inertia_input,1),T}[];
        constrained = false,
        type = :generic
    ) where T
    nr̄ps = length(r̄ps)
    nās = length(ās)
    if !movable
        constrained = true
    end
    mtype = StaticArrays.similar_type(inertia_input)
    r̄type = SVector{size(inertia_input,1)}
    # ātype = SVector{2size(inertia_input,1)-3}	
    ātype = SVector{size(inertia_input,1)}
    return RigidBodyProperty(
        movable,constrained,
        id,type,
        mass,
        mtype(inertia_input),
        r̄type(r̄g),
        r̄type.(r̄ps),
        ātype.(ās)
    )
end

struct NonminimalCoordinatesCache{fT,MT,JT,VT,ArrayT}
    pres_idx::Vector{Int}
    free_idx::Vector{Int}
    Φ_mask::Vector{Int}
    nΦ::Int
    funcs::fT
    M::MT
    M⁻¹::MT
    ∂Mq̇∂q::JT
    ∂M⁻¹p∂q::JT
    Ṁq̇::VT
    ∂T∂qᵀ::VT
    Co::ArrayT
    Cg::ArrayT
    Cps::Vector{ArrayT}
end

function get_CoordinatesCache(prop::RigidBodyProperty{N,T},
                                 lncs::NCF.LNC,
                                 pres_idx=Int[],
                                 Φ_mask=get_Φ_mask(lncs)) where {N,T}
    (;mass,inertia,r̄g,r̄ps) = prop
    nr̄ps = length(r̄ps)
    free_idx = NCF.get_free_idx(lncs,pres_idx)
    nΦ = length(Φ_mask)
    cf = NCF.CoordinateFunctions(lncs,free_idx,Φ_mask)
    mass_matrix = NCF.make_M(cf,mass,inertia,r̄g)
    M⁻¹ = inv(mass_matrix)
    ∂Mq̇∂q = zero(mass_matrix)
    ∂M⁻¹p∂q = zero(mass_matrix)
    Ṁq̇ = @MVector zeros(T,size(mass_matrix,2))
    ∂T∂qᵀ = @MVector zeros(T,size(mass_matrix,2))
    (;C,c) = cf
    Co = C(c(zero(r̄g)))
    Cg = C(c(r̄g))
    Cps = [typeof(Cg)(C(c(r̄ps[i]))) for i in 1:nr̄ps]
    if prop.movable
        if prop.constrained && pres_idx == Int[]
            @error "Rigid body constrained, but no index specified."
        elseif !prop.constrained && !(pres_idx == Int[])
            @error "Rigid body not constrained. No index should be specified."
        end
    end
    NonminimalCoordinatesCache(pres_idx,free_idx,Φ_mask,nΦ,cf,mass_matrix,M⁻¹,∂Mq̇∂q,∂M⁻¹p∂q,Ṁq̇,∂T∂qᵀ,Co,Cg,Cps)
end

function get_CoordinatesCache(prop::RigidBodyProperty{N,T},
                                qcs::QBF.QC,
                                pres_idx=Int[],
                                Φ_mask=get_Φ_mask(qcs)) where {N,T}
    (;mass,inertia,r̄g,r̄ps) = prop
    nr̄ps = length(r̄ps)
    free_idx = deleteat!(collect(1:7),pres_idx)
    nΦ = length(Φ_mask)
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
    Cps = [deepcopy(Co) for i in 1:nr̄ps]
    if prop.movable
        if prop.constrained && pres_idx == Int[]
            @error "Rigid body constrained, but no index specified."
        elseif !prop.constrained && !(pres_idx == Int[])
            @error "Rigid body not constrained. No index should be specified."
        end
    end
    NonminimalCoordinatesCache(pres_idx,free_idx,Φ_mask,nΦ,cf,mass_matrix,M⁻¹,∂Mq̇∂q,∂M⁻¹p∂q,Ṁq̇,∂T∂qᵀ,Co,Cg,Cps)
end

"""
刚体状态mutable类。所有坐标在同一个惯性系中表达。
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
mutable struct RigidBodyState{N,M,T,cacheType} <: AbstractRigidBodyState{N,T}
    "Origin of local frame"
    ro::MVector{N,T}
    "Rotation Matrix of local frame"
    R::MMatrix{N,N,T}
    "Translational velocity of local frame"
    ṙo::MVector{N,T}
    "Angular velocity of local frame (expressed in the global frame)"
    ω::MVector{M,T}
    "Position of centroid in global frame"
    rg::MVector{N,T}
    "Velocity of centroid in global frame"
    ṙg::MVector{N,T}
    "Positions of anchor points in global frame"
    rps::Vector{MVector{N,T}}
    "Velocities of anchor points in global frame"
    ṙps::Vector{MVector{N,T}}
    "Directions of axes in global frame"
    as::Vector{MVector{N,T}}
    "Angular velocities of axes in global frame"
    ȧs::Vector{MVector{N,T}}
    "Resultant force"
    f::MVector{N,T}
    "Resultant torque"
    τ::MVector{M,T}
    "Forces exerted at anchor points"
    fps::Vector{MVector{N,T}}
    "Torques exerted around axes"
    τps::Vector{MVector{M,T}}
    "Other cache"
    cache::cacheType
end

"""
刚体状态构造子
$(TYPEDSIGNATURES)
---
`pres_idx`为约束坐标的索引。
`Φ_mask`为约束方程的索引。
"""
function RigidBodyState(prop::RigidBodyProperty{N,T},
                        lncs,
                        r_input,rotation_input,ṙ_input,ω_input,
                        pres_idx=Int[],
                        Φ_mask=get_Φ_mask(lncs)) where {N,T}
    (;r̄g,r̄ps,ās) = prop
    nr̄ps = length(r̄ps)
    nās = length(ās)
    ro = MVector{N}(r_input)
    ṙo = MVector{N}(ṙ_input)
    if rotation_input isa Number
        rotation = rotation_matrix(rotation_input)
    else
        rotation = Matrix(rotation_input)
    end
    R = MMatrix{N,N,T}(rotation)
    M = 2N-3
    ω = MVector{M}(ω_input)
    rg = MVector{N}(ro+R*r̄g)
    ṙg = MVector{N}(ṙo+ω×(rg-ro))
    f = @MVector zeros(T,N)
    τ = @MVector zeros(T,M)

    rps = MVector{N,T}[(ro+R*r̄ps[i])      for i in 1:nr̄ps]
    ṙps = MVector{N,T}[(ṙo+ω×(rps[i]-ro)) for i in 1:nr̄ps]
    as  = MVector{N,T}[(R*ās[i])          for i in 1:nās]
    ȧs  = MVector{N,T}[(ω×(as[i]))        for i in 1:nās]
    fps = MVector{N,T}[zeros(T,N) for i in 1:nr̄ps]
    τps = MVector{M,T}[zeros(T,M) for i in 1:nās]

    cache = get_CoordinatesCache(prop,lncs,pres_idx,Φ_mask)

    RigidBodyState(ro,R,ṙo,ω,rg,ṙg,rps,ṙps,as,ȧs,f,τ,fps,τps,cache)
end

"""
通用刚体类
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct RigidBody{N,M,T,cacheType,meshType} <: AbstractRigidBody{N,T}
    "刚体属性"
    prop::RigidBodyProperty{N,T}
    "刚体状态"
    state::RigidBodyState{N,M,T,cacheType}
    "可视化网格"
    mesh::meshType
end

RigidBody(prop,state) = RigidBody(prop,state,nothing)


update_rigid!(rb::AbstractBody,q,q̇) = update_rigid!(rb.state,rb.state.cache,rb.prop,q,q̇)
move_rigid!(rb::AbstractBody,q,q̇)	= move_rigid!(rb.state,rb.state.cache,rb.prop,q,q̇)
stretch_rigid!(rb::AbstractBody,c) = stretch_rigid!(rb.state.cache,rb.prop,c)
update_transformations!(rb::AbstractBody,q) = update_transformations!(rb.state.cache,rb.state,rb.prop,q)

function update_rigid!(state::RigidBodyState,
            cache::NonminimalCoordinatesCache{<:NCF.CoordinateFunctions},
            prop::RigidBodyProperty,q,q̇)
    (;cache,ro,R,ṙo,ω,rg,ṙg) = state
    (;funcs,Co,Cg) = cache
    (;nmcs) = funcs
    mul!(ro, Co, q)
    mul!(ṙo, Co, q̇)
    mul!(rg, Cg, q)
    mul!(ṙg, Cg, q̇)
    R .= NCF.find_R(nmcs,q)
    ω .= NCF.find_ω(nmcs,q,q̇)
end

function update_rigid!(state::RigidBodyState,
            cache::NonminimalCoordinatesCache{<:QBF.CoordinateFunctions},
            prop::RigidBodyProperty,x,ẋ)
    (;cache,ro,R,ṙo,ω,rg,ṙg) = state
    (;M,M⁻¹,∂Mq̇∂q,∂M⁻¹p∂q,∂T∂qᵀ,funcs) = cache
    # update inertia
    # @show x[4:7],ẋ[4:7]
    M .= funcs.build_M(x)
    M⁻¹ .= funcs.build_M⁻¹(x)	
    ∂Mq̇∂q .= funcs.build_∂Mẋ∂x(x,ẋ)
    ∂M⁻¹p∂q .= funcs.build_∂M⁻¹y∂x(x,M*ẋ)
    ∂T∂qᵀ .= funcs.build_∂T∂xᵀ(x,ẋ)
    ro .= rg .= x[1:3]
    ṙo .= ṙg .= ẋ[1:3]
    R .= QBF.find_R(x)
    ω .= R*QBF.find_Ω(x,ẋ)
end

function stretch_rigid!(cache::NonminimalCoordinatesCache{<:NCF.CoordinateFunctions},
                    prop::RigidBodyProperty,c)
    (;r̄ps) = prop
    (;Cps,funcs) = cache
    nlocaldim = get_nlocaldim(cache)
    for pid in eachindex(r̄ps)
        Cps[pid] = funcs.C(c[nlocaldim*(pid-1)+1:nlocaldim*pid])
    end
end

function stretch_rigid!(cache::NonminimalCoordinatesCache{<:QBF.CoordinateFunctions},
                prop::RigidBodyProperty,c)
end

function move_rigid!(state::RigidBodyState,
                    cache::NonminimalCoordinatesCache{<:NCF.CoordinateFunctions},
                    prop::RigidBodyProperty,q,q̇)
    (;rps,ṙps) = state
    (;Cps) = cache
    for (i,(rp,ṙp)) in enumerate(zip(rps,ṙps))
        mul!(rp, Cps[i], q)
        mul!(ṙp, Cps[i], q̇)
    end
end

function move_rigid!(state::RigidBodyState,
                    cache::NonminimalCoordinatesCache{<:QBF.CoordinateFunctions},
                    prop::RigidBodyProperty,q,q̇)
    update_rigid!(state,cache,prop,q,q̇)
    (;r̄ps) = prop
    (;ro,R,ṙo,ω,rps,ṙps) = state
    for (i,(rp,ṙp)) in enumerate(zip(rps,ṙps))
        vp = R * r̄ps[i]
        rp .= ro .+ vp
        ṙp .= ṙo .+ ω × vp
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
    (;r̄g,r̄ps) = prop
    (;R,) = state
    (;Cg,Cps) = cache
    L = QBF.Lmat(q[4:7])
    for (Cp,r̄p) in zip(Cps,r̄ps)
        for i = 1:3 
        	Cp[i,i] = 1
        end
        Cp[1:3,4:7] .= -2R*NCF.skew(r̄p)*L
    end 
    for i = 1:3 
    	Cg[i,i] = 1
    end
    Cg[1:3,4:7] .= -2R*NCF.skew(r̄g)*L
end

function generalize_force!(F,state::AbstractRigidBodyState)
    (;cache,f,fps) = state
    (;Cps,Cg,∂T∂qᵀ) = cache
    for (pid,fp) in enumerate(fps)
        # F .+= transpose(Cps[pid])*fp
        mul!(F,transpose(Cps[pid]),fp,1,1)
    end
    # F .+= transpose(Cg)*f
    mul!(F,transpose(Cg),f,1,1)
    F .+= ∂T∂qᵀ
    F
end

# kinematic joint constraints

function get_rbids(rbs)
    ids = mapreduce((rb)->rb.prop.id,vcat,rbs;init=Int[])
    nb = length(ids)
    ids,nb
end

"""
返回约束方程编号。
$(TYPEDSIGNATURES)
"""
function get_Φ_mask(lncs::NCF.LNC)
    nΦ = NCF.get_nconstraints(lncs)
    collect(1:nΦ)
end

get_Φ_mask(::QBF.QC) = collect(1:1)
##

# operations on rigid body
"""
返回刚体平移动能。
$(TYPEDSIGNATURES)
"""
function kinetic_energy_translation(rb::AbstractRigidBody)
    (;mass) = rb.prop
    (;ṙg) = rb.state
    T = 1/2*transpose(ṙg)*mass*ṙg
end

"""
返回刚体旋转动能。
$(TYPEDSIGNATURES)
"""
function kinetic_energy_rotation(rb::AbstractRigidBody)
    (;inertia) = rb.prop
    (;R,ω) = rb.state
    Ω = inv(R)*ω
    T = 1/2*transpose(Ω)*inertia*Ω
end

"""
返回刚体重力势能。
$(TYPEDSIGNATURES)
"""
function potential_energy_gravity(rb::AbstractRigidBody)
    (;mass) = rb.prop
    (;rg) = rb.state
    gravity_acceleration = get_gravity(rb)
    -transpose(rg)*gravity_acceleration*mass
end

"""
返回刚体应变势能。
$(TYPEDSIGNATURES)
"""
function potential_energy_strain(rb::AbstractRigidBody)
    zero(get_numbertype(rb))
end

