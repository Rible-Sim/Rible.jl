abstract type AbstractRigidBody{N,T,C} end
abstract type AbstractRigidBodyProperty{N,T} end
abstract type AbstractRigidBodyState{N,T,C} end
abstract type AbstractTensegrityStructure{N,T} end

struct RigidBodyProperty{N,T,L} <: AbstractRigidBodyProperty{N,T}
    movable::Bool
    constrained::Bool
    name::Symbol
    type::Symbol
    mass::T
    inertia::Union{T,SArray{Tuple{N,N},T,2,L}}
    r̄g::SArray{Tuple{N},T,1,N}
    naps::Int
    aps::Vector{SArray{Tuple{N},T,1,N}}
end

struct RigidBodyCoordinates{N,T}
    q::MArray{Tuple{N},T,1,N}
    q̇::MArray{Tuple{N},T,1,N}
    q̈::MArray{Tuple{N},T,1,N}
    Q::MArray{Tuple{N},T,1,N}
end

struct RigidBodyState{N,T,L,C,cacheType} <: AbstractRigidBodyState{N,T,C}
    ro::MArray{Tuple{N},T,1,N}
    R::MArray{Tuple{N,N},T,2,L}
    ṙo::MArray{Tuple{N},T,1,N}
    ω::Union{T,MArray{Tuple{N},T,1,N}}
    rg::MArray{Tuple{N},T,1,N}
    ṙg::MArray{Tuple{N},T,1,N}
    rps::Vector{MArray{Tuple{N},T,1,N}} # Anchor Points in global frame
    ṙps::Vector{MArray{Tuple{N},T,1,N}} # Anchor Points in global frame
    F::MArray{Tuple{N},T,1,N}
    τ::Union{T,MArray{Tuple{N},T,1,N}}
    Faps::Vector{MArray{Tuple{N},T,1,N}}
    τaps::Union{Vector{T},Vector{MArray{Tuple{N},T,1,N}}}
    coords::C
    cache::cacheType
end


struct RigidBody{N,T,L,C,cacheType} <: AbstractRigidBody{N,T,C}
    prop::RigidBodyProperty{N,T,L}
    state::RigidBodyState{N,T,L,C,cacheType}
end

struct ConstrainedFunctions{ΦT,ΦqT,∂Aᵀλ∂qT}
    Φ::ΦT
    Φq::ΦqT
    ∂Aᵀλ∂q::∂Aᵀλ∂qT
end

struct NaturalCoordinatesCache{ArrayT,MT,fT,cfT}
    Cg::ArrayT
    Co::ArrayT
    Cp::Vector{ArrayT}
    M::MT
    funcs::fT
    nc::Int
    constrained_index::Vector{Int}
    cfuncs::cfT
end

function get_index_Φ_q0(rb)
    q0 = rb.state.coords.q
    index = rb.state.cache.constrained_index
    q0[index]
end

function make_index_Φ(index,q0)
    nΦ = length(index)
    @inline @inbounds inner_Φ(q)   = q[index]-q0[index]
    @inline @inbounds inner_Φ(q,d) = q[index]-d
    inner_Φ
end

function make_index_A(index,nq)
    nΦ = length(index)
    @inline function inner_A(q)
        ret = zeros(eltype(q),nΦ,nq)
        for (iΦ,i) in enumerate(index)
            ret[iΦ,i] = 1
        end
        ret
    end
end

# 2D and 3D
function RigidBodyProperty(i::Integer,movable::Bool,
                            mass::Real,inertia_input,
                            r̄g::AbstractVector,aps;constrained=false)
    name = Symbol("rb"*string(i))
    type = :generic
    naps = length(aps)
    if !movable
        constrained = true
    end
    if typeof(inertia_input)<:Real
        return RigidBodyProperty{2,eltype(r̄g),4}(movable,constrained,name,type,mass,
                                inertia_input,SVector{2}(r̄g),naps,SVector{2}.(aps))
    else
        mtype = similar_type(inertia_input)
        vtype = SVector{size(inertia_input)[1]}
        return RigidBodyProperty(movable,constrained,name,type,mass,
                                mtype(inertia_input),vtype(r̄g),naps,vtype.(aps))
    end
end

function NaturalCoordinatesCache(prop::RigidBodyProperty{N,T,iT},
                                 lncs::NaturalCoordinates.LocalNaturalCoordinates,
                                 q,constrained_index=Vector{Int}()) where {N,T,iT}
    cf = NaturalCoordinates.CoordinateFunctions(lncs)
    @unpack mass,inertia,r̄g,naps,aps = prop
    M = NaturalCoordinates.make_M(cf,mass,inertia,r̄g)
    @unpack C,c = cf
    Cg = C(c(r̄g))
    Co = C(c(zeros(T,N)))
    Cp = [typeof(Cg)(C(c(aps[i]))) for i in 1:naps]
    if prop.movable
        if prop.constrained && constrained_index == Vector{Int}()
            @error "Rigid body constrained, but no index specified."
        elseif !prop.constrained && !(constrained_index == Vector{Int}())
            @error "Rigid body not constrained. No index should be specified."
        end
    end
    nc = length(constrained_index)
    nq = size(Cg)[2]
    q0 = copy(q)
    cfuns = ConstrainedFunctions(
                make_index_Φ(constrained_index,q0),
                make_index_A(constrained_index,nq),
                NaturalCoordinates.∂Φqᵀ∂q_forwarddiff(cf))
    NaturalCoordinatesCache(Cg,Co,Cp,M,cf,nc,constrained_index,cfuns)
end

# 2D
function NaturalCoordinatesCache(prop::RigidBodyProperty,
                                 L::Real,q,constrained_index)
    lncs = NaturalCoordinates.LocalNaturalCoordinates2P(L)
    NaturalCoordinatesCache(prop,lncs,q,constrained_index)
end

function RigidBodyState(prop,ri::AbstractVector{T},
                             θ::T=zero(T),
                             constrained_index=Vector{Int}()) where {T}
    u = [cos(θ),sin(θ)]
    q = MVector{4}(vcat(ri,u))
    q̇ = zero(q)
    q̈ = zero(q)
    Q = zero(q)
    coords = RigidBodyCoordinates(q,q̇,q̈,Q)
    r = MVector{2}(ri)
    ṙ = zero(r)
    R = MMatrix{2,2}(cos(θ), -sin(θ), sin(θ), cos(θ))
    ω = 0.0
    F = MVector(0.0,0.0)
    τ = 0.0
    τaps = MVector(0.0,0.0)
    lncs = NaturalCoordinates.LocalNaturalCoordinates1P1V{T}()
    cache = NaturalCoordinatesCache(prop,lncs,q,constrained_index)
    naps = prop.naps
    rps = [MVector{2}(cache.Cp[i]*q) for i in 1:naps]
    Faps = [zeros(MVector{2}) for i in 1:naps]
    τaps = [zeros(MVector{2}) for i in 1:naps]
    RigidBodyState(r,R,ṙ,ω,rps,F,τ,Faps,τaps,coords,cache)
end

RigidBodyState(prop,ri,rj,constrained_index) = RigidBodyState(prop,ri,rj,zero(ri),zero(rj),constrained_index)

function RigidBodyState(prop,ri::AbstractVector{T},rj::AbstractVector{T},
                             vi::AbstractVector{T},vj::AbstractVector{T},
                             constrained_index=Vector{Int}()) where {T}
    q = MVector{4}(vcat(ri,rj))
    q̇ = MVector{4}(vcat(vi,vj))
    q̈ = zero(q)
    Q = zero(q)
    coords = RigidBodyCoordinates(q,q̇,q̈,Q)
    ro = MVector{2}(ri)
    ṙo = zero(ro)
    u = rj-ri
    θ = atan(u[2],u[1])
    R = MMatrix{2,2}(cos(θ), -sin(θ), sin(θ), cos(θ))
    ω = 0.0
    @unpack r̄g = prop
    rg = ro + R*r̄g
    ṙg = ṙo + ω×rg
    F = MVector(0.0,0.0)
    τ = 0.0
    τaps = MVector(0.0,0.0)
    cache = NaturalCoordinatesCache(prop,norm(u),q,constrained_index)
    naps = prop.naps
    rps = [MVector{2}(cache.Cp[i]*q) for i in 1:naps]
    ṙps = [MVector{2}(cache.Cp[i]*q̇) for i in 1:naps]
    Faps = [zeros(MVector{2}) for i in 1:naps]
    τaps = [zeros(MVector{2}) for i in 1:naps]
    RigidBodyState(ro,R,ṙo,ω,rg,ṙg,rps,ṙps,F,τ,Faps,τaps,coords,cache)
end

function RigidBodyState(prop::RigidBodyProperty{2,T},
                        lncs::NaturalCoordinates.LocalNaturalCoordinates2D6C,
                        r_input,θ_input,ṙ_input,ω_input,
                        q_input,q̇_input=zero(q_input),
                        constrained_index=Vector{Int}()) where {T}
    q = MVector{6}(q_input)
    q̇ = MVector{6}(q̇_input)
    q̈ = zero(q)
    Q = zero(q)
    coords = RigidBodyCoordinates(q,q̇,q̈,Q)
    cache = NaturalCoordinatesCache(prop,lncs,q,constrained_index)

    @unpack r̄g,naps,aps = prop
    @unpack C,c = cache.funcs

    R = MMatrix{2,2,T,4}(rotation_matrix(θ_input))
    rg = MVector{2}(r_input+R*r̄g)
    ω = ω_input ###
    ṙg = MVector{2}(ṙ_input+ω×(rg-r_input))

    Cg = C(c(r̄g))
    Cp = [C(c(ap)) for ap in aps]
    p = [Cp[i]*q for i in 1:naps]
    nap = prop.naps
    rps = [MVector{2}(cache.Cp[i]*q) for i in 1:nap]
    F = @MVector zeros(T,2)
    τ = zero(ω_input) ###
    Faps = [@MVector zeros(T,2) for i in 1:nap]
    τaps = [zero(ω_input) for i in 1:nap]

    RigidBodyState(rg,R,ṙg,ω,rps,F,τ,Faps,τaps,coords,cache)
end

# 3D
# function NaturalCoordinatesCache(prop::RigidBodyProperty{3,T,iT},r::AbstractVector{T}...) where {T,iT}
#     lncs = NaturalCoordinates.LocalNaturalCoordinates3D(r...)
#     NaturalCoordinatesCache(prop,lncs)
# end

function RigidBodyState(prop::RigidBodyProperty{3,T,iT},
                        lncs::NaturalCoordinates.LocalNaturalCoordinates3D,
                        r_input,R_input,ṙ_input,ω_input,
                        q_input,q̇_input=zero(q_input),
                        constrained_index=Vector{Int}()) where {T,iT}
    q = MVector{12}(q_input)
    q̇ = MVector{12}(q̇_input)
    q̈ = zero(q)
    Q = zero(q)
    coords = RigidBodyCoordinates(q,q̇,q̈,Q)
    cache = NaturalCoordinatesCache(prop,lncs,q,constrained_index)
    @unpack r̄g,naps,aps = prop
    @unpack C,c = cache.funcs
    Cg = C(c(r̄g))
    Cp = [C(c(ap)) for ap in aps]
    p = [Cp[i]*q for i in 1:naps]
    nap = prop.naps
    rps = [MVector{3}(cache.Cp[i]*q) for i in 1:nap]
    ṙps = [MVector{3}(cache.Cp[i]*q̇) for i in 1:nap]
    F = @MVector zeros(T,3)
    τ = @MVector zeros(T,3) ###
    Faps = [@MVector zeros(T,3) for i in 1:nap]
    τaps = [@MVector zeros(T,3) for i in 1:nap]
    R = MMatrix{3,3,T,9}(R_input)
    ro = MVector{3}(r_input)
    ω = MVector{3,T}(ω_input) ###
    ṙo = MVector{3}(ṙ_input)
    rg = MVector{3}(ro+R*r̄g)
    ṙg = MVector{3}(ṙo+ω×(rg-ro))
    RigidBodyState(ro,R,ṙo,ω,rg,ṙg,rps,ṙps,F,τ,Faps,τaps,coords,cache)
end

# operations on rigid body
function reset_forces!(rb::RigidBody)
    for f in rb.state.Faps
        f .= 0.0
    end
    rb.state.F .= 0.0
end

function get_initial(rb::AbstractRigidBody{3,T,CT}) where {T,CT}
    @unpack q,q̇ = rb.state.coords
    λ0 = zeros(T,6)
    Array(q),Array(q̇),λ0
end

function kinetic_energy_coords(state::AbstractRigidBodyState,q̇)
    M = state.cache.M
    T = 1/2*transpose(q̇)*M*q̇
end

function kinetic_energy_coords(state::AbstractRigidBodyState)
    M = state.cache.M
    q̇ = state.coords.q̇
    T = 1/2*transpose(q̇)*M*q̇
end

kinetic_energy_coords(rb::AbstractRigidBody) = kinetic_energy_coords(rb.state)
kinetic_energy_coords(rb::AbstractRigidBody,q̇) = kinetic_energy_coords(rb.state,q̇)
