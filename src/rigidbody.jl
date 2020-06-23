abstract type AbstractRigidBody{N,T,CType} end
abstract type AbstractRigidBodyProperty{N,T} end
abstract type AbstractRigidBodyState{N,T,CType} end
abstract type AbstractTensegrityStructure{N,T} end

struct RigidBodyProperty{N,T,iT} <: AbstractRigidBodyProperty{N,T}
    movable::Bool
    constrained::Bool
    name::Symbol
    type::Symbol
    mass::T
    inertia::iT
    CoM::SArray{Tuple{N},T,1,N}
    naps::Int
    aps::Vector{SArray{Tuple{N},T,1,N}}
end

struct RigidBodyState{N,T,L,ωT,CType,cacheType} <: AbstractRigidBodyState{N,T,CType}
    r::MArray{Tuple{N},T,1,N}
    R::MArray{Tuple{N,N},T,2,L}
    ṙ::MArray{Tuple{N},T,1,N}
    ω::ωT
    p::Vector{MArray{Tuple{N},T,1,N}} # Anchor Points in global frame
    F::MArray{Tuple{N},T,1,N}
    τ::ωT
    Faps::Vector{MArray{Tuple{N},T,1,N}}
    τaps::Vector{MArray{Tuple{N},T,1,N}}
    coords::CType
    cache::cacheType
end

struct RigidBodyCoordinates{N,T}
    q::MArray{Tuple{N},T,1,N}
    q̇::MArray{Tuple{N},T,1,N}
    q̈::MArray{Tuple{N},T,1,N}
    Q::MArray{Tuple{N},T,1,N}
end

struct RigidBody{N,T,iT,L,ωT,CType,cacheType} <: AbstractRigidBody{N,T,CType}
    prop::RigidBodyProperty{N,T,iT}
    state::RigidBodyState{N,T,L,ωT,CType,cacheType}
end

struct ConstrainedFunctions{ΦT,ΦqT}
    Φ::ΦT
    Φq::ΦqT
end

struct NaturalCoordinatesCache{ArrayT,MT,fT,cfT}
    CG::ArrayT
    Cp::Vector{ArrayT}
    M::MT
    funcs::fT
    nc::Int
    constrained_index::Vector{Int}
    cfuncs::cfT
end

function make_index_Φ(index,q0)
    nΦ = length(index)
    @inline function inner_Φ(q)
        ret = zeros(eltype(q),nΦ)
        ret .= q[index]-q0[index]
        ret
    end
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
                            CoM::AbstractVector,aps;constrained=false)
    name = Symbol("rb"*string(i))
    type = :generic
    naps = length(aps)
    if !movable
        constrained = true
    end
    if typeof(inertia_input)<:Real
        return RigidBodyProperty(movable,constrained,name,type,mass,
                                inertia_input,CoM,naps,aps)
    else
        @assert size(inertia_input) == (3,3)
        inertia = SMatrix{3,3}(inertia_input)
        return RigidBodyProperty(movable,constrained,name,type,mass,
                                inertia,CoM,naps,aps)
    end
end

function NaturalCoordinatesCache(prop::RigidBodyProperty{N,T,iT},
                                 bps::NaturalCoordinates.BasicPoints,
                                 q,constrained_index=Vector{Int}()) where {N,T,iT}
    cf = NaturalCoordinates.CoordinateFunctions(bps)
    @unpack mass,inertia,CoM,naps,aps = prop
    M = NaturalCoordinates.make_M(cf,mass,inertia,CoM)
    @unpack C,c = cf
    CG = C(c(CoM))
    Cp = [typeof(CG)(C(c(aps[i]))) for i in 1:naps]
    if prop.movable
        if prop.constrained && constrained_index == Vector{Int}()
            @error "Rigid body constrained, but no index specified."
        elseif !prop.constrained && !(constrained_index == Vector{Int}())
            @error "Rigid body not constrained. No index should be specified."
        end
    end
    nc = length(constrained_index)
    nq = size(CG)[2]
    q0 = copy(q)
    cfuns = ConstrainedFunctions(
                make_index_Φ(constrained_index,q0),
                make_index_A(constrained_index,nq))
    NaturalCoordinatesCache(CG,Cp,M,cf,nc,constrained_index,cfuns)
end

# 2D
function NaturalCoordinatesCache(prop::RigidBodyProperty,
                                 L::Real,q,constrained_index)
    bps = NaturalCoordinates.BasicPoints2P(L)
    NaturalCoordinatesCache(prop,bps,q,constrained_index)
end

function RigidBodyState(prop,ri,rj,constrained_index=Vector{Int}())
    q = MVector{4}(vcat(ri,rj))
    q̇ = zero(q)
    q̈ = zero(q)
    Q = zero(q)
    coords = RigidBodyCoordinates(q,q̇,q̈,Q)
    r = MVector{2}(ri)
    ṙ = zero(r)
    u = rj-ri
    θ = atan(u[2],u[1])
    R = MMatrix{2,2}(cos(θ), -sin(θ), sin(θ), cos(θ))
    ω = 0.0
    F = MVector(0.0,0.0)
    τ = 0.0
    τaps = MVector(0.0,0.0)
    cache = NaturalCoordinatesCache(prop,norm(u),q,constrained_index)
    naps = prop.naps
    p = [MVector{2}(cache.Cp[i]*q) for i in 1:naps]
    Faps = [zeros(MVector{2}) for i in 1:naps]
    τaps = [zeros(MVector{2}) for i in 1:naps]
    RigidBodyState(r,R,ṙ,ω,p,F,τ,Faps,τaps,coords,cache)
end

# 3D
# function NaturalCoordinatesCache(prop::RigidBodyProperty{3,T,iT},r::AbstractVector{T}...) where {T,iT}
#     bps = NaturalCoordinates.BasicPoints3D(r...)
#     NaturalCoordinatesCache(prop,bps)
# end

function RigidBodyState(prop::RigidBodyProperty{3,T,iT},bps::NaturalCoordinates.BasicPoints,
                        r_input,R_input,ṙ_input,ω_input,
                        q_input,q̇_input=zero(q_input)) where {T,iT}
    q = MVector{12}(q_input)
    q̇ = MVector{12}(q̇_input)
    q̈ = zero(q)
    Q = zero(q)
    coords = RigidBodyCoordinates(q,q̇,q̈,Q)
    cache = NaturalCoordinatesCache(prop,bps,q,constrained_index)
    @unpack CoM,naps,aps = prop
    @unpack C,c = cache.funcs
    CG = C(c(CoM))
    Cp = [C(c(ap)) for ap in aps]
    p = [Cp[i]*q for i in 1:naps]
    nap = prop.naps
    p = [MVector{3}(cache.Cp[i]*q) for i in 1:nap]
    F = @MVector zeros(T,3)
    τ = @MVector zeros(T,3) ###
    Faps = [@MVector zeros(T,3) for i in 1:nap]
    τaps = [@MVector zeros(T,3) for i in 1:nap]
    R = MMatrix{3,3,T,9}(R_input)
    r = MVector{3}(r_input+R*CoM)
    ω = MVector{3,T}(ω_input) ###
    ṙ = MVector{3}(ṙ_input+ω×(r-r_input))
    RigidBodyState(r,R,ṙ,ω,p,F,τ,Faps,τaps,coords,cache)
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
