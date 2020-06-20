module TensegrityRobotSimulator
include("contact.jl")
include("NC.jl")
using LinearAlgebra
using Parameters
using Quats
using StaticArrays
using Rotations
using .CollisionDetection
using .NC
export BallJoint
export quat_multiply, point_position, tilde, numberofconstraint


struct AnchorPoint{T}
    p::SArray{Tuple{3},T,1,3}
end
function AnchorPoint(sp)
    p = SVector{3,eltype(sp)}(sp[1], sp[2], sp[3])
    AnchorPoint{eltype(p)}(p)
end

function AnchorPoint(p1,p2,p3)
    p = SVector(p1, p2, p3)
    AnchorPoint{eltype(p)}(p)
end

abstract type AbstractBody end

struct RigidBodyProperty{T,NP}
    movable::Bool
    name::Symbol
    type::Symbol
    mass::T
    inertia::SArray{Tuple{3,3},T,2,9}
    CoM::SArray{Tuple{3},T,1,3}
    aps::SArray{Tuple{NP},AnchorPoint{T},1,NP}
end

struct RigidBodyGeometry
    id::Int64
    # To be implemented
end

abstract type RigidBodyCoordinates{T} end

struct RigidBodyEberlyCoordinates{T} <: RigidBodyCoordinates{T}
    x::MArray{Tuple{3},T,1,3}
    q::MArray{Tuple{4},T,1,4}
    p::MArray{Tuple{3},T,1,3}
    L::MArray{Tuple{3},T,1,3}
    ẋ::MArray{Tuple{3},T,1,3}
    q̇::MArray{Tuple{4},T,1,4}
    ṗ::MArray{Tuple{3},T,1,3}
    L̇::MArray{Tuple{3},T,1,3}
end

struct RigidBodyErlebenCoordinates{T} <: RigidBodyCoordinates{T}
    r::MArray{Tuple{3},T,1,3}
    q::MArray{Tuple{4},T,1,4}
    v::MArray{Tuple{3},T,1,3}
    ω::MArray{Tuple{3},T,1,3}
    ṙ::MArray{Tuple{3},T,1,3}
    q̇::MArray{Tuple{4},T,1,4}
    v̇::MArray{Tuple{3},T,1,3}
    ω̇::MArray{Tuple{3},T,1,3}
end

struct RigidBodyNaturalCoordinates{T} <: RigidBodyCoordinates{T}
    q::MArray{Tuple{12},T,1,12}
    q̇::MArray{Tuple{12},T,1,12}
    q̈::MArray{Tuple{12},T,1,12}
end

abstract type RigidBodyAuxiliaries end
struct NaturalCoordinatesAuxiliaries{T,NP} <: RigidBodyAuxiliaries
    M::SArray{Tuple{12,12},T,2,144}
    CG::SArray{Tuple{3,12},T,2,36}
    Cp::SArray{Tuple{NP},SArray{Tuple{3,12},T,2,36},1,NP}
end
function NaturalCoordinatesAuxiliaries(M,CG,Cp)
    StaticM =SMatrix{12,12}(M)
    StaticCG = SMatrix{3,12}(CG)
    numberofpoints = length(Cp)
    StaticCp = SArray{Tuple{numberofpoints},SArray{Tuple{3,12},eltype(M),2,36},1,numberofpoints}(Cp)
    NaturalCoordinatesAuxiliaries(StaticM,StaticCG,StaticCp)
end
abstract type AbstractState end

struct RigidBodyState{T,NP,CoordinatesType,cacheType} <: AbstractState
    r::MArray{Tuple{3},T,1,3}
    R::MArray{Tuple{3,3},T,2,9}
    ṙ::MArray{Tuple{3},T,1,3}
    ω::MArray{Tuple{3},T,1,3}
    p::SArray{Tuple{NP},MArray{Tuple{3},T,1,3},1,NP} # Anchor Points in global frame
    F::MArray{Tuple{3},T,1,3}
    τ::MArray{Tuple{3},T,1,3}
    Fanc::SArray{Tuple{NP},MArray{Tuple{3},T,1,3},1,NP}
    τanc::SArray{Tuple{NP},MArray{Tuple{3},T,1,3},1,NP}
    coords::CoordinatesType
    cache::cacheType
end

function point_position(prop,id,r,R)
    p = prop.aps[id].p
    r + R*p
end


function RigidBodyEberlyCoordinates(prop,r,ṙ,R,ω)
    @unpack mass,inertia = prop
    x = copy(r)
    quaternion = Quat(R)
    q = MVector{4}(quaternion.w,quaternion.x,quaternion.y,quaternion.z)
    p = ṙ*mass
    Jt = R*inertia*transpose(R)
    L = Jt*ω
    ẋ = copy(ṙ)
    q̇ = MVector{4}(0.5*quat_multiply([0.0,ω[1],ω[2],ω[3]],q))
    ṗ = zero(p)
    L̇ = zero(L)
    RigidBodyEberlyCoordinates(x,q,p,L,ẋ,q̇,ṗ,L̇)
end

function RigidBodyErlebenCoordinates(prop,r,ṙ,R,ω)
    @unpack mass,inertia = prop
    copyr = copy(r)
    quaternion = Quat(R)
    q = MVector{4}(quaternion.w,quaternion.x,quaternion.y,quaternion.z)
    v = copy(ṙ)
    copyω = copy(ω)

    copyṙ = copy(ṙ)
    q̇ = MVector{4}(0.5*quat_multiply([0.0,ω[1],ω[2],ω[3]],q))
    v̇ = zero(v)
    ω̇ = zero(ω)
    RigidBodyErlebenCoordinates(copyr,q,v,copyω,copyṙ,q̇,v̇,ω̇)
end

function coords2state_kinetic!(rb)
    @unpack q,q̇ = rb.state.coords
    ri,rj,u,v = [q[3(k-1)+1:3(k-1)+3] for k = 1:4]
    dot_ri,dot_rj,dot_u,dot_v = [q̇[3(k-1)+1:3(k-1)+3] for k = 1:4]
    @unpack r,R,ṙ,ω,p = rb.state
    r .= ri
    R[:,1] .= rj - ri
    R[:,2] .= u
    R[:,3] .= v
    ṙ .= dot_ri
    Ṙ = zero(R)
    Ṙ[:,1] .= dot_rj-dot_ri
    Ṙ[:,2] .= dot_u
    Ṙ[:,3] .= dot_v
    Ω = Ṙ*transpose(R)
    ω .= [Ω[3,2],Ω[1,3],Ω[2,1]]
    for (id,point) = enumerate(rb.state.p)
        point .= r + R*rb.prop.aps[id].p
    end
end

function RigidBodyNaturalCoordinates(prop,r,ṙ,R,ω)
    @unpack mass,CoM,inertia,aps = prop
    ri = r
    rj = R[:,1] + ri
    u = R[:,2]
    v = R[:,3]
    q = vcat(ri,rj,u,v)
    dot_ri = ṙ + ω × (R*[0.0,0.0,0.0])
    dot_rj = ṙ + ω × (R*[1.0,0.0,0.0])
    dot_u = ω × (R*[0.0,1.0,0.0])
    dot_v = ω × (R*[0.0,0.0,1.0])
    q̇ = vcat(dot_ri,dot_rj,dot_u,dot_v)
    q̈ = zero(q̇)
    coords = RigidBodyNaturalCoordinates(MVector{12}(q),MVector{12}(q̇),MVector{12}(q̈))

    Ji = NC.inertia2Ji(inertia)
    M = NC.form_mass_matrix(mass,CoM,Ji)
    CG = NC.C(NC.c(CoM))
    Cp = [NC.C(NC.c(ap.p)) for ap in aps]

    cache = NaturalCoordinatesAuxiliaries(M,CG,Cp)
    coords,cache
end

function RigidBodyState(rb1prop,r1,R1,ṙ1,ω1,::Val{:Eberly})
    T = typeof(rb1prop).parameters[1]
    @assert T == eltype(r1) == eltype(R1) == eltype(ṙ1) == eltype(ω1)
    r = MVector{3}(r1)
    R = MMatrix{3,3}(R1)
    ṙ = MVector{3}(ṙ1)
    ω = MVector{3}(ω1)
    np = length(rb1prop.aps)
    p = SVector{np,MArray{Tuple{3},T,1,3}}(
        [point_position(rb1prop,id,r,R) for id = 1:np])
    F = @MVector zeros(T,3)
    τ = @MVector zeros(T,3)
    Fanc = SVector{np,MArray{Tuple{3},T,1,3}}(
            [@MVector zeros(T,3) for id = 1:np])
    τanc = SVector{np,MArray{Tuple{3},T,1,3}}(
            [@MVector zeros(T,3) for id = 1:np])
    if typeof(coordstype) <: Val{:Eberly}
        coords = RigidBodyEberlyCoordinates(rb1prop,r,ṙ,R,ω)
    elseif typeof(coordstype) <: Val{:Erleben}
        coords = RigidBodyErlebenCoordinates(rb1prop,r,ṙ,R,ω)
    else
        error("No such $coordstype implemented")
    end
    RigidBodyState(r,R,ṙ,ω,p,F,τ,Fanc,τanc,coords)
end
function RigidBodyState(rb1prop,r1,R1,ṙ1,ω1,::Val{:NC})
    T = typeof(rb1prop).parameters[1]
    @assert T == eltype(r1) == eltype(R1) == eltype(ṙ1) == eltype(ω1)
    r = MVector{3}(r1)
    R = MMatrix{3,3}(R1)
    ṙ = MVector{3}(ṙ1)
    ω = MVector{3}(ω1)
    np = length(rb1prop.aps)
    p = SVector{np,MArray{Tuple{3},T,1,3}}(
        [point_position(rb1prop,id,r,R) for id = 1:np])
    F = @MVector zeros(T,3)
    τ = @MVector zeros(T,3)
    Fanc = SVector{np,MArray{Tuple{3},T,1,3}}(
            [@MVector zeros(T,3) for id = 1:np])
    τanc = SVector{np,MArray{Tuple{3},T,1,3}}(
            [@MVector zeros(T,3) for id = 1:np])

    coords,cache = RigidBodyNaturalCoordinates(rb1prop,r,ṙ,R,ω)

    RigidBodyState(r,R,ṙ,ω,p,F,τ,Fanc,τanc,coords,cache)
end
struct RigidBody{T,NP,CoordinatesType,ObjectType}
    prop::RigidBodyProperty{T,NP}
    state::RigidBodyState{T,CoordinatesType}
    object::ObjectType
end
function RigidBody(prop,state;r=1.0)
    sphere = CollisionDetection.sphere_object(r,state.r)
    RigidBody(prop,state,sphere)
end
function multibodystate(rbs)
    nbody = length(rbs)
    nc = 13
    statevector = Vector{Float64}(undef,nc*nbody)
    statedotvector = Vector{Float64}(undef,nc*nbody)
    for i in eachindex(rbs)
        Eb = rbs[i].state.coords
        state = @view statevector[(i-1)*nc+1:(i-1)*nc+13]
        state[1:3] = Eb.x
        state[4:7] = Eb.q
        state[8:10] = Eb.p
        state[11:13] = Eb.L
        statedot = @view statedotvector[(i-1)*nc+1:(i-1)*nc+13]
        statedot[1:3] = Eb.ẋ
        statedot[4:7] = Eb.q̇
        statedot[8:10] = Eb.ṗ
        statedot[11:13] = Eb.L̇
    end
    statevector,statedotvector
end

abstract type AbstractString end

struct LinearString{T} <: AbstractString
    k::T
    s0::T
    ida::Int64
    idb::Int64
    ra::MArray{Tuple{3},T,1,3}
    rb::MArray{Tuple{3},T,1,3}
    τ::MArray{Tuple{3},T,1,3}
    f::MArray{Tuple{3},T,1,3}
    𝐬::MArray{Tuple{3},T,1,3}
end
function LinearString(k,s0,ida,idb)
    ra = @MVector zeros(3)
    rb = @MVector zeros(3)
    τ = @MVector zeros(3)
    f = @MVector zeros(3)
    𝐬 = @MVector zeros(3)
    LinearString(k,s0,ida,idb,ra,rb,τ,f,𝐬)
end
abstract type AbstractJoint end

struct BallJoint <: AbstractJoint
    pid1::Int64
    pid2::Int64
end

numberofconstraint(jo::BallJoint) = 3

struct RBVector <: AbstractArray{RigidBody, 1}
    rigidbodies::Vector{RigidBody}
    mvbodyindex::Vector{Int}
    movables::Vector{RigidBody}
end
Base.IndexStyle(::Type{RBVector}) = IndexLinear()
Base.size(A::RBVector) = size(A.rigidbodies)
Base.getindex(A::RBVector, i::Int) = getindex(A.rigidbodies, i)
Base.setindex!(A::RBVector, rb::RigidBody, i::Int) = setindex!(A.rigidbodies, rb, i)

struct TensegritySystem{RBType,SType,JType,CNType}
    rigidbodies::RBType
    strings::SType
    joints::JType
    connectivity::CNType
end


struct TensegritySimulator{TGSType,tType,SolType}
    tgsys::TGSType
    tspan::Tuple{tType,tType}
    tgsolution::SolType
end


function quat_multiply(q0,q1)
    w0,x0,y0,z0 = q0
    w1,x1,y1,z1 = q1
    [w0*w1 - x0*x1 - y0*y1 - z0*z1,
     w0*x1 + w1*x0 + y0*z1 - z0*y1,
     w0*y1 + w1*y0 + z0*x1 - x0*z1,
     w0*z1 + w1*z0 + x0*y1 - y0*x1]
end
function tilde(x)
    ret = Matrix{eltype(x)}(undef,3,3)
    ret[1,1] = 0.0
    ret[1,2] = -x[3]
    ret[1,3] = x[2]
    ret[2,1] = x[3]
    ret[2,2] = 0.0
    ret[2,3] = -x[1]
    ret[3,1] = -x[2]
    ret[3,2] = x[1]
    ret[3,3] = 0.0
    ret
end

greet() = print("Hello World!")

function potential_energy(st::LinearString)
    @unpack 𝐬,k,s0 = st
    s_norm = norm(𝐬)
    ds = s_norm - s0
    ds⁺ = ifelse(ds > 0.0, ds, 0.0)
    pe_st = 1/2*k*ds⁺^2
end


function kinetic_energy(rb::RigidBody)
    @unpack mass,inertia = rb.prop
    @unpack r,R,ṙ,ω = rb.state
    translational = 1/2*mass*transpose(ṙ)*ṙ
    local_ω = transpose(R)*ω
    rotational = 1/2*transpose(local_ω) * (inertia * local_ω)
    ke_rb = translational + rotational
end

function compute_string_forces!(tgsys)
    rbs = tgsys.rigidbodies
    sts = tgsys.strings
    for st in sts
        rbida,pida = tgsys.connectivity[st.ida]
        rbidb,pidb = tgsys.connectivity[st.idb]
        rba = rbs[rbida]
        rbb = rbs[rbidb]
        @unpack k,s0,ra,rb,τ,f,𝐬 = st
        ra .= rba.state.p[pida]
        rb .= rbb.state.p[pidb]
        𝐬 .= rb - ra
        s_norm = norm(𝐬)
        τ .= 𝐬./s_norm
        f .= ifelse( s_norm > s0, k*(s_norm-s0).*τ, zero(τ))
        # on body a
        Ca = rba.state.cache.Cp[pida]
        Qa = transpose(Ca)*f
        rba.state.coords.Q .+= Qa
        # on body b
        Cb = rbb.state.cache.Cp[pidb]
        Qb = transpose(Cb)*(-f)
        rbb.state.coords.Q .+= Qb
    end
end

function distribute_q_to_rbs!(rbs,q,q̇)
    for (rbid,rb) = enumerate(rbs)
        is = 12*(rbid-1)
        rb.state.coords.q .= q[is+1:is+12]
        rb.state.coords.q̇ .= q̇[is+1:is+12]
        coords2state_kinetic!(rb)
    end
end

function reset_forces!(rbs)
    for (rbid,rb) in enumerate(rbs)
        rb.state.coords.Q .= 0.0
    end
end
end # module
