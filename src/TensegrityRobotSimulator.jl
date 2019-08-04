module TensegrityRobotSimulator
include("contact.jl")
using LinearAlgebra
using Parameters
using Quats
using StaticArrays
using Rotations
using .CollisionDetection
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
    name::Symbol
    type::Symbol
    mass::T
    inertia::SArray{Tuple{3,3},T,2,9}
    CoM::SArray{Tuple{3},T,1,3}
    anchorpoints::SArray{Tuple{NP},AnchorPoint{T},1,NP}
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

struct RigidBodyNaturalCoordinates{T,qT,CpT} <: RigidBodyCoordinates{T}
    q::qT
    q̇::qT
    q̈::qT
    M::SArray{Tuple{12,12},T,2,144}
    Cp::CpT
end

abstract type AbstractState end

struct RigidBodyState{T,NP,CoordinatesType} <: AbstractState
    r::MArray{Tuple{3},T,1,3}
    R::MArray{Tuple{3,3},T,2,9}
    ṙ::MArray{Tuple{3},T,1,3}
    ω::MArray{Tuple{3},T,1,3}
    p::SArray{Tuple{NP},MArray{Tuple{3},T,1,3},1,NP}
    F::MArray{Tuple{3},T,1,3}
    τ::MArray{Tuple{3},T,1,3}
    Fanc::SArray{Tuple{NP},MArray{Tuple{3},T,1,3},1,NP}
    τanc::SArray{Tuple{NP},MArray{Tuple{3},T,1,3},1,NP}
    coords::CoordinatesType
end

function point_position(prop,id,r,R)
    p = prop.anchorpoints[id].p
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

function RigidBodyState(rb1prop,r1,R1,ṙ1,ω1,coordstype=Val(:Eberly))
    T = typeof(rb1prop).parameters[1]
    @assert T == eltype(r1) == eltype(R1) == eltype(ṙ1) == eltype(ω1)
    r = MVector{3}(r1)
    R = MMatrix{3,3}(R1)
    ṙ = MVector{3}(ṙ1)
    ω = MVector{3}(ω1)
    np = length(rb1prop.anchorpoints)
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
    pid1::Int64
    pid2::Int64
    𝐬::MArray{Tuple{3},T,1,3}
end

abstract type AbstractJoint end

struct BallJoint <: AbstractJoint
    pid1::Int64
    pid2::Int64
end

numberofconstraint(jo::BallJoint) = 3

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

end # module
