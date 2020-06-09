module Robot2D

using LinearAlgebra
using Parameters
using StaticArrays

greet() = print("Hello World!")

struct VString{T}
    k::T
    l0::MArray{Tuple{2},T,1,2}
    lengths::MArray{Tuple{2},T,1,2}
    tensions::MArray{Tuple{2},T,1,2}
end

struct DString{T}
    k::T
    original_restlengths::SArray{Tuple{4},T,1,4}
    restlengths::MArray{Tuple{4},T,1,4}
    lengths::MArray{Tuple{4},T,1,4}
    tensions::MArray{Tuple{4},T,1,4}
end

struct NString{T,N}
    k::SArray{Tuple{N},T,1,N}
    original_restlengths::SArray{Tuple{N},T,1,N}
    restlengths::MArray{Tuple{N},T,1,N}
    lengths::MArray{Tuple{N},T,1,N}
    tensions::MArray{Tuple{N},T,1,N}
end

mutable struct SStringState{T}
    restlength::T
    length::T
    lengthdot::T
    tension::T
end

function SStringState(restlength,length,tension)
    SStringState(restlength,length,zero(length),tension)
end

struct SString{T}
    k::T
    c::T
    original_restlength::T
    state::SStringState{T}
end

function SString(k,origin_restlength,state)
    SString(k,zero(k),origin_restlength,state)
end

struct Actuator{T}
    strings::Vector{SString{T}}
end

struct RigidBody2DProperty{T}
    movable::Bool
    name::Symbol
    type::Symbol
    mass::T
    inertia::T
    CoM::SArray{Tuple{2},T,1,2}
    number_aps::Int
    anchorpoints::Vector{SArray{Tuple{2},T,1,2}}
end

struct ID
    rbid::Int
    apid::Int
end

struct Connectivity{BPConnectType,StringConnectType}
    body2q::BPConnectType
    string2bp::StringConnectType
end

struct RigidBody2DState{T,CoordinatesType,AuxiliariesType}
    r::MArray{Tuple{2},T,1,2}
    θ::T
    ṙ::MArray{Tuple{2},T,1,2}
    ω::T
    p::Vector{MArray{Tuple{2},T,1,2}} # Anchor Points in global frame
    F::MArray{Tuple{2},T,1,2}
    τ::T
    Fanc::Vector{MArray{Tuple{2},T,1,2}}
    coords::CoordinatesType
    auxs::AuxiliariesType
end

struct RigidBody2DNaturalCoordinates{T}
    q::MArray{Tuple{4},T,1,4}
    q̇::MArray{Tuple{4},T,1,4}
    q̈::MArray{Tuple{4},T,1,4}
end

function RigidBody2DState(prop,ri,rj)
    q = MVector{4}(vcat(ri,rj))
    q̇ = zero(q)
    q̈ = zero(q)
    coords = RigidBody2DNaturalCoordinates(q,q̇,q̈)
    r = MVector{2}(ri)
    ṙ = zero(r)
    rirj = rj-ri
    θ = atan(rirj[2],rirj[1])
    ω = 0.0
    F = MVector(0.0,0.0)
    τ = 0.0
    f1 = MVector(0.0,0.0)
    f2 = MVector(0.0,0.0)
    τanc = MVector(0.0,0.0)
    aux = NCaux(prop,ri,rj)
    nap = prop.number_aps
    p = [MVector{2}(aux.Cp[i]*q) for i in 1:nap]
    Fanc = [zeros(MVector{2}) for i in 1:nap]
    RigidBody2DState(r,θ,ṙ,ω,p,F,τ,Fanc,coords,aux)
end

struct RigidBody2D{T,CoordinatesType,AuxiliariesType}
    prop::RigidBody2DProperty{T}
    state::RigidBody2DState{T,CoordinatesType,AuxiliariesType}
end

struct NaturalCoordinatesAuxiliaries2D{T,cT,ΦT,ΦqT}
    M::SArray{Tuple{4,4},T,2,16}
    CG::SArray{Tuple{2,4},T,2,8}
    Cp::Vector{SArray{Tuple{2,4},T,2,8}}
    Q::MArray{Tuple{4},T,1,4}
    Lij::T
    c::cT
    Φ::ΦT
    Φq::ΦqT
end

struct Structure2D{BodyType,StringType,ActuatorType,ConnectType}
    nbody::Int
    nmovablebody::Int
    mvbodyindex::Vector{Int}
    nfixbody::Int
    fixbodyindex::Vector{Int}
    nstring::Int
    rigidbodies::Vector{BodyType}
    strings::Vector{StringType}
    actuators::Vector{ActuatorType}
    connectivity::ConnectType
end

function Structure2D(rbs,ss,acs,cnt)
    nbody = length(rbs)
    mvbodyindex = [i for i in eachindex(rbs) if rbs[i].prop.movable]
    nmvbody = length(mvbodyindex)
    fixbodyindex = [i for i in eachindex(rbs) if !rbs[i].prop.movable]
    nfixbody = length(fixbodyindex)
    nstring = length(ss)
    Structure2D(nbody,nmvbody,mvbodyindex,
                      nfixbody,fixbodyindex,
                      nstring,
                rbs,ss,acs,cnt)
end

function inertia2Z(inertia,Lij)
    z = 1/Lij^2*inertia
end

function form_mass_matrix(m,CoM,Lij,z)
    M = zeros(4,4)
    a = CoM./Lij
    M[1,1] = m - 2m*a[1] + z
    M[1,2] = 0.0
    M[1,3] = m*a[1] - z
    M[1,4] = -m*a[2]
    M[2,2] = m - 2m*a[1] + z
    M[2,3] = m*a[2]
    M[2,4] = m*a[1] - z
    M[3,3] = z
    M[3,4] = 0.0
    M[4,4] = z
    SymM = Matrix(Symmetric(M))
end

function C(c)
    ret = Matrix{eltype(c)}(undef,2,4)
    ret[1,1] = 1-c[1]
    ret[1,2] = c[2]
    ret[1,3] = c[1]
    ret[1,4] = -c[2]
    ret[2,1] = -c[2]
    ret[2,2] = 1-c[1]
    ret[2,3] = c[2]
    ret[2,4] = c[1]
    ret
end

function NCaux(prop,ri,rj)
    @unpack mass,CoM,inertia,anchorpoints = prop
    Lij = norm(ri-rj)
    z = inertia2Z(inertia,Lij)
    M = form_mass_matrix(mass,CoM,Lij,z)
    X̄⁻¹ = 1/Lij*Matrix(1.0I,2,2)
    function c(r̄)
        ret = X̄⁻¹*r̄
    end
    CG = C(c(CoM))
    Q = zeros(4)
    function Φ(q)
        xi,yi,xj,yj = q
        (xj-xi)^2 + (yj-yi)^2 - Lij^2
    end
    function Φq(q)
        xi,yi,xj,yj = q
        ret = similar(q)
        ret[1] =  2(xi-xj)
        ret[2] =  2(yi-yj)
        ret[3] =  2(xj-xi)
        ret[4] =  2(yj-yi)
        ret
    end
    nap = prop.number_aps
    Cp = [SMatrix{2,4}(C(c(anchorpoints[i])))
            for i in 1:nap]
    aux = NaturalCoordinatesAuxiliaries2D(
    SMatrix{4,4}(M),
    SMatrix{2,4}(CG),
    Cp,
    MVector{4}(Q),
    Lij,c,Φ,Φq
    )
end

function lengthdir(v)
    l = norm(v)
    τ = v/l
    l,τ
end

function reset_forces!(rbs)
    for (rbid,rb) in enumerate(rbs)
        for f in rb.state.Fanc
            f .= 0.0
        end
    end
end

function update_forces!(st2d)
    rbs = st2d.rigidbodies
    ss = st2d.strings
    cnt = st2d.connectivity
    for (istr,sstring) in enumerate(ss)
        @unpack k,c = sstring
        sstate = sstring.state
        a,b = cnt.string2bp[istr]
        state1 = rbs[a.rbid].state
        p1 = state1.p[a.apid]
        ṗ1 = state1.auxs.Cp[a.apid]*state1.coords.q̇
        f1 = state1.Fanc[a.apid]
        state2 = rbs[b.rbid].state
        p2 = state2.p[b.apid]
        ṗ2 = state2.auxs.Cp[b.apid]*state2.coords.q̇
        f2 = state2.Fanc[b.apid]
        Δr = p2 - p1
        Δṙ = ṗ2 - ṗ1
        l,τ = lengthdir(p2-p1)
        sstate.length = l
        l̇ = 1/l*(Δr[1]*Δṙ[1] + Δr[2]*Δṙ[2])
        f_raw = k*(l - sstate.restlength) +
                c*l̇
        sstate.tension = ifelse(f_raw > 0.0, f_raw, 0.0)
        f = τ*sstate.tension
        f1 .+= f
        f2 .+= -f
    end
end

function q2rbstate!(st2d,globalq,globalq̇)
    rbs = st2d.rigidbodies
    cnt = st2d.connectivity
    for (rbid,rb) in enumerate(rbs)
        pindex = cnt.body2q[rbid]
        @unpack q, q̇ = rb.state.coords
        q .= globalq[pindex]
        q̇ .= globalq̇[pindex]
        @unpack auxs,p = rb.state
        for (i,ap) in enumerate(p)
            ap .= auxs.Cp[i]*q
        end
    end
end

function genforces!(rbs)
    for (rbid,rb) in enumerate(rbs)
        @unpack state = rb
        @unpack Fanc = state
        @unpack Q,Cp = state.auxs
        Q .= 0.0
        for (pid,f) in enumerate(Fanc)
            Q .+= transpose(Cp[pid])*f
        end
    end
end

function kineticenergy(rbs)
    ke = 0.0
    for (rbid,rb) in enumerate(rbs)
        @unpack q̇ = rb.state.coords
        @unpack M = rb.state.auxs
        ke += 1/2*transpose(q̇)*M*q̇
    end
    ke
end

function potentialenergy(ss)
    pe = 0.0
    for (ssid,sstring) in enumerate(ss)
        @unpack k,original_restlength = sstring
        sstate = sstring.state
        Δ1 = sstate.length-original_restlength
        if Δ1 > 0.0
            pe += 1/2*k*Δ1^2
        end
    end
    pe
end

function energy(q,q̇,st2d)
    q2rbstate!(st2d,q,q̇)
    ss = st2d.strings
    rbs = st2d.rigidbodies
    update_forces!(st2d)
    ke = kineticenergy(rbs)
    pe = potentialenergy(ss)
    ke + pe
end

end # module
