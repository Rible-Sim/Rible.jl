module Robot2D

using LinearAlgebra
using Parameters
using StaticArrays
using Makie

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

struct AnchorPoint2D{T}
    p::SArray{Tuple{2},T,1,2}
end
function AnchorPoint2D(sp)
    p = SVector{2,eltype(sp)}(sp[1], sp[2])
    AnchorPoint2D{eltype(p)}(p)
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

struct RigidBody2DState{T,NP,CoordinatesType,AuxiliariesType}
    r::MArray{Tuple{2},T,1,2}
    θ::T
    ṙ::MArray{Tuple{2},T,1,2}
    ω::T
    p::SArray{Tuple{NP},MArray{Tuple{2},T,1,2},1,NP} # Anchor Points in global frame
    F::MArray{Tuple{2},T,1,2}
    τ::T
    Fanc::SArray{Tuple{NP},MArray{Tuple{2},T,1,2},1,NP}
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
    p = SVector{nap}(
        [MVector{2}(aux.Cp[i]*q) for i in 1:nap])
    Fanc = zero(p)
    RigidBody2DState(r,θ,ṙ,ω,p,F,τ,Fanc,coords,aux)
end
struct RigidBody2D{T,NP,CoordinatesType,AuxiliariesType}
    prop::RigidBody2DProperty{T}
    state::RigidBody2DState{T,NP,CoordinatesType,AuxiliariesType}
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

struct Structure2D{BodyType,StringType,ConnectType}
    rigidbodies::Vector{BodyType}
    strings::Vector{StringType}
    connectivity::ConnectType
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
    vss = st2d.strings
    cnt = st2d.connectivity
    for (i,vs) in enumerate(vss)
        @unpack k,restlengths,lengths,tensions = vs
        rb1,rb2,rb3 = rbs[2i-1:2i+1]
        P1 = rb1.state.p[1]
        P2 = rb1.state.p[2]
        O1 = rb2.state.p[1]
        O2 = rb2.state.p[2]
        P3 = rb3.state.p[1]
        P4 = rb3.state.p[2]
        P3P1 = P1-P3
        P3O1 = O1-P3
        P4P2 = P2-P4
        P4O1 = O1-P4
        lP3P1,τP3P1 = lengthdir(P3P1)
        lP3O1,τP3O1 = lengthdir(P3O1)
        lP4P2,τP4P2 = lengthdir(P4P2)
        lP4O1,τP4O1 = lengthdir(P4O1)
        lengths[1] = lP3P1
        lengths[2] = lP3O1
        lengths[3] = lP4O1
        lengths[4] = lP4P2
        for i = 1:4
            Δi = lengths[i] - restlengths[i]
            tensions[i] = ifelse(Δi > 0.0, Δi*k, 0.0)
        end
        fP1 = rb1.state.Fanc[1]
        fP2 = rb1.state.Fanc[2]
        fO1 = rb2.state.Fanc[1]
        fP3 = rb3.state.Fanc[1]
        fP4 = rb3.state.Fanc[2]
        fP1 .+= -τP3P1*tensions[1]
        fP3 .+=  τP3P1*tensions[1] + τP3O1*tensions[2]
        fP2 .+= -τP4P2*tensions[4]
        fP4 .+=  τP4P2*tensions[4] + τP4O1*tensions[3]
        fO1 .+= -τP3O1*tensions[2] - τP4O1*tensions[3]
    end
end
function q2rbstate!(st2d,globalq,globalq̇)
    rbs = st2d.rigidbodies
    cnt = st2d.connectivity
    for (rbid,rb) in enumerate(rbs)
        pindex = cnt[rbid]
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

function potentialenergy(vss)
    pe = 0.0
    for (vsid,vs) in enumerate(vss)
        @unpack k,l0,lengths = vs
        Δ1 = lengths[1]-l0[1]
        Δ2 = lengths[2]-l0[2]
        if Δ1 > 0.0
            pe += 1/2*k*Δ1^2
        end
        if Δ2 > 0.0
            pe += 1/2*k*Δ2^2
        end
    end
    pe
end

function energy(q,q̇,rbs,vss)
    q2rbstate!(rbs,q,q̇)
    update_forces!(vss,rbs)
    ke = kineticenergy(rbs)
    pe = potentialenergy(vss)
    ke + pe
end

end # module
