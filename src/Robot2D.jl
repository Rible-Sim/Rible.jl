module Robot2D
using SparseArrays
using LinearAlgebra
using Parameters
using StaticArrays
using NLsolve
include("PIDController.jl")
include("NC.jl")
greet() = print("Hello World!")

mutable struct SStringState{T}
    restlength::T
    length::T
    lengthdot::T
    tension::T
    direction::MArray{Tuple{2},T,1,2}
end

function SStringState(restlength,length,tension)
    SStringState(restlength,length,zero(length),tension,MVector(1.0,0.0))
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
    ndim::Int
    nbody::Int
    nmovablebody::Int
    mvbodyindex::Vector{Int}
    nfixbody::Int
    fixbodyindex::Vector{Int}
    nstring::Int
    npoints::Int
    rigidbodies::Vector{BodyType}
    strings::Vector{StringType}
    actuators::Vector{ActuatorType}
    connectivity::ConnectType
end

function Structure2D(rbs,ss,acs,cnt)
    ndim = 2
    nbody = length(rbs)
    mvbodyindex = [i for i in eachindex(rbs) if rbs[i].prop.movable]
    nmvbody = length(mvbodyindex)
    fixbodyindex = [i for i in eachindex(rbs) if !rbs[i].prop.movable]
    nfixbody = length(fixbodyindex)
    nstring = length(ss)
    npoints = 0
    for (rbid,rb) in enumerate(rbs)
        npoints += rb.prop.number_aps
    end
    Structure2D(ndim,nbody,nmvbody,mvbodyindex,
                    nfixbody,fixbodyindex,
                    nstring,npoints,
                    rbs,ss,acs,cnt)
end

struct TGRobot2D{ST,CT}
    st2d::ST
    ctrl::CT
end

function lengthdir(v)
    l = norm(v)
    τ = v/l
    l,τ
end

function reset_forces!(st2d::Structure2D)
    reset_forces!(st2d.rigidbodies)
end
function reset_forces!(rbs)
    for (rbid,rb) in enumerate(rbs)
        for f in rb.state.Fanc
            f .= 0.0
        end
        rb.state.F .= 0.0
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
        sstate.direction = τ
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

function generate_forces!(rbs)
    for (rbid,rb) in enumerate(rbs)
        @unpack state = rb
        @unpack Fanc = state
        @unpack Q,Cp,CG = state.auxs
        Q .= 0.0
        for (pid,f) in enumerate(Fanc)
            Q .+= transpose(Cp[pid])*f
        end
        Q .+= transpose(CG)*state.F
    end
end

function assemble_forces!(F,st2d)
    rbs = st2d.rigidbodies
    @unpack body2q = st2d.connectivity
    generate_forces!(rbs)
    F .= 0.0
    for (rbid,rb) in enumerate(rbs)
        pindex = body2q[rbid]
        F[pindex] .+= rb.state.auxs.Q
    end
end

function apply_gravity!(st2d)
    rbs = st2d.rigidbodies
    for (rbid,rb) in enumerate(rbs)
        @unpack state = rb
        g = 9.8
        rb.state.F .= [0.0,-g]
    end
end

include("inverse.jl")
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

function build_body2q(rbs::Vector{RigidBody2D{T,CT,AT}}) where {T,CT,AT}
    bps = Vector{Vector{T}}()
    bp_number = Vector{Int}()
    push!(bp_number,0)
    body2q = Vector{Vector{Int}}()
    for (rbid,rb) in enumerate(rbs)
        @unpack state = rb
        xi,yi,xj,yj = state.coords.q
        bp1 = [xi,yi]
        bp2 = [xj,yj]
        bp1_find = findfirst(x->x==bp1,bps)
        if bp1_find === nothing
            push!(bps,bp1)
            push!(bp_number,bp_number[end]+1)
            bp1_number = bp_number[end]
        else
            bp1_number = bp1_find
        end
        bp2_find = findfirst(x->x==bp2,bps)
        if bp2_find === nothing
            push!(bps,bp2)
            push!(bp_number,bp_number[end]+1)
            bp2_number = bp_number[end]
        else
            bp2_number = bp2_find
        end
        push!(body2q,[2bp1_number-1,2bp1_number,
                      2bp2_number-1,2bp2_number])
    end
    body2q
end

function build_massmatrix(st2d::Structure2D)
    rbs = st2d.rigidbodies
    body2q = st2d.connectivity.body2q
    build_massmatrix(rbs,body2q)
end
function build_massmatrix(rbs,body2q)
    nq = body2q[end][end]
    mass_matrix = zeros(nq,nq)
    for (rbid,rb) in enumerate(rbs)
        pindex = body2q[rbid]
        mass_matrix[pindex,pindex] .+= rbs[rbid].state.auxs.M
    end
    mass_matrix
end

function get_nconstraint(st2d)
    @unpack nbody,nfixbody = st2d
    ninconstraint = nbody
    nexconstraint = 3nfixbody
    nconstraint = ninconstraint + nexconstraint
end

function build_Φ(st2d)
    rbs = st2d.rigidbodies
    q0,q̇0 = get_q(st2d)
    @unpack body2q = st2d.connectivity
    nfixbody = st2d.nfixbody
    nconstraint = get_nconstraint(st2d)
    function inner_Φ(q)
        ret = Vector{eltype(q)}(undef,nconstraint)
        for (rbid,rb) in enumerate(rbs)
            pindex = body2q[rbid]
            ret[3nfixbody+rbid] = rb.state.auxs.Φ(q[pindex])
        end
        for (fixid,rbid) in enumerate(st2d.fixbodyindex)
            pindex = body2q[rbid]
            ret[3(fixid-1)+1:3fixid] .= q[pindex[[1,2,4]]] - q0[pindex[[1,2,4]]]
        end
        ret
    end
end

function build_A(st2d)
    rbs = st2d.rigidbodies
    @unpack body2q = st2d.connectivity
    nfixbody = st2d.nfixbody
    nconstraint = get_nconstraint(st2d)
    nq = body2q[end][end]
    function inner_Φ(q)
        ret = zeros(eltype(q),nconstraint,nq)
        for (rbid,rb) in enumerate(rbs)
            pindex = body2q[rbid]
            ret[3nfixbody+rbid,pindex] .= rb.state.auxs.Φq(q[pindex])
        end
        for (fixid,rbid) in enumerate(st2d.fixbodyindex)
            ret[3(fixid-1)+1,1:4] .= [1.0,0.0,0.0,0.0]
            ret[3(fixid-1)+2,1:4] .= [0.0,1.0,0.0,0.0]
            ret[3(fixid-1)+3,1:4] .= [0.0,0.0,0.0,1.0]
        end
        ret
    end
end

function get_q(st2d)
    rbs = st2d.rigidbodies
    @unpack body2q = st2d.connectivity
    nq = body2q[end][end]
    q = zeros(nq)
    q̇ = zeros(nq)
    for (rbid,rb) in enumerate(rbs)
        pindex = body2q[rbid]
        q[pindex] .= rbs[rbid].state.coords.q
        q̇[pindex] .= rbs[rbid].state.coords.q̇
    end
    return q,q̇
end

function get_initial(st2d)
    q0,q̇0 = get_q(st2d)
    λ0 = zeros(get_nconstraint(st2d))
    q0,q̇0,λ0
end

function actuate!(st2d::Structure2D,u)
    for (actuator,u) in zip(st2d.actuators,u)
        actuate!(actuator,u)
    end
end
function actuate!(actuator::Actuator,u)
    @unpack strings = actuator
    str1 = strings[1]
    str1.state.restlength = str1.original_restlength + u
    str2 = strings[2]
    str2.state.restlength = str2.original_restlength - u
end
function get_u(st2d)
end
end # module
