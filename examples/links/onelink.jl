using LinearAlgebra
using StaticArrays
using Rotations
using Parameters
using GeometryTypes, Makie, AbstractPlotting
#using DifferentialEquations
#using ForwardDiff
using Revise
using SPARK
using TensegrityRobotSimulator
const TRS = TensegrityRobotSimulator

function RigidBody(name,r = [0.0,0.0,0.0],
                        R = Matrix(1.0I,3,3),
                        ṙ = [0.0,0.0,0.0],
                        ω = [0.0,0.0,0.0])
    a = 0.5
    h = 1.0
    θ = 2π/3
    l = √(a^2+h^2)
    b = √3a

    # Configuration
    offset = [0.0, 0.0, 0.0]
    ap1 = TRS.AnchorPoint([a, 0.0, 0.0] + offset)
    ap2 = TRS.AnchorPoint([a*cos(θ), a*sin(θ), 0.0] + offset)
    ap3 = TRS.AnchorPoint([a*cos(θ), -a*sin(θ), 0.0] + offset)
    ap4 = TRS.AnchorPoint([0.0, 0.0, -h] + offset)
    CoM = @SVector zeros(3)
    prop = TRS.RigidBodyProperty(name,:generic,1.0,
                SMatrix{3,3}(1.0I),
                CoM,
                SVector(ap1,ap2,ap3,ap4))

    state = TRS.RigidBodyState(prop,r,R,ṙ,ω,Val(:NC))

    TRS.RigidBody(prop,state)
end

rb1 = RigidBody(:rb1,[0.0,0.0,-1.0])
rb2 = RigidBody(:rb1,[0.0,0.0,0.0])
rb3 = RigidBody(:rb1,[0.0,0.0,1.0])
rbs = [rb1,rb2,rb3]

sts = nothing
# global to local
connectivity = nothing
#tgsys = TRS.TensegritySystem(rbs,[s1,s2,s3],[j1,j2],connectivity)
tgsys = TRS.TensegritySystem(rbs,sts,nothing,connectivity)

function tg_spark(tgsys)
    @unpack rigidbodies,strings,joints,connectivity = tgsys
    nbody = length(rigidbodies)
    ninconstraint = nbody*6
    nexconstraint = 0
    nconstraint = ninconstraint + nexconstraint
    nq = nbody*12
    mass_matrix = zeros(nq,nq)
    for (rbid,rb) in enumerate(rigidbodies)
        is = 12*(rbid - 1)
        mass_matrix[is+1:is+12,is+1:is+12] .= rigidbodies[rbid].state.cache.M
    end
    function M!(mm,q)
        mm .= mass_matrix
    end

    # System functions
    function ∂T∂q̇!(p,M,q,q̇)
        mul!(p,mass_matrix,q̇)
    end

    function compute_string_forces!(rbs,sts)
        for st in sts
            rbida,pida = connectivity[st.ida]
            rbidb,pidb = connectivity[st.idb]
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
            TRS.coords2state_kinetic!(rb)
        end
    end
    function F!(F,q,q̇,t)
        distribute_q_to_rbs!(rigidbodies,q,q̇)
        compute_string_forces!(rigidbodies,strings)
        for (rbid,rb) = enumerate(rigidbodies)
            is = 12*(rbid-1)
            F[is+1:is+12] .= rb.state.coords.Q
        end
    end


    function Φ(q)
        ret = Vector{eltype(q)}(undef,nbody*6)
        for (rbid,rb) in enumerate(rigidbodies)
            is = 6*(rbid-1)
            ks = 12*(rbid-1)
            ret[is+1:is+6] .= TRS.NC.Φ(q[ks+1:ks+12])
        end
        ret
    end

    function A(q)
        ret = zeros(eltype(q),6nbody,12nbody)
        for (rbid,rb) in enumerate(rigidbodies)
            is = 6*(rbid-1)
            ks = 12*(rbid-1)
            ret[is+1:is+6,ks+1:ks+12] .= TRS.NC.Φq(q[ks+1:ks+12])
        end
        ret
    end

    A,Φ,∂T∂q̇!,F!,M!,nothing
end
A,Φ,∂T∂q̇!,F!,M!,jacs = tg_spark(tgsys)
function system_variables(tgsys)
    q = vcat([rb.state.coords.q for rb in tgsys.rigidbodies]...)
    q̇ = vcat([rb.state.coords.q for rb in tgsys.rigidbodies]...)
    q,q̇
end
q0,q̇0 = system_variables(tgsys)
λ0 = zeros(6*3)

s = 1
tab = SPARKTableau(s)
tspan = (0.0,0.1)
cache = SPARKCache(36,18,0.0001,tspan,s,(A,Φ,∂T∂q̇!,F!,M!,jacs))
state = SPARKsolve!(q0,q̇0,λ0,cache,tab)
