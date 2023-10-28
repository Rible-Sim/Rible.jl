using LinearAlgebra
using StaticArrays
using Rotations
using Parameters
#using DifferentialEquations
#using ForwardDiff
using Revise
using SPARK
using Ribleimulator
const TR = Ribleimulator

function RigidBody(id;r = [0.0,0.0,-1.0],
                      R = Matrix(1.0I,3,3),
                      ṙ = [0.0,0.0,0.0],
                      ω = [0.0,0.0,0.0],
                      movable = true)
    a = 0.5 #m
    h = 1.0 #m
    θ = 2π/3
    l = √(a^2+h^2)
    b = √3a

    # Configuration
    offset = [0.0, 0.0, h]
    ap1 = SVector{3}([a, 0.0, 0.0] + offset)
    ap2 = SVector{3}([a*cos(θ), a*sin(θ), 0.0] + offset)
    ap3 = SVector{3}([a*cos(θ), -a*sin(θ), 0.0] + offset)
    ap4 = SVector{3}([0.0, 0.0, -h] + offset)

    prop = RB.RigidBodyProperty(id,movable,mass,
                SMatrix{3,3}(inertia),
                SVector(mass_locus...),
                [ap1,ap2,ap3,ap4])
    cache = RB.NCFCache(prop)
    #state = RB.RigidBodyState(prop,r,R,ṙ,ω,Val(:NC))
    prop,cache
    #RB.RigidBody(prop,state)
end

mass = 1.0 #kg
#inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-8 # N/m^2
inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-1
mass_locus = [0.0, 0.0, 17.56] .* 1e-4 # m
rb1prop,rb1cache = RigidBody(1,mass = mass, inertia = inertia, mass_locus = mass_locus,
                     r = [0.0,0.0, -1.0], movable = false)
rb2 = RigidBody(:rb2,mass = mass, inertia = inertia, mass_locus = mass_locus,
                     r = [0.0,0.0,  -0.1], R = Matrix(RotX(0.45)))
rb3 = RigidBody(:rb3,mass = mass, inertia = inertia, mass_locus = mass_locus,
                     r = [0.0,-0.3,  0.7], R = Matrix(RotX(0.8)))
rbs = [rb1,rb2,rb3]

rb1 = RigidBody(:rb1,mass = mass, inertia = inertia, mass_locus = mass_locus,
                     r = [0.0,0.0, -1.0])
rb2 = RigidBody(:rb2,mass = mass, inertia = inertia, mass_locus = mass_locus,
                     r = [0.0,0.0,  -0.1], R = Matrix(RotX(-0.5)))
rb3 = RigidBody(:rb3,mass = mass, inertia = inertia, mass_locus = mass_locus,
                     r = [0.0,0.3,  0.7], R = Matrix(RotX(-1.0)))
rbs = [rb1,rb2,rb3]

rb1 = RigidBody(:rb1,mass = mass, inertia = inertia, mass_locus = mass_locus,
                     r = [0.0,0.0, -1.0])
rb2 = RigidBody(:rb2,mass = mass, inertia = inertia, mass_locus = mass_locus,
                     r = [0.0,0.0,  -0.1])
rb3 = RigidBody(:rb3,mass = mass, inertia = inertia, mass_locus = mass_locus,
                     r = [0.0,0.0,  0.8])
rbs = [rb1,rb2,rb3]

mvbodyindex = [i for i in eachindex(rbs) if rbs[i].prop.movable]
mvrbs = [rbs[i] for i in mvbodyindex]

rbv = RB.RBVector(rbs,mvbodyindex,mvrbs)


s0_inner =0.2 #m
k_inner = 10.0 #N/m
s0_outer = 0.8 #m
k_outer = 10.0 #N/m
s1 = RB.LinearString(k_inner,s0_inner,1,8)
s2 = RB.LinearString(k_inner,s0_inner,2,8)
s3 = RB.LinearString(k_inner,s0_inner,3,8)
s4 = RB.LinearString(k_outer,s0_outer,4,8)
s5 = RB.LinearString(k_outer,s0_outer,1,5)
s6 = RB.LinearString(k_outer,s0_outer,2,6)
s7 = RB.LinearString(k_outer,s0_outer,3,7)

s8 = RB.LinearString(k_inner,s0_inner,5,12)
s9 = RB.LinearString(k_inner,s0_inner,6,12)
s10 = RB.LinearString(k_inner,s0_inner,7,12)
s11 = RB.LinearString(k_outer,s0_outer,8,12)
s12 = RB.LinearString(k_outer,s0_outer,5,9)
s13 = RB.LinearString(k_outer,s0_outer,6,10)
s14 = RB.LinearString(k_outer,s0_outer,7,11)

sts = [s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14]
# global to local
connectivity = [
    [1,1], [1,2], [1,3], [1,4],
    [2,1], [2,2], [2,3], [2,4],
    [3,1], [3,2], [3,3], [3,4],
    ]
#tgsys = RB.TensegritySystem(rbs,[s1,s2,s3],[j1,j2],connectivity)
tgsys = RB.TensegritySystem(rbv,sts,nothing,connectivity)


function tg_spark(tgsys)
    @unpack rigidbodies,strings,joints,connectivity = tgsys
    ntotalbody = length(rigidbodies)
    mvrigidbodies = rigidbodies.movables
    nbodies = length(mvrigidbodies)
    ninconstraint = nbodies*6
    nexconstraint = 0
    nconstraint = ninconstraint + nexconstraint
    nq = nbodies*12
    mass_matrix = zeros(nq,nq)
    for (bodyid,rb) in enumerate(mvrigidbodies)
        is = 12*(bodyid - 1)
        mass_matrix[is+1:is+12,is+1:is+12] .= mvrigidbodies[bodyid].state.cache.M
    end
    function M!(mm,q)
        mm .= mass_matrix
    end

    # System functions
    function ∂T∂q̇!(p,M,q,q̇)
        mul!(p,mass_matrix,q̇)
    end


    function F!(F,q,q̇,t)
        RB.reset_forces!.(rigidbodies)
        RB.distribute_q_to_rbs!(mvrigidbodies,q,q̇)
        RB.compute_string_forces!(tgsys)
        for (bodyid,rb) = enumerate(mvrigidbodies)
            is = 12*(bodyid-1)
            F[is+1:is+12] .= body.state.coords.Q
        end
    end


    function Φ(q)
        ret = Vector{eltype(q)}(undef,nbodies*6)
        for (bodyid,rb) in enumerate(mvrigidbodies)
            is = 6*(bodyid-1)
            ks = 12*(bodyid-1)
            ret[is+1:is+6] .= RB.NC.Φ(q[ks+1:ks+12])
        end
        ret
    end

    function A(q)
        ret = zeros(eltype(q),6nbody,12nbody)
        for (bodyid,rb) in enumerate(mvrigidbodies)
            is = 6*(bodyid-1)
            ks = 12*(bodyid-1)
            ret[is+1:is+6,ks+1:ks+12] .= RB.NC.Φq(q[ks+1:ks+12])
        end
        ret
    end

    A,Φ,∂T∂q̇!,F!,M!,nothing
end
A,Φ,∂T∂q̇!,F!,M!,jacs = tg_spark(tgsys)
function system_variables(tgsys)
    q = vcat([body.state.coords.q for rb in tgsys.rigidbodies.movables]...)
    q̇ = vcat([body.state.coords.q̇ for rb in tgsys.rigidbodies.movables]...)
    q,q̇
end
q0,q̇0 = system_variables(tgsys)
λ0 = zeros(6*2)

s = 1
tab = SPARKTableau(s)
tspan = (0.0,10.1)
cache = SPARKCache(24,12,0.01,tspan,s,(A,Φ,∂T∂q̇!,F!,M!,jacs))
state = SPARKsolve!(q0,q̇0,λ0,cache,tab)



function potential_energy(tgsys)
    pe = 0.0
    for (stid,st) in enumerate(tgsys.strings)
        pe_st = RB.potential_energy(st)
        pe += pe_st
    end
    pe
end
potential_energy(tgsys)
function kinetic_energy(tgsys)
    ke = 0.0
    for (bodyid,rb) in enumerate(tgsys.rigidbodies)
        ke_rb = RB.kinetic_energy(body)
        ke += ke_rb
    end
    ke
end
kinetic_energy(tgsys)
function energy(tgsys,q,q̇)
    RB.reset_forces!.(tgsys.rigidbodies)
    RB.distribute_q_to_rbs!(tgsys.rigidbodies.movables,q,q̇)
    RB.compute_string_forces!(tgsys)
    ke = kinetic_energy(tgsys)
    pe = potential_energy(tgsys)
    ke,pe,ke + pe
end
energy(tgsys,q0,q̇0)
state
ks,ps,restitution_coefficients = [energy(tgsys,state.qs[it],state.q̇s[it]) for it in eachindex(state.ts)]
energy(tgsys,state.qs[1],state.q̇s[1])
tgsys.rigidbodies[1].state.rpṡ
