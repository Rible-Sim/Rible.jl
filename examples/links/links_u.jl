using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using Parameters
#using DifferentialEquations
#using ForwardDiff
using NLsolve
using Cubature
using Makie
AbstractPlotting.__init__()
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

function compute_offset(n,α,Δr,a1,a2)
    R0 = RotZ(α)
    function make_f(r)
        dr(z) = sqrt(1+(2*a2*z+a1)^2)
        function f(Z)
            int,_ = hquadrature(dr,0,Z)
            int-r
        end
    end
    zs = zeros(n)
    for i = 1:n
        result = nlsolve(make_f(Δr*(i-1)),[Δr*(i-1)])
        zs[i:i] .= result.zero
    end
    ys = a2.*zs.^2 .+ a1.*zs
    rs = [R0*[0,y,z] for (y,z) in zip(ys,zs)]
    dy_dz(z) = 2a2*z + a1

    Rs = [
        begin
        θ = atan(dy_dz(z))
        b = [1,0,0]
        τ = [0,sin(θ),cos(θ)]
        n = τ×b
        R = R0*hcat(b,n,τ)
        end
        for z in zs
    ]
    Rs[1] = Matrix(1.0I,3,3)
    rs,Rs
end

function links(n,a1=0.0,a2=0.0,α=0.0)
    nbody = n
    nbp = 4*n

    a = 0.08 #m
    h = 0.1 #m
    θ = 2π/3
    l = √(a^2+h^2)
    b = √3a

    mass = 0.1 #kg
    #inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-8 # N/m^2
    inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-2
    CoM = [0.0, 0.0, 0.0]


    ap1 = SVector{3}([0.0, 0.0, 0.0])
    ap2 = SVector{3}([a, 0.0, h])
    ap3 = SVector{3}([a*cos(θ), a*sin(θ), h])
    ap4 = SVector{3}([a*cos(θ), -a*sin(θ), h])
    aps = [ap1,ap2,ap3,ap4]

    movable = ones(Bool,n)
    movable[1] = 0

    props = [TR.RigidBodyProperty(i,movable[i],mass,
                SMatrix{3,3}(inertia),
                SVector{3}(CoM),
                aps) for i = 1:n]

    dr = 0.5h
    # r = [ [0.0,0.0,0.0+(i-1)*dr] for i = 1:n]
    # R = [Matrix(1.0I,3,3) for i = 1:n]
    r,R = compute_offset(n,α,dr,a1,a2)
    ṙ = [zeros(3) for i = 1:n]
    ω = [zeros(3) for i = 1:n]

    function rigidbody(prop,aps,r,R,ṙ,ω)
        ri,rj,rk,rl = [r+R*ap for ap in aps]
        bps,q,q̇ = TR.NaturalCoordinates.BP4P(ri,rj,rk,rl,r,R,ṙ,ω)
        state = TR.RigidBodyState(prop,bps,r,R,ṙ,ω,q,q̇)
        @show R
        rb = TR.RigidBody(prop,state)
    end
    rbs = [rigidbody(props[i],aps,r[i],R[i],ṙ[i],ω[i]) for i = 1:n]

    nstrings = 7*(n-1)
    stringlenH = 0.6h
    stringlenR = 1.1a
    stringlens = repeat(vcat(fill(stringlenH,4),fill(stringlenR,3)),n-1)
    kH = 1e1
    kR = 1e1
    ks = repeat(vcat(fill(kH,4),fill(kR,3)),n-1)
    c = 0.0
    cs = repeat(fill(c,7),n-1)
    ss = [TR.SString3D(stringlens[i],ks[i],cs[i]) for i = 1:nstrings]
    acs = [TR.Actuator(ss[7(i-1)+2:7i]) for i = 1:n-1]
    bodynq = TR.get_nbodycoords(rbs[1])
    body2q = [(bodynq*(i-1)+1:bodynq*i) for i = 1:n]


    string2ap = Vector{Tuple{TR.ID,TR.ID}}()
    for i = 1:n-1
        for j = 1:4
            push!(string2ap,(TR.ID(i,j),TR.ID(i+1,j)))
        end
        for j = 2:4
            push!(string2ap,(TR.ID(i,j),TR.ID(i+1,1)))
        end
    end
    cnt = TR.Connectivity(body2q,string2ap)
    tg = TR.TensegrityStructure(rbs,ss,acs,cnt)
    TR.update_strings_apply_forces!(tg)
    tg
end
nlink = 5
linkn = links(nlink)
α = 0.0
a2 = 2
a1 = 0.0
compute_offset(nlink,α,0.1*0.5,a1,a2)

reflinkn = links(nlink,a1,a2)


function inverse(tgstruct,refstruct)
    function ikfuncs(tgstruct)

        A = TR.build_A(tgstruct)

        Q̃=TR.build_Q̃(tgstruct)

        function F!(F,u)
            TR.reset_forces!(tgstruct)
            TR.actuate!(tgstruct,u)
            TR.update_strings_apply_forces!(tgstruct)
            F .= Q̃*TR.fvector(tgstruct)
        end

        A,F!
    end
    q0,q̇0,λ0 = TR.get_initial(refstruct)
    TR.distribute_q_to_rbs!(tgstruct,q0,q̇0)
    nu = length(tgstruct.actuators)
    u0 = zeros(nu)
    ikprob = TS.IKProblem(ikfuncs(tgstruct),q0,u0,λ0)
    TR.iksolve(ikprob)
end
u,_ = inverse(linkn,reflinkn)

# @code_warntype links(nlink)
# link3.ndim
# link3.connectivity.body2q[end][end]

q0,q̇0,λ0 = TR.get_initial(linkn)
function dynfuncs(tg)
    M = TR.build_massmatrix(tg)
    Φ = TR.build_Φ(tg)
    A = TR.build_A(tg)
    Q̃ = TR.build_Q̃(tg)
    function F!(F,q,q̇,t)
        TR.reset_forces!(tg)
        TR.distribute_q_to_rbs!(tg,q,q̇)
        TR.update_strings_apply_forces!(tg)
        F .= Q̃*TR.fvector(tg)
        # TR.assemble_forces!(F,tg)
    end
    M,Φ,A,F!,nothing
end
M,Φ,A,F!,Jacs = dynfuncs(linkn)
Φ(q0)
@code_warntype Φ(q0)
A(q0)*q̇0
@code_warntype A(q0)



dt = 0.01
prob = TS.DyProblem(dynfuncs(linkn),q0,q̇0,λ0,(0.0,1.0))
sol = TS.solve(prob,TS.Wendlandt(),dt=dt,ftol=1e-13,verbose=true)
sol = TS.solve(prob,TS.ConstSPARK(1),dt=dt,ftol=1e-12,verbose=true)







rb1prop,rb1cache = RigidBody(1,mass = mass, inertia = inertia, CoM = CoM,
                     r = [0.0,0.0, -1.0], movable = false)
rb2 = RigidBody(:rb2,mass = mass, inertia = inertia, CoM = CoM,
                     r = [0.0,0.0,  -0.1], R = Matrix(RotX(0.45)))
rb3 = RigidBody(:rb3,mass = mass, inertia = inertia, CoM = CoM,
                     r = [0.0,-0.3,  0.7], R = Matrix(RotX(0.8)))
rbs = [rb1,rb2,rb3]

rb1 = RigidBody(:rb1,mass = mass, inertia = inertia, CoM = CoM,
                     r = [0.0,0.0, -1.0])
rb2 = RigidBody(:rb2,mass = mass, inertia = inertia, CoM = CoM,
                     r = [0.0,0.0,  -0.1], R = Matrix(RotX(-0.5)))
rb3 = RigidBody(:rb3,mass = mass, inertia = inertia, CoM = CoM,
                     r = [0.0,0.3,  0.7], R = Matrix(RotX(-1.0)))
rbs = [rb1,rb2,rb3]

rb1 = RigidBody(:rb1,mass = mass, inertia = inertia, CoM = CoM,
                     r = [0.0,0.0, -1.0])
rb2 = RigidBody(:rb2,mass = mass, inertia = inertia, CoM = CoM,
                     r = [0.0,0.0,  -0.1])
rb3 = RigidBody(:rb3,mass = mass, inertia = inertia, CoM = CoM,
                     r = [0.0,0.0,  0.8])
rbs = [rb1,rb2,rb3]

mvbodyindex = [i for i in eachindex(rbs) if rbs[i].prop.movable]
mvrbs = [rbs[i] for i in mvbodyindex]

rbv = TR.RBVector(rbs,mvbodyindex,mvrbs)



sts = [s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14]
# global to local
connectivity = [
    [1,1], [1,2], [1,3], [1,4],
    [2,1], [2,2], [2,3], [2,4],
    [3,1], [3,2], [3,3], [3,4],
    ]
#tgsys = TR.TensegritySystem(rbs,[s1,s2,s3],[j1,j2],connectivity)
tgsys = TR.TensegritySystem(rbv,sts,nothing,connectivity)


function tg_spark(tgsys)
    @unpack rigidbodies,strings,joints,connectivity = tgsys
    ntotalbody = length(rigidbodies)
    mvrigidbodies = rigidbodies.movables
    nbody = length(mvrigidbodies)
    ninconstraint = nbody*6
    nexconstraint = 0
    nconstraint = ninconstraint + nexconstraint
    nq = nbody*12
    mass_matrix = zeros(nq,nq)
    for (rbid,rb) in enumerate(mvrigidbodies)
        is = 12*(rbid - 1)
        mass_matrix[is+1:is+12,is+1:is+12] .= mvrigidbodies[rbid].state.cache.M
    end
    function M!(mm,q)
        mm .= mass_matrix
    end

    # System functions
    function ∂T∂q̇!(p,M,q,q̇)
        mul!(p,mass_matrix,q̇)
    end


    function F!(F,q,q̇,t)
        TR.reset_forces!.(rigidbodies)
        TR.distribute_q_to_rbs!(mvrigidbodies,q,q̇)
        TR.compute_string_forces!(tgsys)
        for (rbid,rb) = enumerate(mvrigidbodies)
            is = 12*(rbid-1)
            F[is+1:is+12] .= rb.state.coords.Q
        end
    end


    function Φ(q)
        ret = Vector{eltype(q)}(undef,nbody*6)
        for (rbid,rb) in enumerate(mvrigidbodies)
            is = 6*(rbid-1)
            ks = 12*(rbid-1)
            ret[is+1:is+6] .= TR.NC.Φ(q[ks+1:ks+12])
        end
        ret
    end

    function A(q)
        ret = zeros(eltype(q),6nbody,12nbody)
        for (rbid,rb) in enumerate(mvrigidbodies)
            is = 6*(rbid-1)
            ks = 12*(rbid-1)
            ret[is+1:is+6,ks+1:ks+12] .= TR.NC.Φq(q[ks+1:ks+12])
        end
        ret
    end

    A,Φ,∂T∂q̇!,F!,M!,nothing
end
A,Φ,∂T∂q̇!,F!,M!,jacs = tg_spark(tgsys)
function system_variables(tgsys)
    q = vcat([rb.state.coords.q for rb in tgsys.rigidbodies.movables]...)
    q̇ = vcat([rb.state.coords.q̇ for rb in tgsys.rigidbodies.movables]...)
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
        pe_st = TR.potential_energy(st)
        pe += pe_st
    end
    pe
end
potential_energy(tgsys)
function kinetic_energy(tgsys)
    ke = 0.0
    for (rbid,rb) in enumerate(tgsys.rigidbodies)
        ke_rb = TR.kinetic_energy(rb)
        ke += ke_rb
    end
    ke
end
kinetic_energy(tgsys)
function energy(tgsys,q,q̇)
    TR.reset_forces!.(tgsys.rigidbodies)
    TR.distribute_q_to_rbs!(tgsys.rigidbodies.movables,q,q̇)
    TR.compute_string_forces!(tgsys)
    ke = kinetic_energy(tgsys)
    pe = potential_energy(tgsys)
    ke,pe,ke + pe
end
energy(tgsys,q0,q̇0)
state
ks,ps,es = [energy(tgsys,state.qs[it],state.q̇s[it]) for it in eachindex(state.ts)]
energy(tgsys,state.qs[1],state.q̇s[1])
tgsys.rigidbodies[1].state.ṙ
