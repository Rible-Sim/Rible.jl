using LinearAlgebra
using StaticArrays
using Parameters
using Rotations
using Revise
using TensegrityRobots; const TR = TensegrityRobots
import PyPlot; const plt = PyPlot

m = rand()
r̄g = rand(3)
r̄p = rand(3)
inertia = SMatrix{3,3}(Matrix(1.0I,3,3))
r0 = rand(3)
ṙ0 = rand(3)
R0 = rand(RotMatrix{3})
ω0 = rand(3)
prop = TR.RigidBodyProperty(1,true,m,inertia,
                    SVector{3}(r̄g),[SVector{3}(r̄p)])

ri = @SVector rand(3)
rj = @SVector rand(3)
rk = @SVector rand(3)
rl = @SVector rand(3)
u = rj - ri
v = rk - ri
w = rl - ri
lncs1,q1 = TR.NaturalCoordinates.NC1P3V(ri,r0,R0)
lncs2,q2 = TR.NaturalCoordinates.NC2P2V(ri,rj,r0,R0)
lncs3,q3 = TR.NaturalCoordinates.NC3P1V(ri,rj,rk,r0,R0)
lncs4,q4 = TR.NaturalCoordinates.NC4P(ri,rj,rk,rl,r0,R0)
lncs1,q1,q̇1 = TR.NaturalCoordinates.NC1P3V(ri,r0,R0,ṙ0,ω0)
lncs2,q2,q̇2 = TR.NaturalCoordinates.NC2P2V(ri,rj,r0,R0,ṙ0,ω0)
lncs3,q3,q̇3 = TR.NaturalCoordinates.NC3P1V(ri,rj,rk,r0,R0,ṙ0,ω0)
lncs4,q4,q̇4 = TR.NaturalCoordinates.NC4P(ri,rj,rk,rl,r0,R0,ṙ0,ω0)

# @code_warntype TR.NaturalCoordinates.NC1P3V(ri,r0,R0)
# @code_warntype TR.NaturalCoordinates.NC2P2V(ri,rj,r0,R0)
# @code_warntype TR.NaturalCoordinates.NC3P1V(ri,rj,rk,r0,R0)
# @code_warntype TR.NaturalCoordinates.NC4P(ri,rj,rk,rl,r0,R0)
# @code_warntype TR.NaturalCoordinates.NC1P3V(ri,r0,R0,ṙ0,ω0)
# @code_warntype TR.NaturalCoordinates.NC2P2V(ri,rj,r0,R0,ṙ0,ω0)
# @code_warntype TR.NaturalCoordinates.NC3P1V(ri,rj,rk,r0,R0,ṙ0,ω0)
# @code_warntype TR.NaturalCoordinates.NC4P(ri,rj,rk,rl,r0,R0,ṙ0,ω0)
cache4 = TR.NaturalCoordinatesCache(prop,lncs4,q4)
cache3 = TR.NaturalCoordinatesCache(prop,lncs3,q3)
cache2 = TR.NaturalCoordinatesCache(prop,lncs2,q2)
cache1 = TR.NaturalCoordinatesCache(prop,lncs1,q1)
state4 = TR.RigidBodyState(prop,lncs4,r0,R0,ṙ0,ω0,q4,q̇4)
state3 = TR.RigidBodyState(prop,lncs3,r0,R0,ṙ0,ω0,q3,q̇3)
state2 = TR.RigidBodyState(prop,lncs2,r0,R0,ṙ0,ω0,q2,q̇2)
state1 = TR.RigidBodyState(prop,lncs1,r0,R0,ṙ0,ω0,q1,q̇1)
#
# @code_warntype TR.RigidBodyState(prop,lncs4,r0,R0,ṙ0,ω0,q4,q̇4)
# @code_warntype TR.RigidBodyState(prop,lncs3,r0,R0,ṙ0,ω0,q3,q̇3)
# @code_warntype TR.RigidBodyState(prop,lncs2,r0,R0,ṙ0,ω0,q2,q̇2)
# @code_warntype TR.RigidBodyState(prop,lncs1,r0,R0,ṙ0,ω0,q1,q̇1)

TR.kinetic_energy_coords(state1)
TR.kinetic_energy_coords(state2)
TR.kinetic_energy_coords(state3)
TR.kinetic_energy_coords(state4)


rb1 = TR.RigidBody(prop,state1)
rb2 = TR.RigidBody(prop,state2)
rb3 = TR.RigidBody(prop,state3)
rb4 = TR.RigidBody(prop,state4)

function rbfuncs(rb)
    @unpack cache = rb.state
    @unpack M,Cg,funcs = cache
    @unpack Φ,Φq = funcs
    f = [0.0,0.0,-9.8]
    function F!(F,q,q̇,t)
        # F .= transpose(Cg)*f
        # F .= 0.0
    end
    M,Φ,Φq,F!,nothing
end

dt = 0.01

q0,q̇0,λ0 = TR.get_initial(rb1)
prob1 = TR.DyProblem(rbfuncs(rb1),q0,q̇0,λ0,(0.0,5.0))
sol1 = TR.solve(prob1,TR.Zhong06(),dt=dt,ftol=1e-14,verbose=true)

q0,q̇0,λ0 = TR.get_initial(rb2)
prob2 = TR.DyProblem(rbfuncs(rb2),q0,q̇0,λ0,(0.0,5.0))
sol2 = TR.solve(prob2,TR.Zhong06(),dt=dt,ftol=1e-14,verbose=true)

q0,q̇0,λ0 = TR.get_initial(rb3)
prob3 = TR.DyProblem(rbfuncs(rb3),q0,q̇0,λ0,(0.0,5.0))
sol3 = TR.solve(prob3,TR.Zhong06(),dt=dt,ftol=1e-14,verbose=true)

q0,q̇0,λ0 = TR.get_initial(rb4)
prob4 = TR.DyProblem(rbfuncs(rb4),q0,q̇0,λ0,(0.0,5.0))
sol4 = TR.solve(prob4,TR.Zhong06(),dt=dt,ftol=1e-14,verbose=true)


es = [TR.kinetic_energy_coords(rb1,q̇) for q̇ in sol1.q̇s]
plt.plot(es)

plt.plot([q[1] for q in sol1.qs])
plt.plot([q[1] for q in sol2.qs])
plt.plot([q[1] for q in sol3.qs])
plt.plot([q[1] for q in sol4.qs])

plt.plot([q[2] for q in sol1.qs])
plt.plot([q[2] for q in sol2.qs])
plt.plot([q[2] for q in sol3.qs])
plt.plot([q[2] for q in sol4.qs])

plt.plot([q[3] for q in sol1.qs])
plt.plot([q[3] for q in sol2.qs])
plt.plot([q[3] for q in sol3.qs])
plt.plot([q[3] for q in sol4.qs])

m = 1.0
# inertia = Symmetric(rand(3,3))
#inertia = Matrix(Diagonal([0.04,0.0008,0.002]))
inertia = Matrix(Diagonal([0.8,0.4,0.4])).*1e-3
# inertia = Matrix(Diagonal([1.0,1.0,1.0]))
r0 = zeros(3)
ṙ0 = zeros(3)
R0 = [0 1 0;
-0.500000000000000 0 0.866025403784439;
0.866025403784439 0 0.500000000000000]
ω_bar0 = [40π,0,0]
ω0 = R0*ω_bar0

rG = [0,-0.0200000000000000,0.0346410161513775]
r̄G = inv(R0)*(rG-r0)

ri = zeros(3)
rj = [0,-0.0200000000000000,0.0346410161513775]
# rj = R0*[1.0,0.0,0.0]
rk = R0*[0.0,1.0,0.0]
rl = R0*[0.0,0.0,1.0]
# r̄i = @SVector zeros(3)
# r̄j = SVector{3}(inv(R0)*(rj-ri))
# v̄ = SVector(0.0,1.0,0.0)
# w̄ = SVector(0.0,0.0,1.0)
# bps2 = TR.NaturalCoordinates.LocalNaturalCoordinates2P2V(r̄i,r̄j,v̄,w̄)
# q0 = vcat(ri,rj,R0*v̄,R0*w̄)
# q̇0 = vcat(ω0×ri,ω0×rj,ω0×q0[7:9],ω0×q0[10:12])
# cf2 = TR.NaturalCoordinates.CoordinateFunctions(bps2)
# M = TR.NaturalCoordinates.make_M(cf2,m,SMatrix{3,3}(inertia),r̄g)
# bps1,q1,q̇1 = TR.NaturalCoordinates.BP1P3V(ri,r0,R0,ṙ0,ω0)
bps1,q1,q̇1 = TR.NaturalCoordinates.BP2P2V(ri,rj,r0,R0,ṙ0,ω0)
# bps1,q1,q̇1 = TR.NaturalCoordinates.BP3P1V(ri,rj,rk,r0,R0,ṙ0,ω0)
# bps1,q1,q̇1 = TR.NaturalCoordinates.BP4P(ri,rj,rk,rl,r0,R0,ṙ0,ω0)
# q1
# q̇1
# 1/2*transpose(q̇0)*M*q̇0
prop = TR.RigidBodyProperty(1,true,m,inertia,
                    SVector{3}(r̄g),[SVector{3}(r̄g)];constrained=true)
state1 = TR.RigidBodyState(prop,bps1,r0,R0,ṙ0,ω0,q1,q̇1,[1,2,3])
# state1 = TR.RigidBodyState(prop,bps2,r0,R0,ṙ0,ω0,q0,q̇0,[1,2,3])
TR.kinetic_energy_coords(state1)


rb1 = TR.RigidBody(prop,state1)
body2q = [collect(1:12)]
cnt = TR.Connectivity(body2q,nothing)
tgrb1 = TR.TensegrityStructure([rb1],Vector{Int}(),Vector{Int}(),cnt)
q0,q̇0,λ0 = TR.get_initial(tgrb1)
e0 = TR.kinetic_energy_coords(rb1,q̇0)
function dynfuncs(tgstruct,q0)
    M = TR.build_massmatrix(tgstruct)
    Φ = TR.build_Φ(tgstruct,q0)
    A = TR.build_A(tgstruct)
    Cg = Matrix(tgstruct.rigidbodies[1].state.cache.Cg)
    fg = [0.0,0.0,-9.81]
    function F!(F,q,q̇,t)
        F .= transpose(Cg)*fg
    end

    M,Φ,A,F!,nothing
end
M,Φ,Φq,F! = dynfuncs(tgrb1,q0)
M
Φ(q0)
Φq(q0)

F0 = similar(q0)
F!(F0,q0,q̇0,0.0)

# q0,q̇0,λ0  = TR.get_initial(rb1)
# isapprox(Φ(q0),zeros(6);atol=1e-15)
# Φq(q0)*q̇0

dt = 0.002
# prob = TS.DyProblem(rbfuncs(rb1),q0,q̇0,λ0,(0.0,0.0001))
prob = TS.DyProblem(dynfuncs(tgrb1,q0),q0,q̇0,λ0,(0.0,0.2))
# sol = TS.solve(prob,TS.Wendlandt(),dt=dt,ftol=1e-14,verbose=true)
sol = TS.solve(prob,TS.Xu2014(),dt=dt,ftol=1e-14,verbose=true)
sol = TS.solve(prob,TS.Zhong06(),dt=dt,ftol=1e-14,verbose=true)
sol = TS.solve(prob,TS.ConstSPARK(1),dt=dt,ftol=1e-13,verbose=true)

λ1 = sol.λs[2]
sol.qs[2] - refq2
sol.q̇s[2] - refq̇2
sol.ps[2] - refp2
transpose(Φq(q0))*λ1
1/dt*transpose(Φq(sol.qs[2]))*λ1

refq2 = sol.qs[2]
refp2 = sol.ps[2]
refq̇2 = sol.q̇s[2]
refq3 = sol.qs[3]
M*refq̇2



TR.kinetic_energy_coords(rb1.state)
# W_es = [TR.kinetic_energy_coords(rb1,q̇) for q̇ in sol.q̇s]
X_kes = [TR.kinetic_energy_coords(rb1,q̇) for q̇ in sol.q̇s]
X_pes = [TR.gravity_potential_energy(rb1,q) for q in sol.qs]
X_es = X_kes .+ X_pes[1:end-1]
X_es = X_kes .+ X_pes
# plt.plot(es)
# plt.plot(W_es)
X_err = (X_es.-X_es[1])./X_es[1]
plt.plot((X_es.-X_es[1])./X_es[1])
plt.plot(X_kes)
plt.plot(X_err[2:end])

function rb_rotate_funcs(rb)
    @unpack cache = rb.state
    @unpack M,Cg,funcs = cache
    reduceM = M[4:12,4:12]
    f = [0.0,0.0,-9.8]
    function F!(F,q,q̇,t)
        # F .= transpose(Cg)*f
        F .= 0.0
    end

    @unpack r̄i,ū,v̄,w̄ = funcs.lncs
    u_square = ū⋅ū
    v_square = v̄⋅v̄
    w_square = w̄⋅w̄
    uv_dotprod = ū⋅v̄
    uw_dotprod = ū⋅w̄
    vw_dotprod = v̄⋅w̄
    function Φ(q)
        u  = @view q[1:3]
        v  = @view q[4:6]
        w  = @view q[7:9]
        [
        u⋅u - u_square,
        v⋅v - v_square,
        w⋅w - w_square,
        u⋅v - uv_dotprod,
        u⋅w - uw_dotprod,
        v⋅w - vw_dotprod
        ]
    end
    function Φq(q)
        u  = @view q[1:3]
        v  = @view q[4:6]
        w  = @view q[7:9]
        ret = zeros(eltype(q), 6, 9)
        ret[1,1:3] =  2u
        ret[2,4:6] = 2v
        ret[3,7:9] = 2w

        ret[4,1:3] =  v
        ret[4,4:6] =  u

        ret[5,1:3] =  w
        ret[5,7:9] = u

        ret[6,4:6] = w
        ret[6,7:9] = v

        ret
    end

    reduceM,Φ,Φq,F!,nothing
end
function reduce_initial(q0,q̇0,λ0)
    q0[4:12],q̇0[4:12],λ0
end
q0,q̇0,λ0 = reduce_initial(TR.get_initial(rb1)...)
dt = 0.002
prob = TS.DyProblem(rb_rotate_funcs(rb1),q0,q̇0,λ0,(0.0,0.0001))
sol = TS.solve(prob,TS.Wendlandt(),dt=dt,ftol=1e-14,verbose=true)
sol = TS.solve(prob,TS.Xu2014(),dt=dt,ftol=1e-14,verbose=true)

pad_q̇(q̇) = vcat(zeros(3),q̇)
pad_q̇.(sol.q̇s)
W_es = [TR.kinetic_energy_coords(rb1,q̇) for q̇ in pad_q̇.(sol.q̇s)]
X_es = [TR.kinetic_energy_coords(rb1,q̇) for q̇ in pad_q̇.(sol.q̇s)]

plt.plot(es)
plt.plot(W_es)
plt.plot(X_es)
