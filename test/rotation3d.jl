using LinearAlgebra
using StaticArrays
using Parameters
using Rotations
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot
import PyPlot; const plt = PyPlot

m = rand()
CoM = rand(3)
inertia = Symmetric(rand(3,3))
r0 = rand(3)
ṙ0 = rand(3)
R0 = rand(RotMatrix{3})
ω0 = rand(3)

prop = TR.RigidBodyProperty(1,true,m,inertia,
                    SVector{3}(CoM),[SVector{3}(CoM)])

ri = @SVector rand(3)
rj = @SVector rand(3)
rk = @SVector rand(3)
rl = @SVector rand(3)
u = rj - ri
v = rk - ri
w = rl - ri
# Maybe TR.NaturalCoordinates.BasicPoints3D(ri,ro)
bpsv1,q1 = TR.NaturalCoordinates.BP1P3V(ri,r0,R0)
bpsv2,q2 = TR.NaturalCoordinates.BP2P2V(ri,rj,r0,R0)
bpsv3,q3 = TR.NaturalCoordinates.BP3P1V(ri,rj,rk,r0,R0)
bpsv4,q4 = TR.NaturalCoordinates.BP4P(ri,rj,rk,rl,r0,R0)
bpsv1,q1,q̇1 = TR.NaturalCoordinates.BP1P3V(ri,r0,R0,ṙ0,ω0)
bpsv2,q2,q̇2 = TR.NaturalCoordinates.BP2P2V(ri,rj,r0,R0,ṙ0,ω0)
bpsv3,q3,q̇3 = TR.NaturalCoordinates.BP3P1V(ri,rj,rk,r0,R0,ṙ0,ω0)
bpsv4,q4,q̇4 = TR.NaturalCoordinates.BP4P(ri,rj,rk,rl,r0,R0,ṙ0,ω0)
# bps1 = TR.NaturalCoordinates.BasicPoints1P3V(r̄i,ū,v̄,w̄)
# bps2 = TR.NaturalCoordinates.BasicPoints2P2V(r̄i,r̄j,v̄,w̄)
# bps3 = TR.NaturalCoordinates.BasicPoints3P1V(r̄i,r̄j,r̄k,w̄)
# bps4 = TR.NaturalCoordinates.BasicPoints4P(r̄i,r̄j,r̄k,r̄l)
@code_warntype TR.NaturalCoordinates.BP1P3V(ri,r0,R0)
@code_warntype TR.NaturalCoordinates.BP2P2V(ri,rj,r0,R0)
@code_warntype TR.NaturalCoordinates.BP3P1V(ri,rj,rk,r0,R0)
@code_warntype TR.NaturalCoordinates.BP4P(ri,rj,rk,rl,r0,R0)
@code_warntype TR.NaturalCoordinates.BP1P3V(ri,r0,R0,ṙ0,ω0)
@code_warntype TR.NaturalCoordinates.BP2P2V(ri,rj,r0,R0,ṙ0,ω0)
@code_warntype TR.NaturalCoordinates.BP3P1V(ri,rj,rk,r0,R0,ṙ0,ω0)
@code_warntype TR.NaturalCoordinates.BP4P(ri,rj,rk,rl,r0,R0,ṙ0,ω0)
cache4 = TR.NaturalCoordinatesCache(prop,bps4)
cache3 = TR.NaturalCoordinatesCache(prop,bps3)
cache2 = TR.NaturalCoordinatesCache(prop,bps2)
cache1 = TR.NaturalCoordinatesCache(prop,bps1)
state4 = TR.RigidBodyState(prop,bpsv4,r0,R0,ṙ0,ω0,q4,q̇4)
state3 = TR.RigidBodyState(prop,bpsv3,r0,R0,ṙ0,ω0,q3,q̇3)
state2 = TR.RigidBodyState(prop,bpsv2,r0,R0,ṙ0,ω0,q2,q̇2)
state1 = TR.RigidBodyState(prop,bpsv1,r0,R0,ṙ0,ω0,q1,q̇1)

@code_warntype TR.RigidBodyState(prop,bpsv4,r0,R0,ṙ0,ω0,q4,q̇4)
@code_warntype TR.RigidBodyState(prop,bpsv3,r0,R0,ṙ0,ω0,q3,q̇3)
@code_warntype TR.RigidBodyState(prop,bpsv2,r0,R0,ṙ0,ω0,q2,q̇2)
@code_warntype TR.RigidBodyState(prop,bpsv1,r0,R0,ṙ0,ω0,q1,q̇1)
# @code_warntype TR.RigidBodyState(prop,bps1,r0,R0,ṙ0,ω0)
#
# @code_warntype TR.NaturalCoordinatesCache(prop,bps1)

# state4 = TR.RigidBodyState(prop,bps4,r0,R0,ṙ0,ω0)
# state3 = TR.RigidBodyState(prop,bps3,r0,R0,ṙ0,ω0)
# state2 = TR.RigidBodyState(prop,bps2,r0,R0,ṙ0,ω0)
# state1 = TR.RigidBodyState(prop,bps1,r0,R0,ṙ0,ω0)


TR.kinetic_energy_coords(state1)
TR.kinetic_energy_coords(state2)
TR.kinetic_energy_coords(state3)
TR.kinetic_energy_coords(state4)

# TR.kinetic_energy_coords(statev1)
# TR.kinetic_energy_coords(statev2)
# TR.kinetic_energy_coords(statev3)
# TR.kinetic_energy_coords(statev4)

# F0 = zeros(12)
# F!(F0,q0,q̇0,0.0)
# rb1.state.cache.CG
# rb1.state.cache.Cp[1]
# CG*q0
# CoM
# CG - rb1.state.cache.Cp[1]



es = [TR.kinetic_energy_coords(rb1,q̇) for q̇ in sol.q̇s]
errs = [Φ(q) for q in sol.qs]
errAs = [Φq(q)*q̇ for (q,q̇) in zip(sol.qs,sol.q̇s)]

rb2 = TR.RigidBody(prop,state2)
q0,q̇0,λ0 = TR.get_initial(rb2)
M,Φ,Φq,F! = rbfuncs(rb2)

dt = 0.01
prob = TS.DyProblem(rbfuncs(rb2),q0,q̇0,λ0,(0.0,5.0))
sol = TS.solve(prob,TS.Wendlandt(),dt=dt,ftol=1e-14,verbose=true)
es = [TR.kinetic_energy_coords(rb2,q̇) for q̇ in sol.q̇s]



m = 1.0
# inertia = Symmetric(rand(3,3))
inertia = Matrix(1.0I,3,3)
r0 = zeros(3)
R0 = Matrix(1.0I,3,3)
ṙ0 = rand(3)
ω0 = rand(3)
CoM = zeros(3)
ri = zeros(3)
rj = rand(3)
rk = rand(3)
rl = rand(3)
# bps1,q1,q̇1 = TR.NaturalCoordinates.BP1P3V(ri,r0,R0,ṙ0,ω0)
# bps1,q1,q̇1 = TR.NaturalCoordinates.BP2P2V(ri,rj,r0,R0,ṙ0,ω0)
bps1,q1,q̇1 = TR.NaturalCoordinates.BP3P1V(ri,rj,rk,r0,R0,ṙ0,ω0)
# bps1,q1,q̇1 = TR.NaturalCoordinates.BP4P(ri,rj,rk,rl,r0,R0,ṙ0,ω0)
q1
q̇1

prop = TR.RigidBodyProperty(1,true,m,inertia,
                    SVector{3}(CoM),[SVector{3}(CoM)])
state1 = TR.RigidBodyState(prop,bps1,r0,R0,ṙ0,ω0,q1,q̇1)
TR.kinetic_energy_coords(state1)


rb1 = TR.RigidBody(prop,state1)
q0,q̇0,λ0 = TR.get_initial(rb1)
e0 = TR.kinetic_energy_coords(rb1,q̇0)
function rbfuncs(rb)
    @unpack cache = rb.state
    @unpack M,CG,funcs = cache
    @unpack Φ,Φq = funcs
    f = [0.0,0.0,-9.8]
    function F!(F,q,q̇,t)
        F .= transpose(CG)*f
    end
    M,Φ,Φq,F!,nothing
end
M,Φ,Φq,F! = rbfuncs(rb1)
isapprox(Φ(q0),zeros(6);atol=1e-15)
Φq(q0)*q̇0

dt = 0.01
prob = TS.DyProblem(rbfuncs(rb1),q0,q̇0,λ0,(0.0,0.1))
sol = TS.solve(prob,TS.Wendlandt(),dt=dt,ftol=1e-14,verbose=true)
# sol = TS.solve(prob,TS.ConstSPARK(1),dt=dt,ftol=1e-13,verbose=true)

TR.kinetic_energy_coords(rb1.state)
es = [TR.kinetic_energy_coords(rb1,q̇) for q̇ in sol.q̇s]

plt.plot(es)
