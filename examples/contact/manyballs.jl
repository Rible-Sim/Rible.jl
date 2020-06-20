using LinearAlgebra
using StaticArrays
using Rotations
using Parameters
using DifferentialEquations
using ForwardDiff
using Revise
using TensegrityRobotSimulator
TRS = TensegrityRobotSimulator

function ballbody(radius,center)
    ap1 = TRS.AnchorPoint(0.0,0.0,-0.5)
    ap2 = TRS.AnchorPoint(0.0,0.0,0.5)

    CoM1 = @SVector zeros(3)
    rb1prop = TRS.RigidBodyProperty(:rb1,:generic,1.0,
                SMatrix{3,3}(1.0I),
                CoM1,
                SVector(ap1,ap2))

    r1 = center
    R1 = RotX(1.0π)
    ṙ1 = [0.0,0.0,0.0]
    ω1 = [0.0,0.0,0.0]

    rb1state = TRS.RigidBodyState(rb1prop,r1,R1,ṙ1,ω1)
    rb1 = TRS.RigidBody(rb1prop,rb1state,r=radius)
end
ball1 = ballbody(1.0,[0.0,0.0,0.0])
TRS.sphere_object(1.0,[0.0,0.0,0.0])
n = 10
random_radii = rand(n).+0.5
random_positions = [[(rand()-0.5)*10,(rand()-0.5)*10,10.] for i = 1:n]
balls = [ballbody(random_radii[i],random_positions[i]) for i = 1:n]
tgsys = TRS.TensegritySystem(balls,[],[],[])
manyballstate = TRS.multibodystate(tgsys.rigidbodies)
const g = 9.8
function state2bodystate(state,i,ncstate)
    rbstate = @view state[(i-1)*ncstate+1:(i-1)*ncstate+ncstate]
    c1 = @view rbstate[1:3]
    c2 = @view rbstate[4:7]
    c3 = @view rbstate[8:10]
    c4 = @view rbstate[11:13]
    c1,c2,c3,c4
end
function state2bodystate!(c1,c2,c3,c4,state,i,ncstate)
    rbstate = @view state[(i-1)*ncstate+1:(i-1)*ncstate+ncstate]
    c1 .= @view rbstate[1:3]
    c2 .= @view rbstate[4:7]
    c3 .= @view rbstate[8:10]
    c4 .= @view rbstate[11:13]
end

function hat(ω)
    norm_ω = norm(ω)
    ret = Vector{eltype(ω)}(undef,4)
    ret[1] = cos(norm_ω/2)
    ret[2:4] = sin(norm_ω/2).*ω/norm_ω
    ret
end
function step_manyballs!(u1,u,tgsys,dt)
    rbs = tgsys.rigidbodies
    sts = tgsys.strings
    cnt = tgsys.connectivity
    nbody = length(rbs)
    ncstate = 13
    # State pass
    for i in eachindex(rbs)
        @unpack prop,state = rbs[i]
        @unpack mass,inertia,aps = prop
        @unpack r,R,ṙ,ω,coords = state
        @unpack x,q,p,L = coords
        state2bodystate!(x,q,p,L,u,i,ncstate)
        r .= x
        R .= SMatrix(Quat(q[1],q[2],q[3],q[4]))
        ṙ .= p/mass
        invJt = R*inv(inertia)*transpose(R)
        ω .= invJt*L

        for ip in eachindex(aps)
            state.p[ip] .= point_position(prop,ip,r,R)
        end
    end
    # External Force Pass
    for i in eachindex(rbs)
        rbstate = rbs[i].state
        for ip in eachindex(rbs[i].prop.aps)
            rbstate.Fanc[ip] .= 0.0
            rbstate.τanc[ip] .= 0.0
        end
        rbstate.F .= [0.0,0.0,-g]
        rbstate.τ .= 0.0
    end

    # Collision detection and modeling
    # Advance the velocities
    for i in eachindex(rbs)
        @unpack mass,inertia = rbs[i].prop
        @unpack r,R,ṙ,ω,F,τ,coords = rbs[i].state
        @unpack x,q,p,L = coords
        x1,q1,p1,L1 = state2bodystate(u1,i,ncstate)
        v = ṙ
        # Explicit Euler
        v1 = v + F*dt
        p1 .= v1*mass
        L1 .= L + τ*dt
    end
    # Contact resolution
    # Advance the positions
    for i in eachindex(rbs)
        @unpack mass,inertia = rbs[i].prop
        @unpack r,R,ṙ,ω,coords = rbs[i].state
        @unpack x,q,p,L = coords
        x1,q1,p1,L1 = state2bodystate(u1,i,ncstate)
        v1 = p1/mass
        x1 .= x + dt*v1
        invJt = R*inv(inertia)*transpose(R)
        ω1 = invJt*L1
        q1 .= q + dt*0.5*quat_multiply([0.0,ω1[1],ω1[2],ω1[3]],q)
    end
end

u,du = TRS.multibodystate(tgsys.rigidbodies)
function simulate(tgsys,tspan,dt)
    u,u1 = TensegrityRobotSimulator.multibodystate(tgsys.rigidbodies)
    sol = [u]
    t = tspan[1]
    while t < tspan[2]
        step_manyballs!(u1,u,tgsys,dt)
        push!(sol,copy(u1))
        copy!(u,u1)
        t += dt
    end
    sol
end
solution = simulate(tgsys,(0.0,1.0),0.01)
