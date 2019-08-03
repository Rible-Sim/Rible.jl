using LinearAlgebra
using StaticArrays
using Rotations
using Parameters
using DifferentialEquations
using Revise
using TensegrityRobotSimulator
TRS = TensegrityRobotSimulator

ap1 = TRS.AnchorPoint([0.0,0.0,-1.0])
ap2 = TRS.AnchorPoint(0.0,0.0,1.0)

CoM1 = @SVector zeros(3)
rb1prop = TRS.RigidBodyProperty(:rb1,:generic,1.0,
            SMatrix{3,3}(1.0I),
            CoM1,
            SVector(ap1,ap2))

r1 = [0.0,0.0,1.0]
R1 = Matrix(1.0I,3,3)
ṙ1 = [0.0,0.0,0.0]
ω1 = [0.0,0.0,0.0]

rb1state = TRS.RigidBodyState(rb1prop,r1,R1,ṙ1,ω1)
rb1Eb = TRS.RigidBodyEberlyCoordinates(rb1prop,rb1state)
rb1 = TRS.RigidBody(rb1prop,rb1state,rb1Eb)

rb2prop = TRS.RigidBodyProperty(:rb2,:generic,2.0,
            SMatrix{3,3}(1.0I),
            CoM1,
            SVector(ap1,ap2))

r2 = [0.0,1.0,0.0]
R2 = RotX(-π/2)
ṙ2 = [0.0,0.0,0.0]
ω2 = [0.0,0.0,0.0]

rb2state = TRS.RigidBodyState(rb2prop,r2,R2,ṙ2,ω2)
rb2Eb = TRS.RigidBodyEberlyCoordinates(rb2prop,rb2state)
rb2 = TRS.RigidBody(rb2prop,rb2state,rb2Eb)

s1 = TRS.LinearString(1.0,2.0,2,4,@MVector zeros(3))
j1 = TRS.BallJoint(1,3)
rbs = [rb1,rb2]
tgsys = TRS.TensegritySystem(rbs,[s1],[j1],nothing)

twobodystate = TRS.multibodystate(tgsys.rigidbodies)


function Aandb(rb,id,Fext,τext)
    @unpack mass,inertia,anchorpoints = rb.prop
    @unpack R,ω = rb.state
    L = rb.coords.L
    invJt = R*inv(inertia)*transpose(R)
    point = anchorpoints[id]
    x = R*point.p
    x̃ = tilde(x)
    A = I/mass - x̃*invJt*x̃
    b = Fext/mass + (invJt*(τext + L×ω)) × x +
        ω × (ω × x)
    A,b
end
function force2torque(force,rb,pointid)
    x = rb.prop.anchorpoints[pointid].p
    R = rb.state.R
    τ = (R*x × force)
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

function twobody(tgsys)
    function twobody!(du,u,p,t)
        rbs = tgsys.rigidbodies
        sts = tgsys.strings
        system_cstate_dot = du
        system_cstate = u
        nbody = length(rbs)
        nc = 13
        # State pass
        for i in eachindex(rbs)
            @unpack prop,state,coords = rbs[i]
            @unpack mass,inertia,anchorpoints = prop
            @unpack r,R,ṙ,ω = state
            @unpack x,q,p,L = coords
            rb_cstate = @view system_cstate[(i-1)*nc+1:(i-1)*nc+13]
            x .= rb_cstate[1:3]
            q .= rb_cstate[4:7]
            p .= rb_cstate[8:10]
            L .= rb_cstate[11:13]

            r .= x
            R .= SMatrix(Quat(q[1],q[2],q[3],q[4]))
            ṙ .= p/mass
            invJt = R*inv(inertia)*transpose(R)
            ω .= invJt*L

            for ip in eachindex(anchorpoints)
                state.p[ip] .= point_position(prop,ip,r,R)
            end
        end
        # External Force Pass
        for i in eachindex(rbs)
            rbstate = rbs[i].state
            for ip in eachindex(rbs[i].prop.anchorpoints)
                rbstate.Fext[ip] .= 0.0
                rbstate.τext[ip] .= 0.0
            end
            rbstate.F .= sum(rbstate.Fext)
            rbstate.τ .= sum(rbstate.τext)
        end
        #[0.0,0.0,-0.1]
        #force2torque(Fext2,rb2,2)
        # String Pass
        for is in eachindex(sts)
            st = sts[is]
            rb1 = rbs[1]
            rb2 = rbs[2]
            pid1 = 2
            pid2 = 2
            st.𝐬 .= rb2.state.p[pid2] - rb1.state.p[pid1]
            s = norm(s1.𝐬)
            Δs = s - s1.s0
            if Δs > 0.0
                fs = st.k*Δs*st.𝐬/s
            else
                fs = zeros(3)
            end

            rb1.state.Fext[pid1] .+= fs
            rb1.state.τext[pid1] .+= force2torque( fs,rb1,pid1)

            rb2.state.Fext[pid2] .+= -fs
            rb2.state.τext[pid2] .+= force2torque(-fs,rb2,pid2)
        end
        # Joints Pass
        for ij in eachindex(tgsys.joints)
            rb1 = rbs[1]
            rb2 = rbs[2]
            pid1 = 1
            pid2 = 1
            A1,b1 = Aandb(rb1,pid1,
                        sum(rb1.state.Fext),sum(rb1.state.τext))
            A2,b2 = Aandb(rb2,pid2,
                        sum(rb2.state.Fext),sum(rb2.state.τext))
            Fc = (A1+A2)\(b2-b1)
            rb1.state.Fext[pid1] .+=  Fc
            rb2.state.Fext[pid2] .+= -Fc
            rb1.state.τext[pid1] .+= force2torque( Fc,rb1,pid1)
            rb2.state.τext[pid2] .+= force2torque(-Fc,rb2,pid2)
        end
        #@show Fc,τc
        for i in eachindex(rbs)
            @unpack mass,inertia = rbs[i].prop
            @unpack r,R,ṙ,ω = rbs[i].state
            @unpack x,q,p,L = rbs[i].coords
            rb_cstate_dot = @view system_cstate_dot[(i-1)*nc+1:(i-1)*nc+13]
            ẋ = @view rb_cstate_dot[1:3]
            q̇ = @view rb_cstate_dot[4:7]
            ṗ = @view rb_cstate_dot[8:10]
            L̇ = @view rb_cstate_dot[11:13]
            ẋ .= ṙ
            q̇ .= 0.5*quat_multiply([0.0,ω[1],ω[2],ω[3]],q)

            ṗ .= sum(rbs[i].state.Fext)
            L̇ .= sum(rbs[i].state.τext)

        end

    end
end

fequ = twobody(tgsys)

sv,sdotv = TRS.multibodystate(rbs)
u0 = sv
tspan = (0.0,10.0)
prob = ODEProblem(fequ,u0,tspan)
sol = solve(prob,dt=0.01,adaptive = false)


function statevector(rbstate::TRS.RigidBodyEberlyCoordinates)
    state = Vector{Float64}(undef,13)
    state[1:3] = rbstate.x
    state[4:7] = rbstate.q
    state[8:10] = rbstate.p
    state[11:13] = rbstate.L
    statedot = Vector{Float64}(undef,13)
    statedot[1:3] = rbstate.ẋ
    statedot[4:7] = rbstate.q̇
    statedot[8:10] = rbstate.ṗ
    statedot[11:13] = rbstate.L̇
    state, statedot
end
rb1sv,rb1svdot = statevector(rb1Eb)

tmp_Eberly = deepcopy(rb1Eb)
tmp_state = deepcopy(rb1state)

singlebody!(rb1svdot,rb1sv,tmp_state,rb1)

@code_warntype singlebody!(rb1svdot,rb1sv,tmp_state,rb1)

function compute_energy(rb,state::TRS.RigidBodyState)
    @unpack mass, inertia = rb
    @unpack ṙ, ω = state
    translational = 0.5*mass*transpose(ṙ)*ṙ
    rotational = 0.5*transpose(ω)*inertia*ω
    total = translational + rotational
end
compute_energy(rb1,rb1state)

compute_energy(rb1,rb1sv)

using DifferentialEquations


function singlebody(tmp_state,rb)
    function f(du,u,p,t)
        singlebody!(du,u,tmp_state,rb)
    end
end
fequ = singlebody(tmp_state,rb1)

rb1sv,rb1svdot = statevector(rb1Eb)
u0 = rb1sv
tspan = (0.0,1.0)
prob = ODEProblem(fequ,u0,tspan)
sol = solve(prob)

[compute_energy(rb1,sol[:,i]) for i = 1:length(sol)]
function many_random_body(n)
    TRS = TensegrityRobotSimulator
    rbprops = Vector{TRS.RigidBodyProperty{Float64,2}}(undef,n)
    rbstates = Vector{TRS.RigidBodyState{Float64}}(undef,n)
    rbEbs = Vector{TRS.RigidBodyEberlyCoordinates{Float64}}(undef,n)
    for i = 1:n
        ap1 = TRS.AnchorPoint([1.0,0.0,0.0])
        ap2 = TRS.AnchorPoint(0.0,1.0,0.0)

        rbname = Symbol("rb"*string(i))
        mass = 1.0 + rand()
        inertia = SMatrix{3,3}(Symmetric(rand(3,3)))
        CoM = @SVector zeros(3)
        rbprop = TRS.RigidBodyProperty(rbname,:generic,
                    mass,
                    inertia,
                    CoM,
                    SVector(ap1,ap2))

        r = 5*rand(3)
        R = rand(RotMatrix{3})
        ṙ = 1*rand(3)
        ω = 5*rand(3)

        rbstate = TRS.RigidBodyState(r,R,ṙ,ω)
        rbEb = TRS.RigidBodyEberlyCoordinates(rbprop,rbstate)

        rbprops[i] = rbprop
        rbstates[i] = rbstate
        rbEbs[i] = rbEb
    end
    rbprops,rbstates,rbEbs
end

rbprops,rbstates,rbEbs  = many_random_body(10)
@code_warntype many_random_body(10)
rbs = [TRS.RigidBody(rbprops[i],rbstates[i],rbEbs[i]) for i = 1:10]

s1 = TRS.LinearString(:s1,1.0,1.0,1,2)
j1 = TRS.BallJoint(1,2)
tgsys = TRS.TensegritySystem(rbs,[s1],[j1],nothing)

compute_energy(rbs[1],sol[1:13,1])
function compute_energy(rb,state)
    @unpack mass, inertia = rb.prop
    p = @view state[8:10]
    L = @view state[11:13]
    ṙ = p./mass
    translational = 0.5*mass*transpose(ṙ)*ṙ
    q = @view state[4:7]
    R = Quat(q...)
    Jt = R*inertia*transpose(R)
    ω = inv(Jt)*L
    rotational = 0.5*transpose(ω)*Jt*ω
    total = translational + rotational
end
