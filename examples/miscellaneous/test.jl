using LinearAlgebra
using StaticArrays
using Rotations
using Parameters
using Revise
using TensegrityRobotSimulator
TRS = TensegrityRobotSimulator

ap1 = TRS.AnchorPoint([1.0,0.0,0.0])
ap2 = TRS.AnchorPoint(0.0,1.0,0.0)
@code_warntype TRS.AnchorPoint(0.0,1.0,0.0)

CoM1 = @SVector zeros(3)
rb1 = TRS.RigidBody(:rb1,:generic,1.0,
            SMatrix{3,3}(1.0I),
            CoM1,
            SVector(ap1,ap2))

r1 = [1.0,1.0,2.0]
R1 = Matrix(1.0I,3,3)
ṙ1 = [0.1,0.0,0.2]
ω1 = [1.0,2.0,3.0]

rb1state = TRS.RigidBodyState(r1,R1,ṙ1,ω1)
rb1Eb = TRS.RigidBodyEberlyCoordinates(rb1,rb1state)
rb1Er = TRS.RigidBodyErlebenCoordinates(rb1,rb1state)

rb2 = TRS.RigidBody(:rb2,:generic,2.0,
            SMatrix{3,3}(1.0I),
            CoM1,
            SVector(ap1,ap2))
rb2state = deepcopy(rb1state)
rb2Eb = TRS.RigidBodyEberlyCoordinates(rb2,rb2state)
tgsys = TRS.TensegritySystem([rb1,rb2],[rb1state,rb2state],[rb1Eb,rb2Eb],[s1],[j1],nothing)


function singlebody!(statedot,state,tmp_state,rb)
    @unpack ω,R = tmp_state
    @unpack mass,inertia = rb

    x = @view state[1:3]
    q = @view state[4:7]
    p = @view state[8:10]
    L = @view state[11:13]
    R .= SMatrix(Quat(q[1],q[2],q[3],q[4]))
    ω .= R*inv(inertia)*transpose(R)*L

    ẋ = @view statedot[1:3]
    q̇ = @view statedot[4:7]
    ṗ = @view statedot[8:10]
    L̇ = @view statedot[11:13]
    ẋ .= inv(mass)*p
    q̇ .= 0.5*quat_multiply([0.0,ω[1],ω[2],ω[3]],q)
    ṗ .= 0.0#F(t)
    L̇ .= 0.0#τ(t)
end



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


function multibody(rbs)
    function multbody!(du,u,p,t)
        statedotvector = du
        statevector = u
        nbody = length(rbs)
        nc = 13
        for i in eachindex(rbs)
            @unpack mass,inertia = rbs[i].prop
            state    = @view    statevector[(i-1)*nc+1:(i-1)*nc+13]
            statedot = @view statedotvector[(i-1)*nc+1:(i-1)*nc+13]

            x = @view state[1:3]
            q = @view state[4:7]
            p = @view state[8:10]
            L = @view state[11:13]
            R = SMatrix(Quat(q[1],q[2],q[3],q[4]))
            ω = R*inv(inertia)*transpose(R)*L

            ẋ = @view statedot[1:3]
            q̇ = @view statedot[4:7]
            ṗ = @view statedot[8:10]
            L̇ = @view statedot[11:13]
            ẋ .= inv(mass)*p
            q̇ .= 0.5*quat_multiply([0.0,ω[1],ω[2],ω[3]],q)
            ṗ .= 0.0#F(t)
            L̇ .= 0.0#τ(t)
        end
    end
end

fequ = multibody(rbs)


sv,sdotv = multibodystate(rbs)
u0 = sv
tspan = (0.0,1.0)
prob = ODEProblem(fequ,u0,tspan)
sol = solve(prob)

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

function coords2state(rb)
    @unpack mass,inertia = rb.prop
    @unpack r,ṙ,R,ω = rb.state
    @unpack x,q,p,L = rb.coords
    r .= x
    ṙ .= p/mass
    R .= SMatrix(Quat(q[1],q[2],q[3],q[4]))
    Jt = R*inertia*transpose(R)
    ω .= inv(Jt)*L
end
