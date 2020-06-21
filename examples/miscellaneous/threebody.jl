using LinearAlgebra
using StaticArrays
using Rotations
using Parameters
using DifferentialEquations
using ForwardDiff
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

r1 = [0.0,0.0,-1.0]
R1 = RotX(1.0π)
ṙ1 = [0.0,0.1,0.0]
ω1 = [0.1,0.0,0.0]

rb1state = TRS.RigidBodyState(rb1prop,r1,R1,ṙ1,ω1)
rb1 = TRS.RigidBody(rb1prop,rb1state)

rb2prop = TRS.RigidBodyProperty(:rb2,:generic,1.0,
            SMatrix{3,3}(1.0I),
            CoM1,
            SVector(ap1,ap2))

r2 = [0.0,√3/2,1/2]
R2 = RotX(-π/3)
ṙ2 = [0.0,0.0,0.0]
ω2 = [0.0,0.0,0.0]

rb2state = TRS.RigidBodyState(rb2prop,r2,R2,ṙ2,ω2)
rb2 = TRS.RigidBody(rb2prop,rb2state)

rb3prop = TRS.RigidBodyProperty(:rb3,:generic,1.0,
            SMatrix{3,3}(1.0I),
            CoM1,
            SVector(ap1,ap2))
θ3 = π/3
r3 = [0.0,-sin(θ3),cos(θ3)]
R3 = RotX(θ3)
ṙ3 = [0.0,0.0,0.0]
ω3 = [0.0,0.0,0.0]

rb3state = TRS.RigidBodyState(rb3prop,r3,R3,ṙ3,ω3)
rb3 = TRS.RigidBody(rb3prop,rb3state)

s1 = TRS.LinearString(10.0,1.0,2,4,@MVector zeros(3))
s2 = TRS.LinearString(10.0,1.0,4,6,@MVector zeros(3))
s3 = TRS.LinearString(10.0,1.0,6,2,@MVector zeros(3))

j1 = TRS.BallJoint(1,3)
j2 = TRS.BallJoint(3,5)

rbs = [rb1,rb2,rb3]
connectivity = [
    [1,1],
    [1,2],
    [2,1],
    [2,2],
    [3,1],
    [3,2]]
#tgsys = TRS.TensegritySystem(rbs,[s1,s2,s3],[j1,j2],connectivity)
tgsys = TRS.TensegritySystem([rb1,rb2],[s1],[j1],connectivity)

twobodystate = TRS.multibodystate(tgsys.rigidbodies)


function Aandb(rb,id,Faps,τaps)
    @unpack mass,inertia,aps = rb.prop
    @unpack R,ω = rb.state
    L = rb.state.coords.L
    invJt = R*inv(inertia)*transpose(R)
    point = aps[id]
    x = R*point.p
    x̃ = tilde(x)
    A = I/mass - x̃*invJt*x̃
    b = Faps/mass + (invJt*(τaps + L×ω)) × x +
        ω × (ω × x)
    A,b
end
function force2torque(force,rb,pointid)
    x = rb.prop.aps[pointid].p
    R = rb.state.R
    τ = R*x × force
end

function constant_mass_matrix(rbs)
    nb = length(rbs)
    M = zeros(eltype(rbs).parameters[1],6nb,6nb)
    for i = 1:nb
        is = 6*(i-1)
        for j = 1:3
            M[is+j,is+j] = rbs[i].prop.mass
        end
        M[is+4:is+6,is+4:is+6] .= rbs[i].prop.inertia
    end
    M
end
function init_jacobian_matrix(tgsys)
    @unpack joints = tgsys
    nj = length(joints)
    nc_list = zeros(Int,nj)
    c_position = zeros(Int,nj)
    nq = 6
    totalnc = 0
    for ij in eachindex(joints)
        nc = numberofconstraint(joints[ij])
        nc_list[ij] = nc
        c_position[ij] = sum(nc_list[1:ij])-nc
        totalnc += nc
    end
    nb = length(tgsys.rigidbodies)
    J = zeros(eltype(rbs).parameters[1],totalnc,nq*nb)
    Qd = zeros(eltype(rbs).parameters[1],totalnc)
    J,c_position,Qd
end
function threebody(tgsys)
    function threebody!(du,u,p,t)
        rbs = tgsys.rigidbodies
        sts = tgsys.strings
        cnt = tgsys.connectivity
        system_cstate_dot = du
        system_cstate = u
        nbody = length(rbs)
        ncstate = 13
        # State pass
        for i in eachindex(rbs)
            @unpack prop,state = rbs[i]
            @unpack mass,inertia,aps = prop
            @unpack r,R,ṙ,ω,coords = state
            @unpack x,q,p,L = coords
            rb_cstate = @view system_cstate[(i-1)*ncstate+1:(i-1)*ncstate+13]
            x .= rb_cstate[1:3]
            q .= rb_cstate[4:7]
            p .= rb_cstate[8:10]
            L .= rb_cstate[11:13]

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
                rbstate.Faps[ip] .= 0.0
                rbstate.τaps[ip] .= 0.0
            end
            rbstate.F .= 0.0
            rbstate.τ .= 0.0
        end
        #[0.0,0.0,-0.1]
        #force2torque(Faps2,rb2,2)
        # String Pass
        for is in eachindex(sts)
            st = sts[is]
            rbid1,pid1 = cnt[st.pid1]
            rbid2,pid2 = cnt[st.pid2]
            rb1 = rbs[rbid1]
            rb2 = rbs[rbid2]
            st.𝐬 .= rb2.state.p[pid2] - rb1.state.p[pid1]
            s = norm(st.𝐬)
            Δs = s - st.s0
            if Δs > 0.0
                fs = st.k*Δs*st.𝐬/s
            else
                fs = zeros(eltype(Δs),3)
            end

            rb1.state.Faps[pid1] .+= fs
            rb1.state.τaps[pid1] .+= force2torque( fs,rb1,pid1)

            rb2.state.Faps[pid2] .+= -fs
            rb2.state.τaps[pid2] .+= force2torque(-fs,rb2,pid2)
        end
        for i in eachindex(rbs)
            rbs[i].state.F .= sum(rbs[i].state.Faps)
            rbs[i].state.τ .= sum(rbs[i].state.τaps)
        end
        # Joints Pass
        M = constant_mass_matrix(rbs)
        invM = inv(M)
        J,c_position,Qd = init_jacobian_matrix(tgsys)
        for ij in eachindex(tgsys.joints)
            jo = tgsys.joints[ij]
            rbid1,pid1 = cnt[jo.pid1]
            rbid2,pid2 = cnt[jo.pid2]
            rb1 = rbs[rbid1]
            rb2 = rbs[rbid2]
            r1,q1,v1,ω1 = newton_euler_coords(rb1)
            r2,q2,v2,ω2 = newton_euler_coords(rb2)
            s = vcat(r1,q1,r2,q2)
            u = vcat(v1,ω1,v2,ω2)
            Ji,J̇i = constraint_jacobian(jo,s,u,rb1,rb2,pid1,pid2)
            nc = numberofconstraint(jo)
            js = c_position[ij]
            J[js+1:js+nc,6(rbid1-1)+1:6(rbid1-1)+6] .= Ji[:,1:6]
            J[js+1:js+nc,6(rbid2-1)+1:6(rbid2-1)+6] .= Ji[:,7:12]
            Qd[js+1:js+nc] = -J̇i*u
        end
        Qe = zeros(eltype(rbs).parameters[1],6length(rbs))
        for i in eachindex(rbs)
            is = 6(i-1)
            R = rbs[i].state.R
            Qe[is+1:is+3] = rbs[i].state.F
            Qe[is+4:is+6] = transpose(R)*rbs[i].state.τ
        end
        H = J*invM*transpose(J)
        λ = H\(J*invM*Qe-Qd)
        Qc = transpose(J)*λ
        for i in eachindex(rbs)
            is = 6(i-1)
            R = rbs[i].state.R
            Qc[is+4:is+6] = R*Qc[is+4:is+6]
            rbs[i].state.F .+= -Qc[is+1:is+3]
            rbs[i].state.τ .+= -Qc[is+4:is+6]
        end
        # @show -Qc
        # for ij in eachindex(tgsys.joints)
        #     jo = tgsys.joints[ij]
        #     rbid1,pid1 = cnt[jo.pid1]
        #     rbid2,pid2 = cnt[jo.pid2]
        #     rb1 = rbs[rbid1]
        #     rb2 = rbs[rbid2]
        #     A1,b1 = Aandb(rb1,pid1,
        #                 sum(rb1.state.Faps),sum(rb1.state.τaps))
        #     A2,b2 = Aandb(rb2,pid2,
        #                 sum(rb2.state.Faps),sum(rb2.state.τaps))
        #     Fc = (A1+A2)\(b2-b1)
        #     rb1.state.Faps[pid1] .+=  Fc
        #     rb2.state.Faps[pid2] .+= -Fc
        #     rb1.state.τaps[pid1] .+= force2torque( Fc,rb1,pid1)
        #     rb2.state.τaps[pid2] .+= force2torque(-Fc,rb2,pid2)
        #     FF = vcat(Fc,force2torque( Fc,rb1,pid1),-Fc,force2torque(-Fc,rb2,pid2))
        #     @show FF
        # end
        # for i in eachindex(rbs)
        #     rbs[i].state.F .= sum(rbs[i].state.Faps)
        #     rbs[i].state.τ .= sum(rbs[i].state.τaps)
        # end
        for i in eachindex(rbs)
            @unpack mass,inertia = rbs[i].prop
            @unpack r,R,ṙ,ω,coords = rbs[i].state
            @unpack x,q,p,L = coords
            rb_cstate_dot = @view system_cstate_dot[(i-1)*ncstate+1:(i-1)*ncstate+13]
            ẋ = @view rb_cstate_dot[1:3]
            q̇ = @view rb_cstate_dot[4:7]
            ṗ = @view rb_cstate_dot[8:10]
            L̇ = @view rb_cstate_dot[11:13]
            ẋ .= ṙ
            q̇ .= 0.5*quat_multiply([0.0,ω[1],ω[2],ω[3]],q)

            ṗ .= rbs[i].state.F
            L̇ .= rbs[i].state.τ

        end

    end
end

fequ = threebody(tgsys)

sv,sdotv = TRS.multibodystate(rbs)
u0 = sv
tspan = (0.0,10.0)
prob = ODEProblem(fequ,u0,tspan)
sol = solve(prob,dt=0.01,adaptive = false)
function kinetic_energy(rbs)
    translational = sum([0.5*rbs[i].prop.mass*transpose(rbs[i].state.ṙ)*rbs[i].state.ṙ for i in eachindex(rbs)])
    rotational = zero(eltype(rbs).parameters[1])
    for i in eachindex(rbs)
        @unpack inertia = rbs[i].prop
        @unpack R,ω = rbs[i].state
        local_ω = transpose(R)*ω
        rotational += 0.5*transpose(local_ω)*inertia*local_ω
    end
    translational + rotational
end
function potential_energy(sts)
    potential = zero(eltype(sts).parameters[1])
    for i in eachindex(sts)
        @unpack k,s0,𝐬 = sts[i]
        s = norm(𝐬)
        Δs = s - s0
        if Δs > 0.0
            potential += 0.5*k*Δs^2
        end
    end
    potential
end
function energy(tgsys)
    kinetic_energy(tgsys.rigidbodies)+potential_energy(tgsys.strings)
end

E = [begin sol2state!(tgsys,sol[:,i]);energy(tgsys) end for i = 1:length(sol)]
