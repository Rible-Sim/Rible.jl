using ForwardDiff
using Rotations
q = rand(Quat)
function rotation_matrix(q)
    R = Matrix{eltype(q)}(undef,3,3)
    q0,q1,q2,q3 = q
    R[1,1] = 2*q0^2 + 2*q1^2 - 1
    R[1,2] = 2*q1*q2 - 2q0*q3
    R[1,3] = 2*q1*q3 + 2q0*q2
    R[2,1] = 2*q1*q2 + 2q0*q3
    R[2,2] = 2*q0^2 + 2*q2^2 - 1
    R[2,3] = 2*q2*q3 - 2q0*q1
    R[3,1] = 2*q1*q3 - 2q0*q2
    R[3,2] = 2*q2*q3 + 2q0*q1
    R[3,3] = 2*q0^2 + 2*q3^2 - 1
    R
end
unitq = [q.w,q.x,q.y,q.z]
R = rotation_matrix(unitq)
ForwardDiff.jacobian(rotation_matrix,unitq)
function Q_matrix(q)
    Q = Matrix{eltype(q)}(undef,4,3)
    s,x,y,z = q
    Q[1,1] = -0.5x
    Q[1,2] = -0.5y
    Q[1,3] = -0.5z
    Q[2,1] =  0.5s
    Q[2,2] =  0.5z
    Q[2,3] = -0.5y
    Q[3,1] = -0.5z
    Q[3,2] =  0.5s
    Q[3,3] =  0.5x
    Q[4,1] =  0.5y
    Q[4,2] = -0.5x
    Q[4,3] =  0.5s
    Q
end
Q_matrix(unitq)
ForwardDiff.jacobian(Q_matrix,unitq)
w = rand(Quat)
quat_multiply([0.0,1.0,2.0,3.0],unitq)
2Q_matrix(unitq)*[1.0,2.0,3.0]

function S2_matrix(s)
    qi = s[4:7]
    qj = s[11:14]
    S = zeros(eltype(s),2*7,2*6)
    for i = 1:3
        S[i,i] = 1
    end
    S[4:7,4:6] = Q_matrix(qi)
    for i = 7:9
        S[i+1,i] = 1
    end
    S[11:14,10:12] = Q_matrix(qj)
    S
end
S2 = S2_matrix(s)
ForwardDiff.jacobian(S2_matrix,s)
function constraint(jo::BallJoint,rb1,rb2,pid1,pid2)
    p1 = rb1.prop.anchorpoints[pid1].p
    p2 = rb2.prop.anchorpoints[pid2].p
    function inner_Φ(s)
        ri = s[1:3]
        qi = s[4:7]
        rj = s[8:10]
        qj = s[11:14]
        ri + rotation_matrix(qi)*p1 - rj - rotation_matrix(qj)*p2
    end
end

function constraint_jacobian(jo,s,u,rb1,rb2,pid1,pid2)
    Φ = constraint(jo,rb1,rb2,pid1,pid2)
    Φs(s) = ForwardDiff.jacobian(Φ,s)
    function J(s)
        Φs(s)*S2_matrix(s)
    end
    nc = numberofconstraint(jo)
    Js(s) = ForwardDiff.jacobian(J,s)
    J(s),reshape(Js(s)*S2_matrix(s)*u,nc,:)
end

constraint_jacobian(jo,s,u,rb1,rb2,pid1,pid2)
function newton_euler_coords(rb)
    @unpack state,coords = rb
    @unpack ṙ,R,ω = state
    @unpack x,q = coords
    r = x
    v = ṙ
    r,q,v,transpose(R)*ω
end


s = vcat(rand(3),unitq,rand(3),unitq)
Φs = ForwardDiff.jacobian(Φ1,s)
function JΦ(s)
    Φs = ForwardDiff.jacobian(Φ1,s)
    S = S2_matrix(s)
    Φs*S
end
JΦ(s)
Js(s) = ForwardDiff.jacobian(JΦ,s)
Js(s)
*S2_matrix(s)
