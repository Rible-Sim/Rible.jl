using LinearAlgebra, Printf
using StaticArrays
using Rotations
using Parameters
using GeometryTypes, AbstractPlotting
using Makie
#using DifferentialEquations
#using ForwardDiff
using PATHSolver
using Suppressor
options(convergence_tolerance=1e-14, output=:no,
    time_limit=3600)
using Revise
using SPARK
using TensegrityRobotSimulator
const TRS = TensegrityRobotSimulator
mass = 1.0
inertia = Matrix(1.0I,3,3)
CoM = zeros(3)
r = [0,0,2.5-0.0000001]
#R = Matrix(1.0I,3,3)
#R = Matrix(RotXY(0.1,0.3))
R = Matrix(RotX(float(π)))
ṙ = zeros(3)
ω = zeros(3)
aplist = [[0.5,0.5,0.5], [0.5,-0.5,0.5],
          [-0.5,-0.5,0.5], [-0.5,0.5,0.5],
          [0.5,0.5,-0.5], [0.5,-0.5,-0.5],
          [-0.5,-0.5,-0.5], [-0.5,0.5,-0.5]]
apvector = SVector([TRS.AnchorPoint(aplist[i]) for i = 1:8]...)
prop = TRS.RigidBodyProperty(:cube,:generic,mass,
            SMatrix{3,3}(inertia),
            SVector(CoM...),
            apvector)

state = TRS.RigidBodyState(prop,r,R,ṙ,ω,Val(:NC))

cube = TRS.RigidBody(prop,state)

function mlcp_solution(A,B,C,D,a,b)
    invA = inv(A)
    DinvA = D*invA
    M = (B-DinvA*C)
    q = b-DinvA*a
    function lcpfunction(x)
        M*x+q
    end
    n = length(q)
    lb = zeros(n)
    ub = 1e20*ones(n)

    status, v = @suppress solveLCP(lcpfunction, M, lb, ub)
    if status != :Solved
        @error "LCP solution not found."
    else
        u = -invA*(C*v+a)
    end
    return u, v
end
h = 0.01
tspan = (0.0,1.0)
q0 = copy(cube.state.coords.q)
q̇0 = copy(cube.state.coords.q̇)
function onestepping(h,tspan,q0,q̇0,cube)
    @info "Starting new simulation"
    # collision related
    e = fill(1.0,8)
    function gap_function(cube,j)
        C = cube.state.auxs.Cp[j]
        function g(q)
            r = C*q
            ret = zero(r)
            ret[3] = r[3]
        end
    end
    gaps = [gap_function(cube,j) for j = 1:8]

    # simulation
    tstart = tspan[1]
    tend = tspan[2]
    ts = [tstart]
    qs = [copy(q0)]
    q̇s = [copy(q̇0)]
    step = 0
    # time starts
    while ts[end] < tend
        qₖ₊₁ = copy(qs[1])
        vₖ₊₁ = copy(q̇s[1])
        qₖ = qs[end]
        vₖ = q̇s[end]
        @unpack M,Q = cube.state.auxs
        invM = inv(M)
        # External forces (Gravity)
        fG = [0,0,-9.8]
        QG = transpose(cube.state.auxs.CG)*fG
        Q .= QG
        # "Free" velocity
        v_free = vₖ + invM*h*Q

        # predict
        q̃ = qₖ + 1/2*h*vₖ
        active_index = [i for i = 1:8 if gaps[i](q̃) < 0]
        @show active_index
        ν = length(active_index)
        # MLCP
        A = TRS.NC.Φq
        Ã = A(q̃)
        a = Ã*v_free
        AA = Ã*invM*transpose(Ã)
        if ν == 0
            #solve linear
            Pμ = AA\-a
            vₖ₊₁ .= v_free + invM*transpose(Ã)*Pμ
        else
            #solve MLCP
            CN = zeros(ν,12)
            UN = zeros(ν)
            ee = zeros(ν)
            for α = 1:ν
                C = cube.state.auxs.Cp[active_index[α]]
                CN[α,:] .= C[3,:]
                UN[α] = (C*vₖ)[3]
                @info "Pre-contact velocity: $α, $(UN[α])"
                ee[α] = e[active_index[α]]
            end
            C = Ã*invM*transpose(CN)
            b = CN*v_free + Diagonal(ee)*UN
            D = CN*invM*transpose(Ã)
            B = CN*invM*transpose(CN)
            Pμ, PN = mlcp_solution(AA,B,C,D,a,b)
            vₖ₊₁ .= v_free + invM*(transpose(Ã)*Pμ + transpose(CN)*PN)
            UNₖ₊₁ = CN*vₖ₊₁
            Û = CN*vₖ₊₁ + Diagonal(ee)*UN
            for α = 1:ν
                @info "Contact velocity and inpulse: $α, $(Û[α]), $(PN[α])"
            end
            mlcp_err = mlcp_compute_error(AA,B,C,D,a,b,Pμ,PN)
            @info "MLCP error: $mlcp_err"
            @info "$(inv(AA))"
        end
        qₖ₊₁ .= qₖ + 0.5*h*vₖ₊₁
        cube.state.coords.q .= qₖ₊₁
        cube.state.coords.q̇ .= vₖ₊₁
        TRS.coords2state_kinetic!(cube)
        push!(qs,qₖ₊₁)
        push!(q̇s,vₖ₊₁)
        step += 1
        push!(ts,ts[end] + h)
        @printf("Prog.: %5.1f%%, step: %s, time: %f, iterations: %s, Collsions: %s.\n",
                (ts[end]/tspan[end]*100), step, ts[end], iterations, ν)
        ke = TRS.kinetic_energy(cube)
        @info "Current position: r = $(cube.state.r)"
        @info "Current kinetic energy: ke = $(ke)"
    end

    ts,qs,q̇s
end
ts,qs,q̇s = onestepping(h,tspan,q0,q̇0,cube)

TRS.NC.Φ.(qs)
function energy(cube,q,q̇)
    cube.state.coords.q .= q
    cube.state.coords.q̇ .= q̇
    TRS.coords2state_kinetic!(cube)
    @unpack prop, state = cube
    @unpack mass = prop
    @unpack r = state
    pe = mass*9.8*r[3]
    ke = TRS.kinetic_energy(cube)
    ke,pe,pe + ke
end
ks,ps,es = [energy(cube,qs[it],q̇s[it]) for it in eachindex(ts)]
TRS.NC.Φ(q0)
# Pμ = zeros(6)
# PN = zeros(3ν)
# mlcp = MLCP(AA,B,C,D,a,b,length(a),length(b),Pμ,PN)
