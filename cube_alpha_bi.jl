using LinearAlgebra, Printf
using StaticArrays
using Rotations
using Parameters
using GeometryTypes, AbstractPlotting
using Makie
#using DifferentialEquations
using ForwardDiff
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
r = [0,0,10.5-0.0000001]
R = Matrix(1.0I,3,3)
#R = Matrix(RotXY(0.1,0.3))
#R = Matrix(RotX(float(π)))
ṙ = zeros(3)
#ω = zeros(3)
ω = [0.0,0.0,0.1]
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

h = 0.01
tspan = (0.0,0.1)
q0 = copy(cube.state.coords.q)
q̇0 = copy(cube.state.coords.q̇)

function generalized_α(ρ∞)
    am = (2ρ∞-1)/(ρ∞+1)
    af = ρ∞/(ρ∞+1)
    γ = 1/2 + af - am
    β = 1/4*(γ+1/2)^2
    am,af,γ,β
end
function Newmark(ρ∞)
    af = am = 0.0
    γ = 1/2
    β = 1/4
    am,af,γ,β
end

function alphastepping(h,tspan,q0,q̇0,cube)
    function scheme_parameters(ρ∞,h,f)
        am,af,γ,β = f(ρ∞)
        β̃ = (1-am)/((1-af)*h^2*β)
        γ̃ = γ/(h*β)
        am,af,γ,β,β̃,γ̃
    end
    αm,αf,γ,β,β̃,γ̃ = scheme_parameters(0.8,h,generalized_α)
    nq = 12
    ncon = 6
    @info "Starting new simulation"
    # simulation
    tstart = tspan[1]
    tend = tspan[2]
    ts = [tstart]
    qs = [copy(q0)]
    q̇s = [copy(q̇0)]

    # mass matrix
    @unpack M,Q = cube.state.auxs
    invM = inv(M)

    # External forces (Gravity)
    fG = [0,0,-9.8]
    QG = transpose(cube.state.auxs.CG)*fG
    #Q .= QG
    Q .= 0.0

    ṽ̇0 = invM*Q
    a0 = copy(ṽ̇0)
    ṽ̇s = [ṽ̇0]
    as = [a0]

    function Φq_jacobians(x,v,Λ)
        Φqq = reshape(ForwardDiff.jacobian(TRS.NC.Φq, x),ncon,nq,nq)
        retΦq₂ = zeros(ncon,nq)
        for i = 1:ncon
            retΦq₂[i,:] = transpose(v)*Φqq[i,:,:]
        end
        retAᵀΛq = zeros(nq,nq)
        for j = 1:nq
            retAᵀΛq[j,:] = transpose(Λ)*Φqq[:,j,:]
        end
        retΦq₂,retAᵀΛq
    end

    function compute_St(qₖ₊₁,vₖ₊₁,Λᵇₖ₊₁)
        Φq₂, AᵀΛq = Φq_jacobians(qₖ₊₁,vₖ₊₁,Λᵇₖ₊₁)
        Kₜ = AᵀΛq
        Φqₖ₊₁ = TRS.NC.Φq(qₖ₊₁)
        M = cube.state.auxs.M
        St = zeros(nq+ncon,nq+ncon)
        St[1:nq,        1:nq] .= M.*β̃ + Kₜ
        St[nq+1:nq+ncon,1:nq] .= Φqₖ₊₁

        St[1:nq,        nq+1:nq+ncon] .= transpose(Φqₖ₊₁)
        #St[nq+1:nq+ncon,nq+1:nq+ncon] .= 0.0

        St
    end

    DL = zeros(nq+ncon,nq+ncon)
    DR = zeros(nq+ncon,nq+ncon)
    for i = 1:nq
        DL[i,i] = β*h^2
        DR[i,i] = 1.0
    end
    for i = nq+1:nq+ncon
        DL[i,i] = 1.0
        DR[i,i] = 1/(β*h^2)
    end

    step = 0
    # time starts

    while ts[end] < tend
        qₖ₊₁ = copy(qs[1])
        vₖ₊₁ = copy(q̇s[1])
        ṽ̇ₖ₊₁ = copy(ṽ̇s[1])
        aₖ₊₁ = copy(as[1])
        qₖ = qs[end]
        vₖ = q̇s[end]
        ṽ̇ₖ = ṽ̇s[end]
        aₖ = as[end]
        #
        qₖ₊₁ .= qₖ + h*vₖ + h^2*(1/2-β)*aₖ
        vₖ₊₁ .= vₖ + h*(1-γ)*aₖ
       Λᵇₖ₊₁  = zeros(6)
        aₖ₊₁ .= 1/(1-αm)*(αf*ṽ̇ₖ-αm*aₖ)
        qₖ₊₁ .+= h^2*β*aₖ₊₁
        vₖ₊₁ .+= h*γ*aₖ₊₁
        ṽ̇ₖ₊₁ .= 0.0

        res = zeros(nq+ncon)
        tol = 1e-14
        iter_max = 1000
        iterations = 0
        res_total = 0.0

        for iLCP = 1:iter_max
            Φqₖ₊₁ = TRS.NC.Φq(qₖ₊₁)
            resq = M*ṽ̇ₖ₊₁+transpose(Φqₖ₊₁)*Λᵇₖ₊₁-Q
            resΦ = TRS.NC.Φ(qₖ₊₁)
            res_total = sqrt(norm(resq)^2 + norm(resΦ)^2)

            if res_total <= tol
                iterations = iLCP
                break
            elseif iLCP == iter_max
                @error "Reach max iteractions $iLCP, residual: $res_total"
                iterations = iter_max
            end
            St = compute_St(qₖ₊₁,vₖ₊₁,Λᵇₖ₊₁)
            S̄t = DL*St*DR
            invSt = inv(S̄t)

            res[1:nq] .= resq
            res[nq+1:nq+ncon] .= resΦ
            invStres = invSt*(DL*res)
            Δx = DR*(-invStres)
            Δq = Δx[1:nq]
            Δλ = Δx[nq+1:nq+ncon]
            qₖ₊₁ .+= Δq
            vₖ₊₁ .+= γ̃*Δq
            ṽ̇ₖ₊₁ .+= β̃*Δq
           Λᵇₖ₊₁ .+= Δλ
        end
        aₖ₊₁ .+= (1-αf)/(1-αm)*ṽ̇ₖ₊₁

        cube.state.coords.q .= qₖ₊₁
        cube.state.coords.q̇ .= vₖ₊₁
        TRS.coords2state_kinetic!(cube)
        push!(qs,qₖ₊₁)
        push!(q̇s,vₖ₊₁)
        push!(ṽ̇s,ṽ̇ₖ₊₁)
        push!(as,aₖ₊₁)
        step += 1
        push!(ts,ts[end] + h)
        @printf("Prog.: %5.1f%%, step: %s, time: %f, iterations: %s.\n", (
        ts[end]/tspan[end]*100), step, ts[end], iterations)
        ke = TRS.kinetic_energy(cube)
        constraint_error = norm(TRS.NC.Φ(qₖ₊₁))
        #@info "Current position: r = $(cube.state.r)"
        @info "KE: $(ke), Cons. Err.: $constraint_error"
    end

    ts,qs,q̇s
end
ts,qs,q̇s = alphastepping(0.01,(0.0,0.05),q0,q̇0,cube)
ts,qs,q̇s = onestepping(h,tspan,q0,q̇0,cube)
ts,qs,q̇s = onestepping(h,tspan,q289,q̇289,cube)
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
