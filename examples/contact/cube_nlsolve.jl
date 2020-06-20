using LinearAlgebra, Printf
using StaticArrays
using Rotations
using Parameters
using GeometryTypes, AbstractPlotting
using Makie
#using DifferentialEquations
using DiffEqDiffTools
using NLsolve
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
r = [0,0,5.5-0.0000001]
#R = Matrix(1.0I,3,3)
R = Matrix(RotXY(0.1,0.3))
#R = Matrix(RotX(0.3))
#R = Matrix(RotX(float(π)))
ṙ = zeros(3)
ω = [0.1,-0.2,0.3]
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

function lcp_solution(M,q)
    function lcpfunction(x)
        M*x+q
    end
    n = length(q)
    lb = zeros(n)
    ub = 1e20*ones(n)

    status, v = @suppress solveLCP(lcpfunction, M, lb, ub)
    if status != :Solved
        @error "LCP solution not found."
    end
    return v
end
h = 0.01
tspan = (0.0,0.1)
q0 = copy(cube.state.coords.q)
q̇0 = copy(cube.state.coords.q̇)

Ψfunc(a,b) = a - max(0.0,a-1.0b)

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

function unilateral_alphastepping(h,tspan,q0,q̇0,cube)
    function scheme_parameters(ρ∞,h,f)
        am,af,γ,β = f(ρ∞)
        γₜ = (1-am)/((1-af)*γ*h)
        βₜ = β*h/γ - 1/2*h
        am,af,γ,β,γₜ,βₜ
    end
    αm,αf,γ,β,γₜ,βₜ = scheme_parameters(0.8,h,generalized_α)
    nq = 12
    ncon = 6
    @info "Starting new simulation"
    # collision related
    e = fill(0.5,8)
    function gap_function(cube,j)
        C = cube.state.cache.Cp[j]
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

    # mass matrix
    @unpack M = cube.state.cache
    @unpack Q = cube.state.coords
    invM = inv(M)

    # External forces (Gravity)
    fG = [0,0,-9.8]
    QG = transpose(cube.state.cache.CG)*fG
    Q .= QG
    #Q .= 0.0

    ṽ̇0 = invM*Q
    a0 = copy(ṽ̇0)
    ṽ̇s = [ṽ̇0]
    as = [a0]

    function Φq₂(x,v)
        Φqq = reshape(ForwardDiff.jacobian(TRS.NC.Φq, x),ncon,nq,nq)
        ret = zeros(ncon,nq)
        for i = 1:ncon
            ret[i,:] = transpose(v)*Φqq[i,:,:]
        end
        ret
    end
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
        γq = (1-αm)/((1-αf)*h^2*β)
        # Cₜ = [0.0]
        Φq₂, AᵀΛq = Φq_jacobians(qₖ₊₁,vₖ₊₁,Λᵇₖ₊₁)
        Kₜ = γq*M
        Kₜᶜᵒⁿ = 2/h*M + AᵀΛq
        Φqₖ₊₁ = TRS.NC.Φq(qₖ₊₁)

        St = zeros(2nq+ncon,2nq+ncon)
        St[1:nq,          1:nq] .= γₜ*M + βₜ*Kₜ
        St[nq+1:2nq,      1:nq] .=   -M + βₜ*Kₜᶜᵒⁿ
        St[2nq+1:2nq+ncon,1:nq] .=        βₜ*Φq₂

        St[1:nq,          nq+1:2nq] .=         1/2*h*Kₜ
        St[nq+1:2nq,      nq+1:2nq] .=     M + 1/2*h*Kₜᶜᵒⁿ
        St[2nq+1:2nq+ncon,nq+1:2nq] .= Φqₖ₊₁ + 1/2*h*Φq₂

        #St[1:nq,          2nq+1:2nq+ncon] .= 0.0
        St[nq+1:2nq,      2nq+1:2nq+ncon] .= transpose(Φqₖ₊₁)
        #St[2nq+1:2nq+ncon,2nq+1:2nq+ncon] .= 0.0
        St
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
        ṽ̇ₖ₊₁ .= 0.0
        aₖ₊₁ .= 1/(1-αm)*(αf*ṽ̇ₖ-αm*aₖ)
        ṽₖ₊₁ = vₖ + h*(1-γ)*aₖ + h*γ*aₖ₊₁
        vₖ₊₁ .= ṽₖ₊₁
        Λᵇₖ₊₁= zeros(6)

        function alpha!(ret,x)
            ṽ = x[1:nq]
            v = x[nq+1:2nq]
            Λᵇ = x[2nq+1:2nq+ncon]
            W = v - ṽ
            a = 1/(h*γ)*(ṽ - (vₖ + h*(1-γ)*aₖ))
            ṽ̇ = 1/(1-αf)*((1-αm)*a + αm*aₖ - αf*ṽ̇ₖ)
            q = qₖ + h*vₖ + h^2*(1/2-β)*aₖ + h^2*β*a
            Φq = TRS.NC.Φq(q)

            ret[1:nq] = M*ṽ̇ - Q
            ret[nq+1:2nq] = M*W + transpose(Φq)*Λᵇ
            ret[2nq+1:2nq+ncon] = Φq*v
        end
        function alpha(x)
            ret = zeros(2nq+ncon)
            alpha!(ret,x)
            ret
        end
        initial_x = vcat(ṽₖ₊₁,vₖ₊₁,Λᵇₖ₊₁)
        St = compute_St(qₖ₊₁,vₖ₊₁,Λᵇₖ₊₁)
        alpha_jacobian = DiffEqDiffTools.finite_difference_jacobian(alpha,initial_x)
        display(alpha_jacobian[nq+1:2nq,      2nq+1:2nq+ncon])
        display("St")
        display(St[nq+1:2nq,      2nq+1:2nq+ncon])
        nlr = nlsolve(alpha!, initial_x, ftol=1e-14)
        iterations = nlr.iterations
        x = nlr.zero
        ṽₖ₊₁ = x[1:nq]
        vₖ₊₁ = x[nq+1:2nq]
        Λᵇₖ₊₁ = x[2nq+1:2nq+ncon]

        cube.state.coords.q .= qₖ₊₁
        cube.state.coords.q̇ .= vₖ₊₁
        TRS.coords2state_kinetic!(cube)
        push!(qs,qₖ₊₁)
        push!(q̇s,vₖ₊₁)
        push!(ṽ̇s,ṽ̇ₖ₊₁)
        push!(as,aₖ₊₁)
        step += 1
        push!(ts,ts[end] + h)
        @printf("Prog.: %5.1f%%, step: %s, time: %f, iterations: %s.\n",
                (ts[end]/tspan[end]*100), step, ts[end], iterations)
        ke = TRS.kinetic_energy(cube)
        constraint_error = norm(TRS.NC.Φ(qₖ₊₁))
        #@info "Current position: r = $(cube.state.r)"
        @info "KE: $(ke), Cons. Err.: $constraint_error"
    end

    ts,qs,q̇s
end
ts,qs,q̇s = unilateral_alphastepping(0.01,(0.0,0.01),q0,q̇0,cube)

ts,qs,q̇s = alphastepping(0.001,(0.0,1.0),q0,q̇0,cube)
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
