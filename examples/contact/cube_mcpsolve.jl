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
R = Matrix(1.0I,3,3)
#R = Matrix(RotXY(0.1,0.3))
#R = Matrix(RotX(0.3))
#R = Matrix(RotX(float(π)))
ṙ = [0,0,-0.1]
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

function mcpsolve_alphastepping(h,tspan,q0,q̇0,cube)
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

    # constraints
    Λᵇ0 = zeros(6)
    Λᵇs= [Λᵇ0]

    step = 0
    # time starts

    while ts[end] < tend
        qₖ₊₁ = copy(qs[1])
        vₖ₊₁ = copy(q̇s[1])
        ṽ̇ₖ₊₁ = copy(ṽ̇s[1])
        aₖ₊₁ = copy(as[1])
        Λᵇₖ₊₁ = copy(Λᵇs[1])
        qₖ = qs[end]
        vₖ = q̇s[end]
        ṽ̇ₖ = ṽ̇s[end]
        aₖ = as[end]
        Λᵇₖ = Λᵇs[end]

        #
        ṽ̇ₖ₊₁ .= 0.0
        aₖ₊₁ .= 1/(1-αm)*(αf*ṽ̇ₖ-αm*aₖ)
        ṽₖ₊₁ = vₖ + h*(1-γ)*aₖ + h*γ*aₖ₊₁
        vₖ₊₁ .= ṽₖ₊₁

        # predict
        qᵖ = qₖ + 1/2*h*vₖ #qₖ₊₁
        active_index = [i for i = 1:8 if gaps[i](qᵖ) <= 0]
        @show active_index
        ν = length(active_index)
        CN = zeros(ν,12)
        ee = zeros(ν)
        for α = 1:ν
            C = cube.state.cache.Cp[active_index[α]]
            CN[α,:] .= C[3,:]
            ee[α] = e[active_index[α]]
        end
        gqₖ₊₁ = gqₖ = CN
        # unilateral
        Λᴬₖ₊₁= zeros(ν)

        function alpha!(ret,x)
            ṽ = @view x[1:nq]
            v = @view x[nq+1:2nq]
            Λᵇ = @view x[2nq+1:2nq+ncon]
            W = v - ṽ
            a = 1/(h*γ)*(ṽ - (vₖ + h*(1-γ)*aₖ))
            ṽ̇ = 1/(1-αf)*((1-αm)*a + αm*aₖ - αf*ṽ̇ₖ)
            q̃ = qₖ + h*vₖ + h^2*(1/2-β)*aₖ + h^2*β*a
            q = q̃ + 1/2*h*W
            Φq = TRS.NC.Φq(q)

            ret[1:nq] = M*ṽ̇ - Q
            ret[nq+1:2nq] = M*W + transpose(Φq)*Λᵇ
            ret[2nq+1:2nq+ncon] = Φq*v
        end

        function unialpha!(ret,x)
            ṽ = x[1:nq]
            v = x[nq+1:2nq]
            Λᵇ = x[2nq+1:2nq+ncon]
            Λᴬ = x[2nq+ncon+1:2nq+ncon+ν]
            W = v - ṽ
            a = 1/(h*γ)*(ṽ - (vₖ + h*(1-γ)*aₖ))
            ṽ̇ = 1/(1-αf)*((1-αm)*a + αm*aₖ - αf*ṽ̇ₖ)
            q̃ = qₖ + h*vₖ + h^2*(1/2-β)*aₖ + h^2*β*a
            q = q̃ + 1/2*h*W
            Φq = TRS.NC.Φq(q)

            ret[1:nq] = M*ṽ̇ - Q
            ret[nq+1:2nq] = M*W + transpose(Φq)*Λᵇ - transpose(gqₖ)*Λᴬ
            ret[2nq+1:2nq+ncon] = Φq*v
            ret[2nq+ncon+1:2nq+ncon+ν] = gqₖ₊₁*v + Diagonal(ee)*gqₖ*vₖ
        end

        function alpha(x)
            ret = zeros(2nq+ncon)
            alpha!(ret,x)
            ret
        end
        lb = fill(-Inf,2nq+ncon+ν)
        ub = fill( Inf,2nq+ncon+ν)
        if ν > 0
            lb[2nq+ncon+1:2nq+ncon+ν] .= 0
            initial_x = vcat(ṽₖ₊₁,vₖ₊₁,Λᵇₖ₊₁,Λᴬₖ₊₁)
            nlr = mcpsolve(unialpha!, lb, ub, initial_x, ftol=1e-14)

        else
            initial_x = vcat(ṽₖ₊₁,vₖ₊₁,Λᵇₖ₊₁)
            nlr = nlsolve(alpha!, initial_x, ftol=1e-14, autodiff = :forward)
        end


        iterations = nlr.iterations
        x = nlr.zero
        ṽₖ₊₁ .= x[1:nq]
        vₖ₊₁ .= x[nq+1:2nq]
        Λᵇₖ₊₁ .= x[2nq+1:2nq+ncon]
        #Λᴬₖ₊₁ .= x[2nq+ncon+1:2nq+ncon+ν]
        UNₖ₊₁ = CN*vₖ₊₁
        for α = 1:ν
            @info "Contact velocity: $α, $(UNₖ₊₁[α])"
        end

        aₖ₊₁ .= 1/(h*γ)*(ṽₖ₊₁ - (vₖ + h*(1-γ)*aₖ))
        ṽ̇ₖ₊₁ .= 1/(1-αf)*((1-αm)*aₖ₊₁ + αm*aₖ - αf*ṽ̇ₖ)
        q̃ₖ₊₁ = qₖ + h*vₖ + h^2*(1/2-β)*aₖ + h^2*β*aₖ₊₁
        Wₖ₊₁ = vₖ₊₁ - ṽₖ₊₁
        qₖ₊₁ .= q̃ₖ₊₁ + 1/2*h*Wₖ₊₁

        cube.state.coords.q .= qₖ₊₁
        cube.state.coords.q̇ .= vₖ₊₁
        TRS.coords2state_kinetic!(cube)
        push!(qs,qₖ₊₁)
        push!(q̇s,vₖ₊₁)
        push!(ṽ̇s,ṽ̇ₖ₊₁)
        push!(as,aₖ₊₁)
        push!(Λᵇs,Λᵇₖ₊₁)
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
ts,qs,q̇s = mcpsolve_alphastepping(0.01,(0.0,0.2),q0,q̇0,cube)

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
