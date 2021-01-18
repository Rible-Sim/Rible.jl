using LinearAlgebra
using Parameters
using StaticArrays
using Cthulhu
import PyPlot; const plt = PyPlot
using Plots
using DifferentialEquations
using Revise
import TensegrityRobots; const TR = TensegrityRobots
const NC = TR.NaturalCoordinates

force(t) = 0.1*[sin(t),cos(t)]
torque(t) = 0.1cos(t)

function freefall()
    m = 0.9
    r̄g = [0.1,0.2]
    a = 0.5
    ami_g = SMatrix{2,2}([1/12*a^4  -1/24*a^4;-1/24*a^4 1/12*a^4])
    # ami_g = SMatrix{2,2}([1/12*a^4  0.0;0.0 1/12*a^4])
    θ = 0.1
    ω = 0.3
    polar_g = sum(diag(ami_g))
    ro = [0.2,0.3]
    ṙo = [0.1,0.4]
    rg = ro + TR.rotation_matrix(θ)*r̄g
    ṙg = ṙo + ω×(rg-ro)
    T_polar_g = 1/2*m*ṙg⋅ṙg + 1/2*ω*polar_g*ω
    @show rg,ṙg
    ri = @SVector zeros(2)
    rj = @SVector [0.6,0.4]
    rk = @SVector [-0.4,0.2]
    r = rj - ri
    r = rk - ri

    # lncs = NC.LocalNaturalCoordinates1P2V(r̄i,ū,v̄)
    # lncs = NC.LocalNaturalCoordinates2P1V(r̄i,r̄j,v̄)
    # lncs = NC.LocalNaturalCoordinates3P(r̄i,r̄j,r̄k)

    lncs,q0,q̇0 = NC.NC1P2V(ri,ro,θ,ṙo,ω)
    # lncs,q0,q̇0 = NC.NC2P1V(ri,rj,ro,θ,ṙo,ω)
    # lncs,q0,q̇0 = NC.NC3P(ri,rj,rk,ro,θ,ṙo,ω)
    cf = NC.CoordinateFunctions(lncs)

    M = NC.make_M(cf,m,ami_g,r̄g)
    T_nc = 1/2*transpose(q̇0)*M*q̇0
    @show T_polar_g,T_nc
    @show rg,cf.C(cf.c(r̄g))*q0
    @show ṙg,cf.C(cf.c(r̄g))*q̇0
    prop = TR.RigidBodyProperty(1,true,m,ami_g,r̄g,[r̄g])
    state = TR.RigidBodyState(prop,lncs,ro,θ,ṙo,ω,q0,q̇0)
    rb = TR.RigidBody(prop,state)
    # @show cf.Φ(q0)
    # lncs,q,q̇ = NC.NC2P1V(ri,rj,ro,θ,ṙo,ω)
    # lncs,q,q̇ = NC.NC3P(ri,rj,rk,ro,θ,ṙo,ω)
    λ0 = zeros(eltype(q0),3)
    function dynfuncs()
        @unpack CoM = rb.prop
        @unpack C,c = rb.state.cache.funcs
        Cg = rb.state.cache.CG
        ū1 = [1.0,0.0]
        C1 = C(c(CoM.+ū1))
        M = Array(NC.make_M(cf,m,ami_g,r̄g))
        Φ = cf.Φ
        A = cf.Φq
        function F!(F,q,q̇,t)
            F .= zero(q)
            F .+= transpose(Cg)*force(t)
            rg = Cg*q
            r1 = C1*q
            f_torque = torque(t)×(r1-rg)
            F .+= transpose(C1)*f_torque
            F .-= transpose(Cg)*f_torque
        end
        M,Φ,A,F!,nothing
    end
    prob = TR.DyProblem(dynfuncs(),Array(q0),Array(q̇0),Array(λ0),(0.0,100.0))
    dt = 0.01
    sol = TR.solve(prob,TR.Zhong06(),dt=dt,ftol=1e-14)
    rb, sol
end
rb1, sol1 = freefall()
function get_state(rb,q,q̇)
    @unpack C,c = rb.state.cache.funcs
    r̄g = rb.prop.CoM
    Cg = C(c(r̄g))
    rg = Cg*q
    ṙg = Cg*q̇
    Co = C(c([0.0,0.0]))
    ro = Co*q
    ṙo = Co*q̇
    C1 = C(c([1.0,0.0]))
    r1 = C1*q
    ṙ1 = C1*q̇
    u1 = r1 - ro
    θ = atan(u1[2],u1[1])
    θ = mod2pi(θ+2π)
    ṙ1o = ṙ1-ṙo
    # ω×u1 = [-ω*u1[2], ω*u1[1]]
    # @assert ṙ1o[1]/(-u1[2])≈ṙ1o[2]/(u1[1]) "$(ṙ1o[1]/(-u1[2]))≠$(ṙ1o[2]/(u1[1]))"
    ω = ṙ1o[2]/(u1[1])
    rg,θ,ṙg,ω
end
function plot_trajectory(rb,sol)
    @unpack ts,qs,q̇s = sol
    T = TR.get_numbertype(rb)
    xs = Vector{T}()
    ys = Vector{T}()
    θs = Vector{T}()
    ẋs = Vector{T}()
    ẏs = Vector{T}()
    ωs = Vector{T}()
    for (q,q̇) in zip(qs,q̇s)
        rg,θ,ṙg,ω = get_state(rb,q,q̇)
        push!(xs,rg[1])
        push!(ys,rg[2])
        push!(θs,θ)
        push!(ẋs,ṙg[1])
        push!(ẏs,ṙg[2])
        push!(ωs,ω)
    end
    ts = sol.ts
    # fig,axs = plt.subplots(1,3)
    # axs[1].plot(ts,xs,label="x")
    # axs[1].plot(ts,ys,label="y")
    # axs[2].plot(ts,θs,label="θ")
    # axs[3].plot(ts,ωs,label="ω")
    # axs[3].plot(ts,ẋs,label="ẋ")
    # axs[3].plot(ts,ẏs,label="ẏ")
    # fig.legend()
    fig,ax = plt.subplots(1,1)
    ax.plot(ts,θs,label="θ")
    fig
end
fig = plot_trajectory(rb1,sol1)
fig
function odesolve(rb,force,torque)
    m = rb.prop.mass
    Ig = sum(diag(rb.prop.inertia))
    function rbd!(ddu,du,u,p,t)
        ddu[1:2] .= 0.0./m
        ddu[3] = torque(t)/Ig
    end
    θ = 0.1
    ω = rb.state.ω
    u0 = vcat(rb.state.r,θ)
    du0 = vcat(rb.state.ṙ,ω)
    prob = SecondOrderODEProblem(rbd!,du0,u0,(0.0,100.0))
    sol = solve(prob,DPRKN6())
end
sol_ode = odesolve(rb1,force,torque)
pyplot()
plot(sol_ode,vars = [1, 2, 3])
plot(sol_ode,vars = [4,5])
plot(sol_ode,vars = [((t,θ)->(t,mod2pi(θ)),0,6)])
