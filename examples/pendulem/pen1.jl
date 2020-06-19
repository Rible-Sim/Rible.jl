using LinearAlgebra
using BenchmarkTools
using StaticArrays
using ForwardDiff
using NLsolve
using Revise
using SPARK
using Robot2D
const R2 = Robot2D

function penfunc()
    m = 1.0
    L = 1.0
    prop = R2.RigidBody2DProperty(true,:rb1,:generic,
                m,m*L,
                SVector(L,0.0),
                [SVector(0.0,0.0),SVector(L,0.0)]
                )
    θ = -π/3
    state = R2.RigidBody2DState(prop,[0.0,0.0],[cos(θ),sin(θ)])
    rb = R2.RigidBody2D(prop,state)
    const_mass_matrix = @SMatrix [m/3   0 m/6   0;
                   0 m/3   0 m/6;
                 m/6   0 m/3   0;
                   0 m/6   0 m/3]
    #M = const_mass_matrix
    M = rb.state.auxs.M
    # CG = [1/2   0 1/2   0;
    #        0 1/2   0 1/2;]
    CG = rb.state.auxs.CG
    g = 9.8
    fg = [0, -g]
    const_F = SVector{4}(transpose(CG)*fg)
    function F(q,q̇,t)
       const_F
    end


    function Φ(q)
      xi,yi,xj,yj = q
      [xi,
       yi,
       (xj-xi)^2 + (yj-yi)^2 - L^2]
    end

    function A(q)
      xi,yi,xj,yj = q
      ret = zeros(eltype(q),3,4)
      ret[1,1] = 1
      ret[2,2] = 1
      ret[3,1] = -2(xj-xi)
      ret[3,2] = -2(yj-yi)
      ret[3,3] = 2(xj-xi)
      ret[3,4] = 2(yj-yi)
      ret
    end
    M,Φ,A,F,nothing
end
M,Φ,A,F,Jacs = penfunc()
θ = -π/3
q0 = [0,0,cos(θ),sin(θ)]
q̇0 = zeros(4)
λ0 = zeros(3)
nq = length(q0)
nλ = length(λ0)
dt = 0.01
function R1!(R,x)
    h = dt
    myq1 = @view x[   1:nq]
    myλ0 = @view x[nq+1:nq+nλ]
    R[   1:nq]    .= M*((myq1-q0)/h - q̇0)/(h/2) -
                     F(1,2,1) -
                     transpose(A(q0))*myλ0
    R[nq+1:nq+nλ] .= Φ(myq1)
end

initial_x = vcat(q0,λ0)
initial_R = similar(initial_x)
R1_result = nlsolve(R1!,initial_x,ftol=1e-13)
@show R1_result.zero
q1 = R1_result.zero[1:nq]
λ0 = R1_result.zero[nq+1:nq+nλ]

qs = [copy(q0),copy(q1)]
λs = [copy(λ0)]

for timestep = 1:400
    qᵏ = qs[end]
    qᵏ⁻¹ = qs[end-1]
    function R2!(R,x)
        h = dt
        qᵏ⁺¹ = @view x[   1:nq]
        λᵏ   = @view x[nq+1:nq+nλ]
        R[   1:nq]    .= M*(qᵏ⁺¹.-2 .*qᵏ.+qᵏ⁻¹)./h^2 .-
                         F(1,2,1) .-
                         transpose(A(qᵏ))*λᵏ
        R[nq+1:nq+nλ] .= Φ(qᵏ⁺¹)
    end

    initial_x = vcat(qᵏ,λs[end])
    R2_result = nlsolve(R2!,initial_x,ftol=1e-14, autodiff = :forward)
    push!(qs,R2_result.zero[1:nq])
    @show R2_result.zero[3:4]
    push!(λs,R2_result.zero[nq+1:nq+nλ])
    @show R2_result.iterations
end
using Plots
plot([q[3] for q in qs])
plot!([q[4] for q in qs])
prob = SPARK.DyProblem(penfunc(),q0,q̇0,λ0,(0.0,0.02))

state = SPARK.solve(prob,dt=dt,ftol=1e-12)
@code_warntype SPARK.solve(prob,dt=dt,ftol=1e-12)
using DifferentialEquations
using Sundials
function f(out,du,u,p,t)
    q = @view u[1:4]
    q̇ = @view u[5:8]
    λ = @view u[9:11]
    v = @view du[1:4]
    v̇ = @view du[5:8]
    out[1:4] = v - q̇
    out[5:8] = v̇ - inv(M)*transpose(A(q))*λ + inv(M)*F(1,1,1)
    out[9:11] = Φ(q)
end
u₀ = vcat(q0,q̇0,zero(λ0))
du₀ = vcat(q̇0,zero(q̇0),zero(λ0))
tspan = (0.0,400dt)
differential_vars = vcat(ones(Bool,8),zeros(Bool,3))
prob = DAEProblem(f,du₀,u₀,tspan,differential_vars=differential_vars)
sol = solve(prob,IDA(),dt=dt)
