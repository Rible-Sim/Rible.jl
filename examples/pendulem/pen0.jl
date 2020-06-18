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
                2,
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
    function F!(F,q,q̇,t)
       F .=const_F
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
    M,Φ,A,F!,nothing
end
M,Φ,A,F!,Jacs = penfunc()
θ = -π/3
q0 = [0,0,cos(θ),sin(θ)]
q̇0 = zeros(4)
λ0 = zeros(3)
nq = length(q0)
nλ = length(λ0)
dt = 0.01
prob = SPARK.DyProblem(penfunc(),q0,q̇0,λ0,(0.0,0.02))

state = SPARK.solve(prob,dt=dt,ftol=1e-13)
using Traceur
@trace SPARK.solve(prob,dt=dt,ftol=1e-13)
