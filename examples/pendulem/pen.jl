using LinearAlgebra
using Parameters
using StaticArrays
using Revise
using NLsolve
using SPARK
using Robot2D
const R2 = Robot2D
function npen(n)
    Ls = zeros(n)
    Ls .= 1.0
    A = zeros(2,n+1)
    for i = 1:n
        A[:,i+1] .= [0.0,-i*Ls[i]]
    end
    ms = ones(n)
    Is = (ms.*Ls.^2)./12
    function rigidbody(i,m,L,I,ri,rj)
        name = Symbol("rb"*string(i))
        type = :generic
        if i == 1
            movable = false
        else
            movable = true
        end

        CoM_x = L/2
        CoM_y = 0.0
        CoM = SVector{2}([CoM_x,CoM_y])

        nap = 2
        ap1 = SVector{2}([0.0,0.0])
        ap2 = SVector{2}([L,0.0])
        anchorpoints = [ap1,ap2]

        prop = R2.RigidBody2DProperty(movable,name,type,
                    m,I,
                    CoM,
                    anchorpoints
                    )
        state = R2.RigidBody2DState(prop,ri,rj)
        rb = R2.RigidBody2D(prop,state)
    end
    rbs = [rigidbody(i,ms[i],Ls[i],
            Is[i],A[:,i],A[:,i+1]) for i = 1:n]

    body2q = R2.build_body2q(rbs)
    cnt = R2.Connectivity(body2q,nothing)
    R2.Structure2D(rbs,Vector{Int}(),Vector{Int}(),cnt)
end

function npen_wend(st2d)
    rbs = st2d.rigidbodies
    cnt = st2d.connectivity
    @unpack body2q = cnt
    M = R2.build_massmatrix(rbs,body2q)
    Φ = R2.build_Φ(st2d)
    A = R2.build_A(st2d)
    function F(q,q̇,t)
        R2.reset_forces!(st2d)
        R2.q2rbstate!(st2d,q,q̇)
        R2.apply_gravity!(st2d)
        ret = similar(q)
        R2.assemble_forces!(ret,st2d)
        ret
    end
    M,Φ,A,F,nothing
end
pen2 = npen(4)
M,Φ,A,F,Jacs = npen_wend(pen2)
q0,q̇0,λ0 = R2.get_initial(pen2,Φ)
# Φ(q0)
# @code_warntype Φ(q0)
# @time Φ(q0)
# A(q0)
# @time A(q0)
# F(q0,q̇0,0.0)
q̇0[end-1] = 0.01
prob = SPARK.DyProblem(npen_wend(deepcopy(pen2)),q0,q̇0,λ0,(0.0,10.0))
dt = 0.01
state = SPARK.solve(prob,dt=dt,ftol=1e-9)
using Plots
pyplot()
plot([state.qs[i][5] for i = 1:length(state.qs)])
plot!([state.qs[i][6]+2 for i = 1:length(state.qs)])
