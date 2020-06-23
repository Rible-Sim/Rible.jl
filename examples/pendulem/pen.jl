using LinearAlgebra
using Parameters
using StaticArrays
using Revise
using NLsolve
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot
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
        if i == 1
            movable = false
        else
            movable = true
        end

        if i == 2
            constrained = true
            constrained_index = [1,2]
        else
            constrained = false
            constrained_index = Vector{Int}()
        end

        CoM = SVector(L/2,0.0)

        nap = 2
        ap1 = SVector(0.0,0.0)
        ap2 = SVector(L,  0.0)
        aps = [ap1,ap2]

        prop = TR.RigidBodyProperty(i,movable,
                    m,I,
                    CoM,
                    aps;constrained=constrained
                    )
        state = TR.RigidBodyState(prop,ri,rj,constrained_index)
        rb = TR.RigidBody(prop,state)
    end
    rbs = [rigidbody(i,ms[i],Ls[i],
            Is[i],A[:,i],A[:,i+1]) for i = 1:n]

    body2q = TR.filter_body2q(TR.build_body2q(rbs),rbs)
    cnt = TR.Connectivity(body2q,nothing)
    TR.TensegrityStructure(rbs,Vector{Int}(),Vector{Int}(),cnt)
end

pen2 = npen(2)
q0,q̇0,λ0 = TR.get_initial(pen2)
function dynfuncs(tgstruct,q0)
    M = TR.build_massmatrix(tgstruct)
    Φ = TR.build_Φ(tgstruct,q0)
    A = TR.build_A(tgstruct)
    function F!(F,q,q̇,t)
        TR.reset_forces!(tgstruct)
        TR.distribute_q_to_rbs!(tgstruct,q,q̇)
        TR.apply_gravity!(tgstruct)
        F .= 0.0
        TR.assemble_forces!(F,tgstruct)
    end
    M,Φ,A,F!,nothing
end
M,Φ,A,F!,Jacs = dynfuncs(pen2,q0)
Φ(q0)
# @code_warntype Φ(q0)
# @time Φ(q0)
A(q0)
# @time A(q0)
F = similar(q0)
F!(F,q0,q̇0,0.0)
q̇0[end-1] = 0.01
prob = TS.DyProblem(dynfuncs(pen2,q0),q0,q̇0,λ0,(0.0,10.0))
dt = 0.01
sol = TS.solve(prob,TS.Wendlandt(),dt=dt,ftol=1e-14,verbose=true)
using Plots
pyplot()
plot([sol.qs[i][4] for i = 1:length(sol.qs)])
plot!([sol.qs[i][4]+2 for i = 1:length(sol.qs)])
[sol.qs[i][4] for i = 1:length(sol.qs)]
