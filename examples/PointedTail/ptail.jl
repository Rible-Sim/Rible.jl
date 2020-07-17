using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
# using BenchmarkTools
# import PyPlot; const plt = PyPlot
# using LaTeXStrings
# using NLsolve
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot

function make_tail(n)
    nbodies = 2 + 1 + 3n
    nver = n+1
    nhor = 2(n+1)
    nb = nver + nhor
    @assert nb == nbodies
    ver_lengths = zeros(nver)
    hor_lengths = zeros(nhor)
    ver_lengths .= 0.04
    hor_lengths .= 0.04
    ver_index = [i for i=1:nb if mod(i,3)==0]
    hor_index = [i for i=1:nb if mod(i,3)!=0]
    A = zeros(2,3(n+2))
    for i = 1:n+1
        for j = 1:3
            js = 3i+j
            A[:,js] = [0.0,-sum(ver_lengths[1:i])]
        end
    end
    for i = 1:n+1
        js = 3(i-1)+2
        A[:,js-1] = A[:,js] + [-hor_lengths[2i-1],0.0]
        A[:,js+1] = A[:,js] + [ hor_lengths[2i],0.0]
    end

    m = fill(1.3e-1,nb)
    ri = [zeros(2) for i = 1:nb]
    rj = [zeros(2) for i = 1:nb]
    for (i,k) in enumerate(hor_index)
        level = div(k,3)
        js = 3level+2
        if isodd(i)
            ri[k] .= A[:,js]
            rj[k] .= A[:,js-1]
        else
            ri[k] .= A[:,js]
            rj[k] .= A[:,js+1]
        end
    end
    for (i,k) in enumerate(ver_index)
        js = 3(i-1)+2
        ri[k] .= A[:,js]
        rj[k] .= A[:,js+3]
    end
    ri,rj

    movable = ones(Bool,nb)
    constrained = zeros(Bool,nb)
    movable[1:2] .= false
    constrained[1:3] .= true
    function rigidbody(i,movable,constrained,m,ri,rj)
        L = norm(ri-rj)
        CoM = SVector(L/2,0.0)
        inertia = m*L^2/3
        ap1 = SVector(0.0,0.0)
        ap2 = SVector(L,0.0)
        aps = [ap1,ap2]
        prop = TR.RigidBodyProperty(i,movable,
                    m,inertia,
                    CoM,aps;constrained=constrained
                    )
        if constrained
            constrained_index = [1,2]
        else
            constrained_index = Vector{Int}()
        end
        state = TR.RigidBodyState(prop,ri,rj,constrained_index)
        rb = TR.RigidBody(prop,state)
    end
    rbs = [rigidbody(i,movable[i],constrained[i],
                    m[i],ri[i],rj[i]) for i = 1:nb]

    nstrings = 4(n)+2
    original_restlens = sqrt(2)*0.04
    ks = fill(1000.0,nstrings)
    ss = [TR.SString2D(original_restlens,ks[i],0.0) for i = 1:nstrings]

    acs = [TR.Actuator(SVector{1}(ss[i])) for i = 1:nstrings]

    body2q = TR.filter_body2q(TR.build_body2q(rbs),rbs)

    string2ap = Vector{Tuple{TR.ID,TR.ID}}()
    for i = 1:n
        is = 3(i-1)
        push!(string2ap,(TR.ID(is+1,2),TR.ID(is+4,2)))
        push!(string2ap,(TR.ID(is+3,1),TR.ID(is+4,2)))
        push!(string2ap,(TR.ID(is+3,1),TR.ID(is+5,2)))
        push!(string2ap,(TR.ID(is+2,2),TR.ID(is+5,2)))
    end
    push!(string2ap,(TR.ID(nb,2),TR.ID(nb-1,2)))
    push!(string2ap,(TR.ID(nb,2),TR.ID(nb-2,2)))

    cnt = TR.Connectivity(body2q,string2ap)
    tg = TR.TensegrityStructure(rbs,ss,acs,cnt)
    TR.update_strings_apply_forces!(tg)
    tg
end
n = 3
tail = make_tail(n)
q0,q̇0,λ0 = TR.get_initial(tail)
# @code_warntype make_tail(n)

function dynfuncs(tgstruct,q0)

    M = TR.build_massmatrix(tgstruct)
    Φ = TR.build_Φ(tgstruct,q0)
    A = TR.build_A(tgstruct)

    #Q̃=TR.build_Q̃(tgstruct)

    function F!(F,q,q̇,t)
        TR.reset_forces!(tgstruct)
        TR.distribute_q_to_rbs!(tgstruct,q,q̇)
        TR.update_strings_apply_forces!(tgstruct)
        # F .= Q̃*TR.fvector(tgstruct)
        F .= 0.0
        TR.assemble_forces!(F,tgstruct)
    end

    M,Φ,A,F!,nothing
end

M,Φ,A,F!,Jacs = dynfuncs(tail,q0)
q̇0[end-1:end] .= [0.01,0.0]
Φ(q0)
A(q0)*q̇0
dt = 0.01
prob = TS.DyProblem(dynfuncs(tail,q0),q0,q̇0,λ0,(0.0,10.0))
sol = TS.solve(prob,TS.Wendlandt(),dt=dt,ftol=1e-14,verbose=true)
sol = TS.solve(prob,TS.ConstSPARK(1),dt=dt,ftol=1e-12,verbose=true)
