using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using Parameters
#using DifferentialEquations
#using ForwardDiff
using NLsolve
using Cubature
using Makie
AbstractPlotting.__init__()
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using Robot
const TR = Robot


function compute_offset(n,α,Δr,a1,a2)
    R0 = RotZ(α)
    function find_zs()
        function make_f(r)
            dr(z) = sqrt(1+(2*a2*z+a1)^2)
            function f(Z)
                int,_ = hquadrature(dr,0,Z)
                int-r
            end
        end
        ret = zeros(n)
        for i = 1:n
            result = nlsolve(make_f(Δr*(i-1)),[Δr*(i-1)])
            ret[i:i] .= result.zero
        end
        ret
    end
    zs = find_zs()
    ys = a2.*zs.^2 .+ a1.*zs
    rs = [R0*[0,y,z] for (y,z) in zip(ys,zs)]
    #
    function make_R1(z)
        dy_dz(z) = 2a2*z + a1
        θ = atan(dy_dz(z))
        b = [1,0,0]
        τ = [0,sin(θ),cos(θ)]
        η = τ×b
        R1 = hcat(b,η,τ)
    end
    Rs = [
        R0*make_R1(z)
        for z in zs
    ]
    Rs[1] = Matrix(1.0I,3,3)
    rs,Rs
end
function links(n,a1=0.0,a2=0.0,α=0.0)
    nbodies = n
    nbp = 4*n

    a = 0.08 #m
    h = 0.1 #m
    θ = 2π/3
    l = √(a^2+h^2)
    b = √3a

    mass = 0.1 #kg
    #inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-8 # N/m^2
    inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-2
    mass_locus = [0.0, 0.0, 0.0]


    ap1 = SVector{3}([0.0, 0.0, 0.0])
    ap2 = SVector{3}([a, 0.0, h])
    ap3 = SVector{3}([a*cos(θ), a*sin(θ), h])
    ap4 = SVector{3}([a*cos(θ), -a*sin(θ), h])
    aps = [ap1,ap2,ap3,ap4]

    contactable = ones(Bool,n)
    contactable[1] = 0

    props = [RB.RigidBodyProperty(i,contactable[i],mass,
                SMatrix{3,3}(inertia),
                SVector{3}(mass_locus),
                aps) for i = 1:n]

    dr = 0.5h
    # r = [ [0.0,0.0,0.0+(i-1)*dr] for i = 1:n]
    # R = [Matrix(1.0I,3,3) for i = 1:n]
    r,R = compute_offset(n,α,dr,a1,a2)
    ṙ = [zeros(3) for i = 1:n]
    ω = [zeros(3) for i = 1:n]

    function rigidbody(prop,aps,r,R,ṙ,ω)
        ri,rj,rk,rl = [r+R*ap for ap in aps]
        nmcs,q,q̇ = RB.NCF.BP4P(ri,rj,rk,rl,r,R,ṙ,ω)
        state = RB.RigidBodyState(prop,r,R,ṙ,ω,q,q̇)
        @show R
        body = RB.RigidBody(prop,state)
    end
    rbs = [rigidbody(props[i],aps,r[i],R[i],ṙ[i],ω[i]) for i = 1:n]

    nstrings = 7*(n-1)
    stringlenH = 0.5h
    stringlenR = 0.09433981132056604
    stringlens = repeat(vcat(fill(stringlenH,4),fill(stringlenR,3)),n-1)
    kH = 1e1
    kR = 1e1
    ks = repeat(vcat(fill(kH,4),fill(kR,3)),n-1)
    c = 1.0
    cs = repeat(fill(c,7),n-1)
    ss = [RB.SString3D(stringlens[i],ks[i],cs[i]) for i = 1:nstrings]
    acs = [RB.Actuator(SVector{1}(ss[7(i-1)+j])) for j = 2:7 for i = 1:n-1]
    bodynq = RB.get_num_of_coords(rbs[1])
    bodyid2q = [(bodynq*(i-1)+1:bodynq*i) for i = 1:n]


    string2ap = Vector{Tuple{RB.ID,RB.ID}}()
    for i = 1:n-1
        for j = 1:4
            push!(string2ap,(RB.ID(i,j),RB.ID(i+1,j)))
        end
        for j = 2:4
            push!(string2ap,(RB.ID(i,j),RB.ID(i+1,1)))
        end
    end
    cnt = RB.Connectivity(bodyid2q,string2ap)
    st = RB.Structure(rbs,ss,acs,cnt)
    RB.update_strings_apply_forces!(st)
    st
end
nlink = 2
linkn = links(nlink)
q0,q̇0,λ0 = RB.get_initial(linkn)
# @code_warntype links(nlink)
α = 0.0
a2 = 1
a1 = 0.0
# compute_offset(nlink,α,0.1*0.5,a1,a2)

reflinkn = links(nlink,a1,a2)
# @code_warntype links(nlink,a1,a2)

function inverse(tgstruct,refstruct)
    function ikfuncs(tgstruct)

        A = RB.build_A(tgstruct)

        Q̃=RB.build_Q̃(tgstruct)

        function F!(F,u)
            RB.reset_forces!(tgstruct)
            RB.actuate!(tgstruct,u)
            RB.update_strings_apply_forces!(tgstruct)
            F .= Q̃*RB.fvector(tgstruct)
        end

        A,F!
    end
    q0,q̇0,λ0 = RB.get_initial(refstruct)
    RB.distribute_q_to_rbs!(tgstruct,q0,q̇0)
    nu = length(tgstruct.actuators)
    u0 = zeros(nu)
    ikprob = TS.IKProblem(ikfuncs(tgstruct),q0,u0,λ0)
    RB.iksolve(ikprob)
end
u,_ = inverse(linkn,reflinkn)
RB.actuate!(linkn,zero(u))
# RB.actuate!(linkn,zeros(12))
# @code_warntype links(nlink)
# link3.ndim
# link3.connectivity.bodyid2q[end][end]

function dynfuncs(st)
    M = RB.build_massmatrix(st)
    Φ = RB.build_Φ(st)
    A = RB.build_A(st)
    Q̃ = RB.build_Q̃(st)
    function F!(F,q,q̇,t)
        RB.reset_forces!(st)
        RB.distribute_q_to_rbs!(st,q,q̇)
        RB.update_strings_apply_forces!(st)
        F .= Q̃*RB.fvector(st)
        # RB.assemble_forces!(F,st)
    end
    M,Φ,A,F!,nothing
end

function linearload(st,u,tend)
    M = RB.build_massmatrix(st)
    Φ = RB.build_Φ(st)
    A = RB.build_A(st)
    Q̃ = RB.build_Q̃(st)
    function F!(F,q,q̇,t)
        ut = u/tend*t
        RB.reset_forces!(st)
        #RB.actuate!(st,ut)
        RB.distribute_q_to_rbs!(st,q,q̇)
        RB.update_strings_apply_forces!(st)
        F .= Q̃*RB.fvector(st)
        #RB.assemble_forces!(F,st)
    end
    M,Φ,A,F!,nothing
end
#
# M,Φ,A,F!,Jacs = dynfuncs(linkn)
# Φ(q0)
# @code_warntype Φ(q0)
# A(q0)*q̇0
# @code_warntype A(q0)



dt = 0.01
tend = 10.0
prob = TS.DyProblem(linearload(linkn,u,tend),q0,q̇0,λ0,(0.0,tend))
sol = TS.solve(prob,TS.Wendlandt(),dt=dt,ftol=1e-13,verbose=true)
sol = TS.solve(prob,TS.ConstSPARK(1),dt=dt,ftol=1e-12,verbose=true)


M,Φ,A,F!,Jacs = dynfuncs(linkn)
dt = 0.01
prob = TS.DyProblem(dynfuncs(linkn),q0,q̇0,λ0,(0.0,10.0))
sol = TS.solve(prob,TS.Wendlandt(),dt=dt,ftol=1e-12,verbose=true)



q0,q̇0,λ0 = RB.get_initial(reflinkn)
M,Φ,A,F!,Jacs = dynfuncs(reflinkn)
dt = 0.01
prob = TS.DyProblem(dynfuncs(reflinkn),q0,q̇0,λ0,(0.0,3.0))
sol = TS.solve(prob,TS.Wendlandt(),dt=dt,ftol=1e-12,verbose=true)
