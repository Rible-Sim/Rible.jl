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
using TensegrityRobot
const TR = TensegrityRobot


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

    mass = 0.025 #kg
    #inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-8 # N/m^2
    inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-3
    r̄g = [0.0, 0.0, 0.0]


    ap1 = SVector{3}([0.0, 0.0, 0.0])
    ap2 = SVector{3}([a, 0.0, h])
    ap3 = SVector{3}([a*cos(θ), a*sin(θ), h])
    ap4 = SVector{3}([a*cos(θ), -a*sin(θ), h])
    aps = [ap1,ap2,ap3,ap4]

    movable = ones(Bool,n)
    movable[1] = 0

    props = [TR.RigidBodyProperty(i,movable[i],mass,
                SMatrix{3,3}(inertia),
                SVector{3}(r̄g),
                aps) for i = 1:n]

    dr = 0.5h
    # r = [ [0.0,0.0,0.0+(i-1)*dr] for i = 1:n]
    # R = [Matrix(1.0I,3,3) for i = 1:n]
    r,R = compute_offset(n,α,dr,a1,a2)
    ṙ = [zeros(3) for i = 1:n]
    ω = [zeros(3) for i = 1:n]

    function rigidbody(prop,aps,r,R,ṙ,ω)
        ri,rj,rk,rl = [r+R*ap for ap in aps]
        lncs,q,q̇ = TR.NCF.BP4P(ri,rj,rk,rl,r,R,ṙ,ω)
        state = TR.RigidBodyState(prop,lncs,r,R,ṙ,ω,q,q̇)
        @show R
        rb = TR.RigidBody(prop,state)
    end
    rbs = [rigidbody(props[i],aps,r[i],R[i],ṙ[i],ω[i]) for i = 1:n]

    nstrings = 7*(n-1)
    stringlenH = 0.5h
    stringlenR = 0.09433981132056604-0.006163534339610306
    stringlens = repeat(vcat(fill(stringlenH,4),fill(stringlenR,3)),n-1)
    kH = 1e3
    kR = 1e3
    ks = repeat(vcat(fill(kH,4),fill(kR,3)),n-1)
    c = 1000.0
    cs = repeat(fill(c,7),n-1)
    ss = [TR.SString3D(stringlens[i],ks[i],cs[i]) for i = 1:nstrings]
    acs = [TR.Actuator(SVector{1}(ss[7(i-1)+j])) for j = 2:7 for i = 1:n-1]
    bodynq = TR.get_nbodycoords(rbs[1])
    body2q = [(bodynq*(i-1)+1:bodynq*i) for i = 1:n]


    string2ap = Vector{Tuple{TR.ID,TR.ID}}()
    for i = 1:n-1
        for j = 1:4
            push!(string2ap,(TR.ID(i,j),TR.ID(i+1,j)))
        end
        for j = 2:4
            push!(string2ap,(TR.ID(i,j),TR.ID(i+1,1)))
        end
    end
    cnt = TR.Connectivity(body2q,string2ap)
    tg = TR.TensegrityStructure(rbs,ss,acs,cnt)
    TR.update_strings_apply_forces!(tg)
    tg
end
nlink = 2
linkn = links(nlink)
q0,q̇0,λ0 = TR.get_initial(linkn)
#
# function dynfuncs(tg,q0)
#     M = TR.build_massmatrix(tg)
#     Φ = TR.build_Φ(tg,q0)
#     A = TR.build_A(tg)
#     Q̃ = TR.build_Q̃(tg)
#     function F!(F,q,q̇,t)
#         TR.reset_forces!(tg)
#         TR.distribute_q_to_rbs!(tg,q,q̇)
#         TR.update_strings_apply_forces!(tg)
#         F .= Q̃*TR.fvector(tg)
#         # TR.assemble_forces!(F,tg)
#     end
#     M,Φ,A,F!,nothing
# end
# M,Φ,A,F!,Jacs = dynfuncs(linkn,q0)
# dt = 0.01
# prob = TS.DyProblem(dynfuncs(linkn,q0),q0,q̇0,λ0,(0.0,10.0))
# sol = TS.solve(prob,TS.Wendlandt(),dt=dt,ftol=1e-13,verbose=true)
#

a2 = 1.1
a1 = 0.00
reflinkn = links(nlink,a1,a2)
qr,q̇r,λr = TR.get_initial(reflinkn)
#
# q0,q̇0,λ0 = TR.get_initial(reflinkn)
# M,Φ,A,F!,Jacs = dynfuncs(reflinkn,q0)
# dt = 0.01
# prob = TS.DyProblem(dynfuncs(reflinkn,q0),q0,q̇0,λ0,(0.0,3.0))
# sol = TS.solve(prob,TS.Wendlandt(),dt=dt,ftol=1e-12,verbose=true)
#

function compensate_gravity(tgstruct)
    # reset velocity
    q0,q̇0,λ0 = TR.get_initial(tgstruct)
    TR.distribute_q_to_rbs!(tgstruct,q0,zero(q̇0))
    nu = length(tgstruct.actuators)
    u0 = zeros(nu)
    ikprob = TS.IKProblem(TR.compensate_gravity_funcs(tgstruct),q0,u0,λ0)
    TR.iksolve(ikprob)
end

# ug,_ = compensate_gravity(linkn)
# ug[4:6]
ur,_ = compensate_gravity(reflinkn)
lensr = TR.get_strings_len(reflinkn,qr)
lens0 = TR.get_strings_len(linkn,q0)
dlens = lensr[2:7]-lens0[2:7]

function no_more_extending!(ur,dlens)
    for i in eachindex(ur)
        if ur[i] > 0.0
            if ur[i] > dlens[i]
                ur[i] = dlens[i]
            end
        end
    end
    ur
end
no_more_extending!(ur,dlens)

function inverse(tgstruct,refstruct)
    function ikfuncs(tgstruct)

        A = TR.build_A(tgstruct)

        Q̃=TR.build_Q̃(tgstruct)

        function F!(F,u)
            TR.reset_forces!(tgstruct)
            TR.actuate!(tgstruct,u)
            TR.update_strings_apply_forces!(tgstruct)
            TR.apply_gravity!(tgstruct)
            F .= 0
            #F .= Q̃*TR.fvector(tgstruct)
            TR.assemble_forces!(F,tgstruct)
        end

        A,F!
    end
    q0,q̇0,λ0 = TR.get_initial(refstruct)
    TR.distribute_q_to_rbs!(tgstruct,q0,zero(q̇0))
    nu = length(tgstruct.actuators)
    u0 = zeros(nu)
    ikprob = TS.IKProblem(ikfuncs(tgstruct),q0,u0,λ0)
    TR.iksolve(ikprob)
end
u,_ = inverse(linkn,reflinkn)
TR.actuate!(reflinkn,u)
TR.update_strings_apply_forces!(reflinkn)
[s.state.tension for s in reflinkn.strings]

TR.actuate!(linkn,ug)
TR.actuate!(linkn,ur)
function linearload(tg,q0,initial_u,target_u,target_t)
    M = TR.build_massmatrix(tg)
    Φ = TR.build_Φ(tg,q0)
    A = TR.build_A(tg)
    Q̃ = TR.build_Q̃(tg)
    function F!(F,q,q̇,t)
        if t < target_t
            ut = initial_u + t/target_t*(target_u-initial_u)
        else
            ut = target_u
        end
        TR.actuate!(tg,ut)
        TR.reset_forces!(tg)
        TR.distribute_q_to_rbs!(tg,q,q̇)
        TR.update_strings_apply_forces!(tg)
        # tensions = [s.state.tension for s in tg.strings]
        # @show tensions
        TR.apply_gravity!(tg)
        #F .= Q̃*TR.fvector(tg)
        F .= 0
        TR.assemble_forces!(F,tg)
    end
    M,Φ,A,F!,nothing
end
dt = 0.01
target_t = 5.0
prob = TS.DyProblem(linearload(linkn,q0,zero(ur),ur,target_t),q0,q̇0,λ0,(0.0,15.0))
prob = TS.DyProblem(linearload(linkn,qr,u,tend),qr,q̇r,λr,(0.0,tend))

sol = TS.solve(prob,TS.Wendlandt(),dt=dt,ftol=1e-14,verbose=true)
TR.distribute_q_to_rbs!(linkn,qr,q̇r)
TR.actuate!(linkn,u)
TR.update_strings_apply_forces!(linkn)
[s.state.tension for s in linkn.strings]
