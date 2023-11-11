using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
# using Parameters
using Makie
import GLMakie
GLMakie.activate!()
using TypeSortedCollections
#using DifferentialEquations
#using ForwardDiff
using DynamicPolynomials
# using HomotopyContinuation
using Printf
using CoordinateTransformations
using Revise
import Rible as RB
cd("examples\\spine3d")
includet("define.jl")
# includet("plotting.jl")
# includet("../analysis.jl")
includet("../vis.jl")
spine = spine3d(2;c=10.0)

plot_traj!(spine)
RB.undamped_eigen!

RB.build_Ǩ(spine.st)

λ0,u0 = RB.inverse_for_restlength(spine,spine)

RB.set_restlen!(spine.st,u0)
_,λ0 = RB.check_static_equilibrium_output_multipliers(spine.st)

function N(q)
    rI = @view q[1:3]
    uI = @view q[4:6]
    vI = @view q[7:9]
    wI = @view q[10:12]
    rJ = @view q[13:15]
    uJ = @view q[16:18]
    vJ = @view q[19:21]
    wJ = @view q[22:24]
    o3 = zero(uJ)
    (-√2/2)wI
    N1 = [
        uI vI wI       o3       o3         o3;
        o3 o3 o3       o3 -√2/2*wI   -√2/2*vI;
        o3 o3 o3  √2/2*wI       o3   -√2/2*uI;
        o3 o3 o3 -√2/2*vI  √2/2*vI         o3;
    ]
    N2 = [
        uJ vJ wJ       o3       o3       o3;
        o3 o3 o3       o3 -√2/2*wJ -√2/2*vJ;
        o3 o3 o3  √2/2*wJ       o3 -√2/2*uJ;
        o3 o3 o3 -√2/2*vJ  √2/2*vJ       o3;
    ]

    N = [N1  zero(N1) ;
        zero(N2) N2;]
end
fullq0,_ = RB.get_fullq(spine.st)
N(fullq0)

function locate_bifurcation_point(bot,γ0)
    (;st) = bot
    (;ncoords,ncables,num_of_dof,num_of_cstr) = bot.st
    (;mvidx,fixidx) = bot.st.connectivity
    nq = ncoords
    ns = ncables
    nd = num_of_dof
    nλ = num_of_cstr
    c_val = RB.get_local_coords(st)
    s_val = RB.get_s(st)
    ℓ_val = RB.get_cables_len(st)
    k_val = RB.get_cables_stiffness(st)
    q_val,_ = RB.get_fullq(st)
    μ_val = RB.get_cables_restlen(st)

    @var c[1:length(c_val)]
    @var s[1:length(s_val)]
    @var k[1:length(k_val)]
    @var q[1:length(q_val)]
    @var μ[1:length(μ_val)]
    eprc = 1c
    eprs = 1s
    eprk = 1k
    eprμ = 1μ
    eprq = subs(1q, q[fixidx]=>q_val[fixidx])
    Φ = RB.build_Φ(st)
    A = RB.build_A(st)
    𝚽_val = Φ(q_val[mvidx])
    𝐀_val = A(q_val[mvidx])
    Q = RB.build_Q(st)
    𝐐 = Q(eprq,eprs,eprμ,eprk,eprc)
    𝐐_val = Q(q_val,s_val,μ_val,k_val,c_val)[mvidx]
    # display(𝐀_val*transpose(𝐀_val))
    @assert 𝐀_val*transpose(𝐀_val) ≈ I
    λ_val=𝐀_val*𝐐_val

    @var λ[1:length(λ_val)]
    eprλ = 1λ
    KE = RB.build_KE(st)
    KG = RB.build_KG(st)
    Aᵀλ = System(transpose(A(eprq[mvidx]))*eprλ, variables = q[mvidx], parameters = λ)
    # jacobian(Aᵀλ)
    S = RB.build_S(st)
    𝐒 = S(eprq,eprs,eprc)
    S_val = S(q_val,s_val,c_val)

    𝐊 = KE(eprq,eprs,eprk,eprc)[mvidx,mvidx] +
         KG(eprq,eprs,eprμ,eprk,eprc)[mvidx,mvidx] +
         jacobian(Aᵀλ)
    𝐍 = N(eprq)[mvidx,7:12]
    𝐍_val = N(q_val)[mvidx,7:12]
    # display(𝐀_val*𝐍_val)
    𝑲 = transpose(𝐍)*𝐊*𝐍
    𝑲_val = subs.(𝑲,q[mvidx]=>q_val[mvidx],s=>s_val,λ=>λ_val,c=>c_val,k=>k_val,μ=>μ_val) .|> to_number

    return eigen(𝑲_val).values

    𝐞 = [Matrix(1I,2,2)[:,i] for i in 1:2]
    i = 1; j = 2
    𝑲ij_val = (I-𝐞[i]*transpose(𝐞[i]))*𝑲_val+𝐞[i]*transpose(𝐞[j])
    ξ = 𝑲ij_val\𝐞[i]

    function Sx(qx,sx)
        q_full = subs.(eprq,q[mvidx]=>qx) .|> to_number
        S(q_full,sx,c_val)
    end

    function Qx(qx,sx,γ)
        q_full = subs.(eprq,q[mvidx]=>qx) .|> to_number
        Q(q_full,sx,γ.*ℓ_val,k_val,c_val)[mvidx]
    end

    function 𝑲x(qx,sx,λx,γ)
        subs.(𝑲,q[mvidx]=>qx,s=>sx,λ=>λx,c=>c_val,k=>k_val,μ=>(γ.*ℓ_val)) .|> to_number
    end

    function bp!(f,x)
        qx = x[            1:nq]
        sx = x[         nq+1:nq+ns]
        λx = x[      nq+ns+1:nq+ns+nλ]
        ξx = x[   nq+ns+nλ+1:nq+ns+nλ+nd]
        γ =  x[nq+ns+nλ+nd+1:end]
        f[         1:nq]          = transpose(A(qx))*λx - Qx(qx,sx,γ)
        f[      nq+1:nq+ns]       = Sx(qx,sx)
        f[   nq+ns+1:nq+ns+nλ]    = Φ(qx)
        f[nq+ns+nλ+1:nq+ns+nλ+nd] = 𝑲x(qx,sx,λx,γ)*ξx
        # @show f[nq+ns+nλ+1:nq+ns+nλ+nd], γ
        f[end]                    = transpose(𝐞[j])*ξx-1
    end
    f_holder = zeros(2nq+ns+1)
    x_initial = vcat(q_val[mvidx],s_val,λ_val,ξ,γ0)
    # eigen(to_number.(𝐾_val))
    bp!(f_holder,x_initial)
    # @show f_holder[nq+ns+nλ+1:nq+ns+nλ+nd],x_initial[nq+ns+nλ+1:nq+ns+nλ+nd]
    bp = nlsolve(bp!,x_initial,ftol=1e-12,iterations=100,method=:newton)
    bp!(f_holder,bp.zero)
    q_bp = bp.zero[         1:nq]
    s_bp = bp.zero[      nq+1:nq+ns]
    λ_bp = bp.zero[   nq+ns+1:nq+ns+nλ]
    ξ_bp = bp.zero[nq+ns+nλ+1:nq+ns+nλ+nd]
    γ_bp = bp.zero[end]
    # @show f_holder[nq+ns+nλ+1:nq+ns+nλ+nd],ξ_bp
    # @show eigen(𝑲x(q_bp,s_bp,λ_bp,γ_bp)).values
    (q=q_bp,s=s_bp,λ=λ_bp,ξ=ξ_bp,γ=γ_bp)
end

locate_bifurcation_point(spine,0.0)

function dynfuncs(bot)
    @unpack st = bot
    M = RB.build_massmatrix(st)
    Φ = RB.build_Φ(st)
    A = RB.build_A(st)

    #Q̃=RB.build_Q̃(st)

    function F!(F,q,q̇,t)
        # du = 0.01*sin(t)
        # RB.actuate!(st.actuators[1],du)
        RB.reset_forces!(st)
        RB.distribute_q_to_rbs!(st,q,q̇)
        RB.update_cables_apply_forces!(st)
        # RB.apply_gravity!(st,factor=0.001)
        # F .= Q̃*RB.fvector(st)
        F .= 0.0
        RB.assemble_forces!(F,st)
    end

    M,Φ,A,F!,nothing
end


dt = 0.01
prob = RB.SimProblem(spine,dynfuncs,(0.0,5.0))
RB.solve!(prob,RB.Zhong06(),dt=dt,ftol=1e-14,verbose=true)
plotstructure(spine)
plotstructure(spine;do_record=true)
RB.get_cables_tension(spine)

spine.st.clustecable[1].segs[1].state.length
