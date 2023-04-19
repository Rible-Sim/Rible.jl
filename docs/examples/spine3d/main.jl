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
import TensegrityRobots as TR
cd("examples\\spine3d")
includet("define.jl")
# includet("plotting.jl")
# includet("../analysis.jl")
includet("../vis.jl")
spine = spine3d(2;c=10.0)

plot_traj!(spine)
TR.undamped_eigen!

TR.build_Ǩ(spine.tg)

λ0,u0 = TR.inverse_for_restlength(spine,spine)

TR.set_restlen!(spine.tg,u0)
_,λ0 = TR.check_static_equilibrium_output_multipliers(spine.tg)

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
fullq0,_ = TR.get_fullq(spine.tg)
N(fullq0)

function locate_bifurcation_point(bot,γ0)
    (;tg) = bot
    (;ncoords,ncables,ndof,nconstraint) = bot.tg
    (;mvindices,fixindices) = bot.tg.connectivity
    nq = ncoords
    ns = ncables
    nd = ndof
    nλ = nconstraint
    c_val = TR.get_c(tg)
    s_val = TR.get_s(tg)
    ℓ_val = TR.get_cables_len(tg)
    k_val = TR.get_cables_stiffness(tg)
    q_val,_ = TR.get_fullq(tg)
    μ_val = TR.get_cables_restlen(tg)

    @var c[1:length(c_val)]
    @var s[1:length(s_val)]
    @var k[1:length(k_val)]
    @var q[1:length(q_val)]
    @var μ[1:length(μ_val)]
    eprc = 1c
    eprs = 1s
    eprk = 1k
    eprμ = 1μ
    eprq = subs(1q, q[fixindices]=>q_val[fixindices])
    Φ = TR.build_Φ(tg)
    A = TR.build_A(tg)
    𝚽_val = Φ(q_val[mvindices])
    𝐀_val = A(q_val[mvindices])
    Q = TR.build_Q(tg)
    𝐐 = Q(eprq,eprs,eprμ,eprk,eprc)
    𝐐_val = Q(q_val,s_val,μ_val,k_val,c_val)[mvindices]
    # display(𝐀_val*transpose(𝐀_val))
    @assert 𝐀_val*transpose(𝐀_val) ≈ I
    λ_val=𝐀_val*𝐐_val

    @var λ[1:length(λ_val)]
    eprλ = 1λ
    KE = TR.build_KE(tg)
    KG = TR.build_KG(tg)
    Aᵀλ = System(transpose(A(eprq[mvindices]))*eprλ, variables = q[mvindices], parameters = λ)
    # jacobian(Aᵀλ)
    S = TR.build_S(tg)
    𝐒 = S(eprq,eprs,eprc)
    S_val = S(q_val,s_val,c_val)

    𝐊 = KE(eprq,eprs,eprk,eprc)[mvindices,mvindices] +
         KG(eprq,eprs,eprμ,eprk,eprc)[mvindices,mvindices] +
         jacobian(Aᵀλ)
    𝐍 = N(eprq)[mvindices,7:12]
    𝐍_val = N(q_val)[mvindices,7:12]
    # display(𝐀_val*𝐍_val)
    𝑲 = transpose(𝐍)*𝐊*𝐍
    𝑲_val = subs.(𝑲,q[mvindices]=>q_val[mvindices],s=>s_val,λ=>λ_val,c=>c_val,k=>k_val,μ=>μ_val) .|> to_number

    return eigen(𝑲_val).values

    𝐞 = [Matrix(1I,2,2)[:,i] for i in 1:2]
    i = 1; j = 2
    𝑲ij_val = (I-𝐞[i]*transpose(𝐞[i]))*𝑲_val+𝐞[i]*transpose(𝐞[j])
    ξ = 𝑲ij_val\𝐞[i]

    function Sx(qx,sx)
        q_full = subs.(eprq,q[mvindices]=>qx) .|> to_number
        S(q_full,sx,c_val)
    end

    function Qx(qx,sx,γ)
        q_full = subs.(eprq,q[mvindices]=>qx) .|> to_number
        Q(q_full,sx,γ.*ℓ_val,k_val,c_val)[mvindices]
    end

    function 𝑲x(qx,sx,λx,γ)
        subs.(𝑲,q[mvindices]=>qx,s=>sx,λ=>λx,c=>c_val,k=>k_val,μ=>(γ.*ℓ_val)) .|> to_number
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
    x_initial = vcat(q_val[mvindices],s_val,λ_val,ξ,γ0)
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
    @unpack tg = bot
    M = TR.build_massmatrix(tg)
    Φ = TR.build_Φ(tg)
    A = TR.build_A(tg)

    #Q̃=TR.build_Q̃(tg)

    function F!(F,q,q̇,t)
        # du = 0.01*sin(t)
        # TR.actuate!(tg.actuators[1],du)
        TR.reset_forces!(tg)
        TR.distribute_q_to_rbs!(tg,q,q̇)
        TR.update_cables_apply_forces!(tg)
        # TR.apply_gravity!(tg,factor=0.001)
        # F .= Q̃*TR.fvector(tg)
        F .= 0.0
        TR.assemble_forces!(F,tg)
    end

    M,Φ,A,F!,nothing
end


dt = 0.01
prob = TR.SimProblem(spine,dynfuncs,(0.0,5.0))
TR.solve!(prob,TR.Zhong06(),dt=dt,ftol=1e-14,verbose=true)
plotstructure(spine)
plotstructure(spine;do_record=true)
TR.get_cables_tension(spine)

spine.tg.clustecable[1].segs[1].state.length
