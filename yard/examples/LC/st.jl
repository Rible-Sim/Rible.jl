using LinearAlgebra
using SparseArrays
using StaticArrays
using TypeSortedCollections
# using HomotopyContinuation
using DynamicPolynomials
using BenchmarkTools
using NLsolve
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using GeometryBasics
using LaTeXStrings
# using CSV
using Unitful, Printf, Match
using EponymTuples
using Revise
import Rible as RB
cd("examples/LC")
includet("mydef.jl")
includet("../analysis.jl")
includet("../vis.jl")
γ0 = 0.858
twodofbot = man_nd1(2e3;ratio=γ0); twodofbot_plot = deepcopy(twodofbot)

# for setting c
# c = RB.get_c(twodofbot)
# RB.set_C!(twodofbot.st,c)

function N(q)
    rI = @view q[1:2]
    uI = @view q[3:4]
    vI = @view q[5:6]
    rJ = @view q[7:8]
    uJ = (rJ - rI)./0.2
    vJ = (@view q[9:10])./0.2
    uK = @view q[11:12]
    vK = @view q[13:14]
    o2 = zero(uJ)
    n1 = vcat(√2/2*uJ,o2,o2,√2/2*uJ,o2,o2,o2)
    n2 = vcat(2vJ/√6,o2,o2,vJ/√6,uJ/√6,o2,o2)
    n3 = vcat(o2,√2/2*vI,-√2/2*uI,o2,o2,o2,o2)
    n4 = vcat(o2,o2,o2,√2/2*vJ,-√2/2*uJ,o2,o2)
    n5 = vcat(o2,o2,o2,o2,o2,√2/2*vK,-√2/2*uK)
    N = hcat(n1,n2,n3,n4,n5)
end

function locate_bifurcation_point(bot_input,γ0)
    bot = deepcopy(bot_input)
    (;st) = bot
    (;ndof,nconstraints,connectivity) = bot.st
    (;cables) = st.tensiles
    (;nfull,nfree) = connectivity.indexed
    ncables = length(cables)
    nλ = nconstraints
    q̌0 = RB.get_q̌(st)
    q0 = RB.get_q(st)
    c0 = RB.get_c(st)
    ℓ0 = RB.get_cables_len(st)
    s0 = 1 ./ℓ0
    k0 = RB.get_cables_stiffness(st)
    μ0 = RB.get_cables_restlen(st)
    _,λ0 = RB.check_static_equilibrium_output_multipliers(st)
    Φ = RB.make_constraints_function(st,q0)
    A = RB.make_constraints_jacobian(st,q0)
    Q̌ = RB.make_Q̌(st,q0)
    S = RB.make_S(st,q0)
    Ǩm_Ǩg = RB.make_Ǩm_Ǩg(st,q0)
    # state variables
    @polyvar pvq̌[1:length(q̌0)]
    @polyvar pvs[1:length(s0)]
    @polyvar pvλ[1:length(λ0)]
    # parameters
    @polyvar pvc[1:length(c0)]
    @polyvar pvk[1:length(k0)]
    @polyvar pvμ[1:length(μ0)]

    pnq̌ = 1.0pvq̌ .+ 0.0
    pns = 1.0pvs .+ 0.0
    pnλ = 1.0pvλ .+ 0.0
    pnc = 1.0pvc .+ 0.0
    pnk = 1.0pvk .+ 0.0
    pnμ = 1.0pvμ .+ 0.0
    polyΦ = Φ(pnq̌)
    polyA = A(pnq̌)
    polyQ̌ = Q̌(pnq̌,pns,pnμ,pnk,pnc)
    polyS = S(pnq̌,pns,pnc)
    polyQ̌c = transpose(polyA)*pnλ
    polyǨc = reduce(hcat,differentiate.(polyQ̌c,Ref(pvq̌))) |> transpose
    polyǨm, polyǨg = Ǩm_Ǩg(pnq̌,pns,pnμ,pnk,pnc)
    polyǨ = polyǨc .- (polyǨm .+ polyǨg)
    # initial
    nq̌ = length(q̌0)
    ndof = 2
    Ǩ0 = RB.build_Ǩ(bot.st,λ0)
    Ň0 = nullspace(A(q̌0))
    𝒦 = transpose(Ň0)*Ǩ0*Ň0
    𝐞 = [1:ndof .== i for i in 1:ndof]
    i = 2; j = 1
    𝒦ij = (I-𝐞[i]*transpose(𝐞[i]))*𝒦+𝐞[i]*transpose(𝐞[j])
    ξ0 = 𝒦ij\𝐞[i]
    ξ0 /= norm(ξ0)
    # ncables = 0
    function make_bf()
        function inner_bp!(f,x)
            q̌x = @view x[                      1:nfree]
            sx = @view x[                nfree+1:nfree+ncables]
            λx = @view x[        nfree+ncables+1:nfree+ncables+nλ]
            ξx = @view x[     nfree+ncables+nλ+1:nfree+ncables+nλ+ndof]
            γ =        x[nfree+ncables+nλ+ndof+1:end]
            # RB.set_restlen!(st,γ.*ℓ0)
            # RB.update!(st)
            # Qx = Q̌(q̌x,sx,γ.*ℓ0,k0,c0)
            Qx = map(polyQ̌) do z
                    z(
                        pvq̌=>q̌x,
                        pvs=>sx,
                        pvμ=>γ.*ℓ0,
                        pvk=>k0,
                        pvc=>c0
                    )
                end
            # Ax = A(q̌x)
            Ax = map(polyA) do z
                    z(pvq̌=>q̌x)
                end
            # Sx = S(q̌x,sx,c0)
            Sx = map(polyS) do z
                    z(pvq̌=>q̌x,pvs=>sx,pvc=>c0)
                end
            # Φx = Φ(q̌x)
            Φx = map(polyΦ) do z
                    z(pvq̌=>q̌x)
                 end
            Nx = nullspace(Ax)
            # Ǩx = RB.build_Ǩ(st,λx)
            Ǩx = map(polyǨ) do z
                    z(
                        pvq̌=>q̌x,
                        pvs=>sx,
                        pvλ=>λx,
                        pvμ=>γ.*ℓ0,
                        pvk=>k0,
                        pvc=>c0
                    )
                end
            f[                 1:nfree]                 = transpose(Ax)*λx - Qx
            f[           nfree+1:nfree+ncables]         = Sx
            f[   nfree+ncables+1:nfree+ncables+nλ]      = Φx
            f[nfree+ncables+nλ+1:nfree+nλ+ncables+ndof] = transpose(Nx)*Ǩx*Nx*ξx
            f[end]                                      = transpose(ξx)*ξx-1
        end
    end
    f_holder = zeros(nfree+nλ+ncables+ndof+1)
    x_initial = vcat(q̌0,s0,λ0,ξ0,γ0)
    bp! = make_bf()
    # bp!(f_holder,x_initial)
    # return f_holder
    bp = nlsolve(bp!,x_initial,ftol=1e-10,iterations=100,method=:newton)
    # @show
    bp!(f_holder,bp.zero)

    @show f_holder[1:nfree+ncables+nλ] |> norm
    @show f_holder[nfree+ncables+nλ+1:nfree+ncables+nλ+ndof] |> norm
    @show  bp.zero[nfree+ncables+nλ+1:nfree+ncables+nλ+ndof]
    @show  bp.zero[end]
    q_bp = bp.zero[                 1:nfree]
    s_bp = bp.zero[           nfree+1:nfree+ncables]
    λ_bp = bp.zero[   nfree+ncables+1:nfree+ncables+nλ]
    ξ_bp = bp.zero[nfree+ncables+nλ+1:nfree+ncables+nλ+ndof]
    γ_bp = bp.zero[end]
    (q̌=q_bp,s=s_bp,λ=λ_bp,ξ=ξ_bp,γ=γ_bp,isconverged=converged(bp))
end

γ0 = 0.9
pf = locate_bifurcation_point(twodofbot,γ0)
pf.isconverged
pf.γ

pfrc = RB.recover(pf,twodofbot.st)
RB.set_new_initial!(twodofbot_plot,pfrc.q)
plot_traj!(twodofbot_plot)
@code_warntype locate_bifurcation_point(twodofbot,γ0)


@btime locate_bifurcation_point(twodofbot,γ0)

plotstructure(twodofbot)

##


##

q0,q̇0,λ0 = RB.get_initial(st)
M = RB.build_massmatrix(st)
A = RB.build_A(st)
u0 = [s.original_restlen for s in st.cables]
K = diagm([s.k for s in st.cables])
J = [RB.build_Ji(st,i) for i = 1:4]
l = [sqrt(q0'*Ji'*Ji*q0) for Ji in J]
s = 1.0./l

Q̃ = RB.build_Q̃(st)
G = K*([(1-u0[i]*s[i])*Ji*q0 for (i,Ji) in enumerate(J)])
Γ = zeros(2*length(G))
for i = 1:length(G)
    Γ[[2i-1,2i]] = G[i]
end

λ = A(q0)'\(Q̃*Γ)

@show A(q0)'*λ-(Q̃*Γ)


ind = st.connectivity.body2q
spe = diagm([2.0 for i = 1:4])
spe[1,3] = -2.0
spe[2,4]= -2.0
spe[3,1] = -2.0
spe[4,2]= -2.0

Φqq1 = zeros(8,8)
Φqq1[ind[1], ind[1]] = spe
Φqq2 = zeros(8,8)
Φqq2[ind[2], ind[2]] = spe
Φqq3 = zeros(8,8)
Φqq3[ind[3], ind[3]] = spe

Φqq = [zeros(8,8),zeros(8,8),zeros(8,8),Φqq1,Φqq2,Φqq2]
λ'*Φqq

p = K*[(1-u0[i]*s[i])*Ji+0.5u0[i]*s[i]^3*Ji*q0*q0'*Ji'*Ji for (i,Ji) in enumerate(J)]

∂Γ∂q=zeros(8,8)
for i = 1:length(p)
    ∂Γ∂q[[2i-1,2i],:] = p[i]
end

K̄11=λ'*Φqq-Q̃*∂Γ∂q

M̄ = zeros(14,14)
ind = collect(1:8)
ind2 = collect(9:14)
M̄[ind,ind] = M

K̄ = zeros(14,14)
K̄[ind,ind]=K̄11
K̄[ind,ind2]=A(q0)'
K̄[ind2,ind]=A(q0)



v,d = eigen(K̄,M̄)
v,K̄,M̄
#end

U0 = zeros(28,28)
ind = collect(1:14)
U0[ind.+14,ind] = M̄
U0[ind,ind.+14] = M̄

V0 = zeros(28,28)
V0[ind,ind] = -K̄
V0[ind.+14,ind.+14] = M̄
v2,d2 = eigen(V0,U0)
