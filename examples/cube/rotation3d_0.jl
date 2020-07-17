using LinearAlgebra
using StaticArrays
using SparseArrays
using Parameters
using Rotations
using Cthulhu
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot
import PyPlot; const plt = PyPlot

m = 1.0
CoM = zeros(3)
inertia = Matrix(1.0I,3,3)
r0 = [0.0,0.0,2.0]
ṙ0 = zeros(3)
R0 = one(RotMatrix{3})
ω0 = rand(3)

aps = [[cos(i*π/2+π/4),sin(i*π/2+π/4),h] for h in [0.5,-0.5] for i = 0:3 ]
prop = TR.RigidBodyProperty(1,true,m,inertia,
                    CoM,aps)
ri = copy(r0)
bps1,q1,q̇1 = TR.NaturalCoordinates.BP1P3V(ri,r0,R0,ṙ0,ω0)
state1 = TR.RigidBodyState(prop,bps1,r0,R0,ṙ0,ω0,q1,q̇1)
rb1 = TR.RigidBody(prop,state1)
body2q = [collect(1:12)]
contacts = [TR.ID(1,i) for i = 1:8]
cnt = TR.Connectivity(body2q,nothing,contacts)
tgrb1 = TR.TensegrityStructure([rb1],Vector{Int}(),Vector{Int}(),cnt)
q0,q̇0,λ0 = TR.get_initial(tgrb1)
function dynfuncs(tgstruct,q0)
    M = TR.build_massmatrix(tgstruct)
    Φ = TR.build_Φ(tgstruct,q0)
    A = TR.build_A(tgstruct)
    function F!(F,q,q̇,t)
        TR.reset_forces!(tgstruct)
        TR.apply_gravity!(tgstruct)
        TR.assemble_forces!(F,tgstruct)
    end

    M,Φ,A,F!,nothing
end

M,Φ,Φq,F! = dynfuncs(tgrb1,q0)
M
Φ(q0)
Φq(q0)
function confuncs(tgstruct)
    ncontact = 8
    ϵ = 0.1
    μ = 0.5
    afcs = [TR.ApproxFrictionalContact(ϵ,μ,24) for i = 1:ncontact]
    function get_gaps(q)
        TR.distribute_q_to_rbs!(tgstruct,q,zero(q))
        gaps = Vector{Float64}()
        for ap in tgstruct.rigidbodies[1].state.p
            push!(gaps,ap[3])
        end
        gaps
    end
    function update_contacts(q,active_index)
        for i = 1:ncontact
            if i in active_index
                afcs[i].impulse.active = true
            else
                afcs[i].impulse.active = false
            end
        end
    end
    function get_HNb(afcs,active_index,q̇)
        active_afcs = afcs[active_index]
        nactive = length(active_index)
        ncoords = tgstruct.ncoords
        nactive_βs = sum(afc.m for afc in active_afcs)
        ncomplementarity = nactive+nactive_βs+nactive
        H = zeros(ncoordinates,ncomplementarity)
        b = spzeros(ncomplementarity)
        is = nactive
        for i = 1:nactive
            icontact = tgstruct.connectivity.contacts[active_index[i]]
            rbid = icontact.rbid
            apid = icontact.apid
            Ci = tgstruct.rigidbodies[rbid].state.cache.Cp[apid]
            ni = active_afcs[i].frame.n
            Di = active_afcs[i].D
            H[:,i] = transpose(Ci)*ni
            mi = active_afcs[i].m
            H[:,is+1:is+mi] = transpose(Ci)*Di
            is += mi
            b[i] = active_afcs[i].ϵ*transpose(ni)*Ci*q̇
        end
        N = spzeros(ncomplementarity,ncomplementarity)
        js = nactive
        for i = 1:nactive
            is = nactive+nactive_βs
            N[is+i,i] = active_afcs[i].μ
            mi = active_afcs[i].m
            ei = active_afcs[i].e
            N[is+i,js+1:js+mi] = -ei
            js += mi
        end
        N[nactive+1:nactive+nactive_βs,
          nactive+nactive_βs+1:ncomplementarity] =
         -transpose(N[nactive+nactive_βs+1:ncomplementarity,
            nactive+1:nactive+nactive_βs])
        H,N,b
    end
    afcs,get_gaps,update_contacts,get_HNb
end
afcs,get_gaps,update_contacts,get_HNb = confuncs(tgrb1)
# @code_warntype  confuncs(tgrb1)
# update_contacts(q0,get_gaps(q0))
# @code_warntype update_contacts(q0,get_gaps(q0))
# afcs
active_index = findall(get_gaps(q0).<=0.0)
H,N,b = get_HNb(afcs,active_index,q̇0)
dt = 0.001
prob = TS.DyProblem(dynfuncs(tgrb1,q0),q0,q̇0,λ0,(0.0,3.0))
sol = TS.solve(prob,TS.Linear(confuncs(tgrb1)),dt=dt,ftol=1e-14,verbose=true)

@descend  TS.solve(prob,TS.Linear(confuncs(tgrb1)),dt=dt,ftol=1e-14,verbose=true)
TR.kinetic_energy_coords(rb1.state)

X_kes = [TR.kinetic_energy_coords(rb1,q̇) for q̇ in sol.q̇s]
X_pes = [TR.gravity_potential_energy(rb1,q) for q in sol.qs]

X_es = X_kes .+ X_pes

plt.plot((X_es.-X_es[1])./X_es[1])

TR.ApproxFrictionalContact(1.0,0.5,12)
fc = TR.ApproxFrictionalContact(1.0,0.5,12)
