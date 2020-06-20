using LinearAlgebra, Printf
using StaticArrays
using Rotations
using Parameters
using GeometryTypes, AbstractPlotting
using Makie
#using DifferentialEquations
using ForwardDiff
using PATHSolver
using Suppressor
options(convergence_tolerance=1e-14, output=:no,
    time_limit=3600)
using Revise
using SPARK
using TensegrityRobotSimulator
const TRS = TensegrityRobotSimulator
mass = 1.0
inertia = Matrix(1.0I,3,3)
CoM = zeros(3)
r = [0,0,10.5-0.0000001]
R = Matrix(1.0I,3,3)
#R = Matrix(RotXY(0.1,0.3))
#R = Matrix(RotX(float(π)))
ṙ = zeros(3)
#ω = zeros(3)
ω = [100.0,0.0,0.0]
aplist = [[0.5,0.5,0.5], [0.5,-0.5,0.5],
          [-0.5,-0.5,0.5], [-0.5,0.5,0.5],
          [0.5,0.5,-0.5], [0.5,-0.5,-0.5],
          [-0.5,-0.5,-0.5], [-0.5,0.5,-0.5]]
apvector = SVector([TRS.AnchorPoint(aplist[i]) for i = 1:8]...)
prop = TRS.RigidBodyProperty(:cube,:generic,mass,
            SMatrix{3,3}(inertia),
            SVector(CoM...),
            apvector)

state = TRS.RigidBodyState(prop,r,R,ṙ,ω,Val(:NC))

cube = TRS.RigidBody(prop,state)

h = 0.01
tspan = (0.0,0.1)
q0 = copy(cube.state.coords.q)
q̇0 = copy(cube.state.coords.q̇)

function cube_dynamics()
    function M!(mass_matrix,q)
        mass_matrix .= cube.state.cache.M
    end
    function ∂T∂q̇!(p,mass_matrix,q,q̇)
        M!(mass_matrix,q)
        mul!(p,mass_matrix,q̇)
    end

    # External forces (Gravity)
    fG = [0,0,-9.8]
    QG = transpose(cube.state.cache.CG)*fG
    function F!(F,q,q̇,t)
        F .= QG
    end
    function Φ(q)
        TRS.NC.Φ(q)
    end
    function A(q)
        TRS.NC.Φq(q)
    end
    jacs = ()
    A,Φ,∂T∂q̇!,F!,M!,jacs
end
λ0 = zeros(6)
tab = SPARKTableau(1)
cache = SPARKCache(12,6,0.01,(0.0,0.1),1,cube_dynamics())
solution_state = SPARKsolve!(q0,q̇0,λ0,cache,tab)
ts,qs,q̇s = alphastepping(0.01,(0.0,0.03),q0,q̇0,cube)
ts,qs,q̇s = onestepping(h,tspan,q0,q̇0,cube)
ts,qs,q̇s = onestepping(h,tspan,q289,q̇289,cube)
TRS.NC.Φ.(qs)
function energy(cube,q,q̇)
    cube.state.coords.q .= q
    cube.state.coords.q̇ .= q̇
    TRS.coords2state_kinetic!(cube)
    @unpack prop, state = cube
    @unpack mass = prop
    @unpack r = state
    pe = mass*9.8*r[3]
    ke = TRS.kinetic_energy(cube)
    ke,pe,pe + ke
end
es = [energy(cube,solution_state.qs[it],solution_state.q̇s[it]) for it in eachindex(solution_state.ts)]
TRS.NC.Φ(q0)
# Pμ = zeros(6)
# PN = zeros(3ν)
# mlcp = MLCP(AA,B,C,D,a,b,length(a),length(b),Pμ,PN)
