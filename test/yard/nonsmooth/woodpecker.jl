include("deps.jl")
using AbbreviatedStackTraces #jl
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true #jl
ENV["JULIA_STACKTRACE_MINIMAL"] = true #jl
import Rible as RB
figdir::String = joinpath(pathof(RB),"../../tmp")
if Sys.iswindows() #src
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\IMSD 2025\LaTex_Abstract" #src
elseif Sys.isapple() #src
    figdir::String = raw"/Users/jacob/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/Papers/IMSD 2024/LaTex_Abstract" #src
end #src
include("../../vis.jl")
includet("../../vis.jl") #jl

tw = 455.8843 #pt |> pt2px

#--  slider crank 
include("../../../examples/robots/woodpecker.jl")
includet("../../../examples/robots/woodpecker.jl")#jl

wt = woodpecker()

q = RB.get_coords(wt.structure)
N_in = RB.intrinsic_nullspace(wt.structure,q)
A = RB.make_cstr_jacobian(wt.structure)(q)
N_ex = RB.extrinsic_nullspace(wt.structure,q)
N = N_in*N_ex
A*N
GM.activate!(;scalefactor=1);plot_traj!(wt;showmesh=false,showground=false)

# Contact Surfaces
halfspaces = RB.StaticContactSurfaces(
    [
        RB.HalfSpace([0,-1.0,0],[0,-0.0025e3,0]),
        RB.HalfSpace([0, 1.0,0],[0, 0.0025e3,0])
    ]
)

# Frictional Contact Dynamics

dt = 1e-4
tspan = (0.0,1.0)

prob = RB.DynamicsProblem(
    wt,
    halfspaces,
    RB.RestitutionFrictionCombined(
        RB.NewtonRestitution(),
        RB.CoulombFriction(),
    )
)

# Zhong06
RB.solve!(
    prob,
    RB.DynamicsSolver(
        RB.Zhong06(),
        # RB.MonolithicContactSolver(
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        )
    );
    dt,tspan,ftol=1e-10,maxiters=50,verbose=true,
    exception=true,progress=false,
    max_restart=1
)

b1g = RB.get_trajectory!(wt,1,0)
lines(b1g[3,:])

b1vg = RB.get_velocity!(wt,1,0)
lines(b1vg[3,:])

b1v1 = RB.get_velocity!(wt,1,1)
lines(b1v1[2,:])

GM.activate!(;scalefactor=1);plot_traj!(wt;showmesh=false,showground=false)

# Moreau
RB.solve!(
    prob,
    RB.DynamicsSolver(
        RB.Moreau(0.5),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        )
    );
    dt,tspan,ftol=1e-12,maxiters=50,verbose=true,
    exception=false,progress=false,
)

GM.activate!(;px_per_unit=1,scalefactor = 1);plot_traj!(wt;showmesh=false,showground=false)


b1g = RB.get_trajectory!(wt,1,0)
lines(b1g[3,:])

b1v1 = RB.get_velocity!(wt,1,1)
lines(b1v1[2,:])


me = RB.mechanical_energy!(wt;gravity=false)
Makie.lines(me.E)

using FileIO
using GLMakie
brain = load(assetpath("brain.stl"))

GLMakie.activate!(;scalefactor = 1) 
GM.activate!(;scalefactor = 2) #bug

fig = Figure()
ax1 = Axis3(fig[1, 1], aspect = :equal, title = "aspect = :equal")
ax2 = Axis3(fig[1, 2], aspect = :data, title = "aspect = :data")
for ax in [ax1, ax2]
    mesh!(ax, brain, color = :gray80)
end
DataInspector(fig)
fig

using ForwardDiff

ForwardDiff.gradient((x)->atan(x[2],x[1]),[0,1.0])
using FiniteDiff
FiniteDiff.finite_difference_gradient((x)->atan(x[2],x[1]),[0,1.0])

R = exp(RB.skew(rand(3)))

rR = rand(RotMatrix3)
rotation2angles(vec(R))
RotZYX(rR)

function myfun(x)
    [x[1]^2,x[2]^2,x[4]*x[2]]
end
x = rand(4)
angles_jacobian = ForwardDiff.jacobian(myfun,x)
angles_hessians = reshape(ForwardDiff.jacobian(x -> ForwardDiff.jacobian(myfun, x), x),3,4,4)

