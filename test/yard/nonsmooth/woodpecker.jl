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

GM.activate!(;scalefactor=1);plot_traj!(wt;showmesh=false,showground=false)

# Contact Surfaces
halfspaces = RB.StaticContactSurfaces(
    [
        RB.HalfSpace([0,-1.0,0],[0,-0.0025,0]),
        RB.HalfSpace([0, 1.0,0],[0, 0.0025,0])
    ]
)

# Frictional Contact Dynamics

dt = 1e-3
tspan = (0.0,0.2)

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
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        )
    );
    dt,tspan,ftol=1e-14,maxiters=50,verbose_contact=true,exception=true,progress=false,
)

GM.activate!(;scalefactor=1);plot_traj!(wt;showmesh=false,showground=false)

# two-layer
RB.solve!(
    prob,
    RB.DynamicsSolver(
        RB.Moreau(1.0),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        )
    );
    dt,tspan,ftol=1e-14,maxiters=50,verbose=true,exception=true,progress=false,
)

GM.activate!(;px_per_unit=2,scalefactor = 2);plot_traj!(sc;showground=false)

me = RB.mechanical_energy!(sc)
Makie.lines(me.E)
