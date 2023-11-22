figdir::String = ""
if Sys.iswindows()
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\FrictionalContact\CMAME"
elseif Sys.isapple()
    figdir::String = raw"/Users/jacob/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/Papers/FrictionalContact/CMAME"
end
#-- end

include("deps.jl")
using AbbreviatedStackTraces #jl
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true #jl
ENV["JULIA_STACKTRACE_MINIMAL"] = true #jl
import Rible as RB
include("../../vis.jl")
includet("../../vis.jl") #jl

#--  slider crank 
include("../../../examples/robots/slider_crank.jl")
includet("../../../examples/robots/slider_crank.jl")#jl

# natural coordinates
sc = slider_crank(;coordsType=RB.NCF.NC)

# Quaternion coordinates
sc = slider_crank(;coordsType=RB.QCF.QC)

# Contact Surfaces
planes = RB.StaticContactSurfaces(
    [
        RB.Plane([0,0, 1.0],[0,0,-0.026]),
        RB.Plane([0,0,-1.0],[0,0, 0.026])
    ]
)

plot_traj!(sc;showground=false)

RB.has_constant_mass_matrix(sc)

dt = 1e-3
tspan = (0.0,683dt)

# No Contact Dynamics
prob = RB.DynamicsProblem(sc,)

RB.solve!(
    prob,
    RB.DynamicsSolver(RB.Zhong06());
    dt,tspan,ftol=1e-14,maxiters=50,verbose=true,exception=true,progress=false,
)

plot_traj!(sc;showground=false)

# Frictionless Contact Dynamics

prob = RB.DynamicsProblem(
    sc,
    planes,
    RB.RestitutionFrictionCombined(
        RB.NewtonRestitution(),
        RB.Frictionless(),
    )
)

# two-layer
RB.solve!(
    prob,
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        )
    );
    dt,tspan,ftol=1e-10,maxiters=50,verbose=true,exception=true,progress=false,
)

# Mono
RB.solve!(
    prob,
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.MonolithicContactSolver(
            RB.InteriorPointMethod()
        )
    );
    dt,tspan,ftol=1e-11,maxiters=50,verbose=true,exception=true,progress=false,
)

# Frictional Contact Dynamics
prob = RB.DynamicsProblem(
    sc,
    planes,
    RB.RestitutionFrictionCombined(
        RB.NewtonRestitution(),
        RB.CoulombFriction(),
    )
)

# two-layer
RB.solve!(
    prob,
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        )
    );
    dt,tspan,ftol=1e-14,maxiters=50,verbose=true,exception=true,progress=false,
)

plot_traj!(sc;showground=false)

me = RB.mechanical_energy!(sc)
lines(me.E)
lines(me.V)
lines(me.T)

using PreallocationTools
using ForwardDiff


