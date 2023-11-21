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

plot_traj!(sc;showground=false)

RB.has_constant_mass_matrix(sc)

dt = 1e-3
tspan = (0.0,1.0)

# No Contact 
prob = RB.DynamicsProblem(sc,)

RB.solve!(
    prob,
    RB.DynamicsSolver(RB.Zhong06());
    dt,tspan,ftol=1e-14,maxiters=50,verbose=true,exception=true,progress=false,
)

plot_traj!(sc;showground=false)

# Contact Surfaces
planes = RB.StaticContactSurfaces(
    [
        RB.Plane([0,0, 1.0],[0,0,-0.026]),
        RB.Plane([0,0,-1.0],[0,0, 0.026])
    ]
)


# Quaternion coordinates
sc = slider_crank(;coordsType=RB.QCF.QC)

# Frictionless Contact Dynamics

prob = RB.DynamicsProblem(
    sc,
    planes,
    RB.RestitutionFrictionCombined(
        RB.NewtonRestitution(),
        RB.Frictionless(),
    )
)

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

RB.solve!(
    prob,
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.MonolithicContactSolver(
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

η = @SVector rand(3)
q = normalize(@SVector rand(4))
ω = @SVector rand(3)
q̇ = RB.QCF.Lᵀmat(q)*ω
RB.QCF.∂Rη∂q(q,η)*q̇
-2RB.QCF.Rmat(q)*RB.skew(η)*RB.QCF.Lmat(q)*q̇



H = RB.QCF.∂²Rη∂qᵀ∂q(η)