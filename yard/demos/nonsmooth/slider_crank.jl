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

plot_traj!(sc;showground=false)

RB.has_constant_mass_matrix(sc)

dt = 1e-3
tspan = (0.0,1.0)

# No Contact Dynamics
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

q = RB.get_coords(sc.structure)
RB.make_cstr_jacobian(sc.structure)(q)

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
using PreallocationTools
using ForwardDiff

randmat = rand(5, 3)
sto = similar(randmat)
stod = DiffCache(sto)

sto = get_tmp(sto, 0.0)

function claytonsample!(sto, τ, α; randmat = randmat)
    sto = get_tmp(sto, τ)
    sto .= randmat
    τ == 0 && return sto

    n = size(sto, 1)
    for i in 1:n
        v = sto[i, 2]
        u = sto[i, 1]
        sto[i, 1] = (1 - u^(-τ) + u^(-τ) * v^(-(τ / (1 + τ))))^(-1 / τ) * α
        sto[i, 2] = (1 - u^(-τ) + u^(-τ) * v^(-(τ / (1 + τ))))^(-1 / τ)
    end
    return sto
end

ForwardDiff.derivative(τ -> claytonsample!(stod, τ, 0.0), 0.3)
ForwardDiff.jacobian(x -> claytonsample!(stod, x[1], x[2]), [0.3; 0.0])
ForwardDiff.hessian(τ -> claytonsample!(stod, τ, 0.0), [0.3])
ForwardDiff.derivative(τ -> claytonsample!(stod, τ, 0.0), 0.3)

