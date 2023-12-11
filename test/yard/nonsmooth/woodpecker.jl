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
    max_restart=3
)

b1g = RB.get_trajectory!(wt,1,0)
lines(b1g[3,:])

b1vg = RB.get_velocity!(wt,1,0)
b1vg = RB.get_mid_velocity!(wt,1,0)
lines(b1vg[3,:])

b1ω = RB.get_mid_angular_velocity!(wt,1)
b1θ = RB.get_mid_orientation!(wt,1)
lines(b1θ[3,:],b1ω[1,:])

b2ω = RB.get_mid_angular_velocity!(wt,2)
b2θ = RB.get_mid_orientation!(wt,2)
lines(b2θ[3,:],b2ω[1,:])

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
    dt,tspan,ftol=1e-10,maxiters=50,verbose=true,
    exception=false,progress=false,
)

GM.activate!(;px_per_unit=1,scalefactor = 1);plot_traj!(wt;showmesh=false,showground=false)


b1vg = RB.get_velocity!(wt,1,0)
lines(b1vg[3,:])

b1ω = RB.get_mid_angular_velocity!(wt,1)
lines(b1ω[1,:])

b1θ = RB.get_mid_orientation!(wt,1)
lines(b1θ[3,:],b1ω[1,:])

b2ω = RB.get_mid_angular_velocity!(wt,2)
lines(b2ω[1,:])

b2θ = RB.get_mid_orientation!(wt,2)
lines(b2θ[3,:],b2ω[1,:])



dts = [1e-3,5e-4,2e-4,1e-4,1e-5]
tspan = (0.0,0.08)
wt_dts = [
    begin
        wt = woodpecker()
        # Zhong06
        RB.solve!(
            RB.DynamicsProblem(
                wt,
                halfspaces,
                RB.RestitutionFrictionCombined(
                    RB.NewtonRestitution(),
                    RB.CoulombFriction(),
                )
            ),
            RB.DynamicsSolver(
                RB.Zhong06(),
                # RB.MonolithicContactSolver(
                RB.InnerLayerContactSolver(
                    RB.InteriorPointMethod()
                )
            );
            dt,tspan,ftol=1e-11,maxiters=50,verbose=false,
            exception=true,progress=true,
            max_restart=3
        ).prob.bot
    end
    for dt in dts
]


_,err_avg = RB.get_err_avg(wt_dts;bid=2,pid=2,di=2,field=:traj)
fig = Figure()
ax = Axis(fig[1,1])
plot_convergence_order!(ax,dts[begin:end-1],err_avg)
ax2 = Axis(fig[1,2])
for bot in wt_dts
    b1r2 = RB.get_trajectory!(bot,2,2)
    lines!(ax2,bot.traj.t,b1r2[2,:])
end
fig