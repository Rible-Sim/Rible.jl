using Revise #jl
import Rible as RB
include(joinpath(pathof(RB),"../../yard/adjoint.jl"))

using AbbreviatedStackTraces #jl
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true #jl
ENV["JULIA_STACKTRACE_MINIMAL"] = true #jl

figdir::String = joinpath(pathof(RB),"../../tmp")
include(joinpath(pathof(RB),"../../test/vis.jl"))
includet(joinpath(pathof(RB),"../../test/vis.jl")) #jl
tw = 455.8843 #pt |> pt2px
scalefactor = 4


#--  cube
include(joinpath(pathof(RB),"../../examples/robots/cube.jl"))
includet(joinpath(pathof(RB),"../../examples/robots/cube.jl"))#jl

include(joinpath(pathof(RB),"../../examples/robots/pointmass.jl"));
includet(joinpath(pathof(RB),"../../examples/robots/pointmass.jl")) #jl

#-- cube
cubebot = make_cube()
plot_traj!(cubebot;showwire=true)

plane_env = RB.StaticContactSurfaces(
    [
        RB.HalfSpace([0,0,1.0],[0,0,0.0])
    ]
)

ros = [[0,0,z] for z = 10.0:-0.1:0.0]
qs = [
    begin
        cube = make_cube(origin_position,RotXY(π/4,π/4))
        q = RB.get_coords(cube.structure)
    end
    for origin_position in ros
]
# p = RB.generalized_alpha(1.0)
Ro = [
    0.368869  -0.824063   0.429949;
    0.769524   0.011314  -0.638518;
    0.521314   0.566385   0.63831
]
origin_velocity = [0.0,0.0,0.0]
ωo = [3.0,4.0,5.0]
# cube = make_cube()
# cube = make_cube(ros[1])
cube = make_cube(ros[2],Ro,origin_velocity,0.1.*ωo)
tspan = (0.0,1.91)
tspan = (0.0,5.0)
dt = 1e-3

prob = RB.DynamicsProblem(
    cube,
    plane_env,
    RB.RestitutionFrictionCombined(
        RB.NewtonRestitution(),
        RB.CoulombFriction(),
    )
)
RB.solve!(prob,
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        )
    );
    tspan,dt,ftol=1e-10,maxiters=50,exception=false
)
RB.solve!(prob,
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.MonolithicContactSolver(
            RB.InteriorPointMethod()
        )
    );
    tspan,dt,ftol=1e-10,maxiters=50,exception=false
)

plot_traj!(cube;
    xlims = [-10,10],
    ylims = [-10,10],
    zlims = [-1e-3,10],
    showinfo=false,
    showwire=true
)

ME = RB.mechanical_energy!(cube)
lines(cube.traj.t,ME.E)



# ground plane
θ=0.0
ground_plane = RB.StaticContactSurfaces(
    [
        RB.HalfSpace([-tan(θ),0,1],zeros(3)),
    ]
);

# time
tspan = (0.0,1.5)
h = 1e-3;

# parameters and initial conditions
restitution_coefficients = [0.5]
v0s = [1.0];

# pointmass
pm = new_pointmass(;
    e = restitution_coefficients[1],
    μ=0.1,
    origin_velocity = [v0s[1],0,0]
);

# frictional contact dynamics problem 
prob = RB.DynamicsProblem(
    pm,
    ground_plane,
    RB.RestitutionFrictionCombined(
        RB.NewtonRestitution(),
        RB.CoulombFriction(),
    )
);

# simulate
RB.solve!(
    prob,
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.MonolithicContactSolver(
            RB.InteriorPointMethod()
        )
    );
    tspan,dt=h,
    ftol=1e-14,
    maxiters=50,
    exception=false
);

GM.activate!(;scalefactor=2); plot_traj!(
    pm,
    xlims=(-1e-3,1.0),
    ylims=(-0.4,0.4),
    zlims=(-1e-3,1.0),
    showinfo =false,
    showmesh=false,
    showwire=false,
    showlabels=false,
    showcables=false,
    showpoints=true,
)