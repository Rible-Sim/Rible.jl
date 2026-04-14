# # Rible.jl Quick Start

# Rible.jl is a multibody dynamics simulation platform for rigid-flexible robots.
# This page walks through a complete workflow — from modeling to solving — using a rolling superball tensegrity robot.

# ## Installation

# In the Julia REPL, run:

# ```julia
# import Pkg
# Pkg.activate(; temp=true)
# Pkg.add(url="https://github.com/Rible-Sim/Rible.jl")
# Pkg.add(url="https://github.com/Rible-Sim/Rible.jl", subdir="RibleTensegrity")
# ```

# > **Note**: Rible.jl requires Julia ≥ 1.11.

# ## Import

using Rible
import Rible as RB
using RibleTensegrity
import RibleTensegrity as RT
using CircularArrays
using StaticArrays
using Rotations
using LinearAlgebra
using TypeSortedCollections

# ## Build the Model

# Rible.jl ships reusable example definitions under `examples/`. The following loads the superball tensegrity robot together with its rigid-bar helper:

include(joinpath(pathof(Rible), "../../examples/bodies/rigidbar_nonsmooth_repro.jl"));
include(joinpath(pathof(Rible), "../../examples/robots/superball.jl"));

# Create a superball instance with moderate rolling motion and disable external mesh loading to keep the example self-contained:

l = 1.7 / 2
d = l / 2

bot = superball(0.0;
    origin_velocity = SVector(2.0, 1.0, 0.0),
    ω = SVector(0.0, 0.0, 1.0),
    μ = 0.05,
    e = 0.0,
    l,
    d,
    z0 = l^2 / (sqrt(5) * d) - 1e-3,
    loadmesh = false,
    visible = true,
    constrained = false,
)

# ## Define the Simulation Problem

# `DynamicsProblem` bundles the robot model, physical environment, and contact model into a complete simulation problem.

# Define the ground collision surface — a half-space with normal pointing in the positive z direction:

ground_plane = RB.StaticContactSurfaces([
    RB.HalfSpace([0.0, 0.0, 1.0], zeros(3))
]);

# Choose a contact model — Newton restitution + Coulomb friction:

contact_model = RB.RestitutionFrictionCombined(
    RB.NewtonRestitution(),
    RB.CoulombFriction(),
);

# Assemble the `DynamicsProblem`:

prob = RB.DynamicsProblem(bot; env=ground_plane, contact_model=contact_model)

# ## Solve

# Use the `Zhong06` integrator with an inner-layer contact solver, and run a short simulation:

solver = RB.DynamicsSolver(
    RB.Zhong06(),
    RB.InnerLayerContactSolver(RB.InteriorPointMethod()),
)

result = RB.solve!(prob, solver; tspan=(0.0, 1.0), dt=1e-2, ftol=1e-12, maxiters=200, exception=false);


# The simulation trajectory is stored inside `bot` and can be accessed via `result`.

# ## Visualize

# Use the trajectory plotting helper to inspect the rolling motion:

using CairoMakie
plot_traj!(bot; do_slide = false, show_loci = false, show_mesh = true, show_ground = false)

# ## Next Steps

# You have completed your first Rible.jl dynamics simulation. Continue with:

# - **Modeling basics**: [Bodies and Apparatuses](body_and_apparatus.md), [Structure](structure.md), and [Control Hub](hub.md)
# - **Dynamics**: [Robot Dynamics](dynamics.md) and [Errors, Costs, and Sensitivity](adjoint/index.md)
# - **Visualization**: [Visualization](vis.md)
