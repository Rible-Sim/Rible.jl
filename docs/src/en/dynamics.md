```@meta
CurrentModule = Rible
ShareDefaultModule = true
```

```@setup
using Rible
import Rible as RB
using StaticArrays
using LinearAlgebra
using Rotations
using TypeSortedCollections
using Makie
import CairoMakie
import GeometryBasics as GB
import RungeKutta

CairoMakie.activate!()
Makie.inline!(true)
const Point3 = GB.Point3
```

# Robot, Dynamics, and Simulation

## Robot

The `Robot` struct is the digital twin of the physical system. It stores all static properties and dynamic history of the robot.

Its core components include:

- **`structure`** — maintains the multibody topology (see [Structure](@ref)), containing all rigid bodies, joints, and apparatuses.
- **`hub`** — manages control-related actuators and error gauges in one place (see [Control Hub](@ref)).
- **Trajectory containers** — persist all data generated during simulation, including `traj` (state trajectory), `contacts_traj` (contact-related trajectory), and `control_traj` (control-related trajectory).

---

## Dynamics and Simulation

Rible provides a extensible dynamics simulation framework supporting multiple integration algorithms, contact mechanics, and control strategies.

> **Core Architecture**: In Rible, the **physical model** (`DynamicsProblem`) and the **solution strategy** (`DynamicsSolver`) are decoupled. You can seamlessly swap between different solver algorithms for the same physical system without modifying any model definition.

## Core Components

### DynamicsProblem

`DynamicsProblem` is the data container for simulation tasks. It is an immutable struct that binds all static properties of the physical system and its environment, providing clear computation boundaries for the solver.

It contains:

- **`Robot`** — the mechanical system's topology, mass distribution, and runtime state (see [Robot](@ref)).
- **`Objective`** — the cost function for adjoint analysis, typically built from `ErrorGauge` instances.
- **`Environment`** — provides spatial context (`Geometry`) and passive external forcing (`Field`), avoiding hard-coding global forces into the robot. Default to (`EmptyEnv`).
  - **`Field`** — provides global external forces, e.g., gravity field (`Gravity`).
  - **`Geometry`** — provides surfaces or mesh terrain for collision detection.
- **`Policy`** — maps the current system state ``(t, q, v)`` to active control forces `u`. Encapsulating it separately enables free switching between open-loop and closed-loop control. Default to `NoPolicy`.
  - **Open-loop**: e.g., time-based curves.
  - **Closed-loop feedback**: e.g., proportional-derivative (`DiscretePD`) or neural network policies.
- **`Models`** — constrains the mathematical form of physical interactions, e.g., contact models (subtypes of `AbstractContactModel`).

### DynamicsSolver

If `DynamicsProblem` defines "what to compute," then `DynamicsSolver` handles "how to compute." Due to the complexity of rigid-flexible coupled systems, it is designed as a pluggable dispatch center. It does not store physical state itself — instead, it composes different numerical integrators and sub-solvers to specify the behavior of the system's evolution.

- **`Integrator`** — To specify the scheme which computes the next system state ``(q_{k+1}, v_{k+1})`` at discrete time steps, handling stiffness and algebraic constraints.
- **Sub-Solvers** — includes the body solver, apparatus solver, and contact solver. For contact-intensive scenarios, advanced nonlinear complementarity solvers can be configured independently.

---

## Configuring Integrators

Integrators are the numerical engines of the `DynamicsSolver`. Selecting the right integrator is crucial for simulation stability and accuracy.

### Zhong 06
The default recommended integrator for most constrained multibody systems, originally proposed by Prof. Wanxie Zhong, and has been modified and developed for non-smooth contact dynamics [luoNonsmoothModifiedSymplectic2024](@cite). It is a second-order, constraint-preserving, and symplectic-like scheme that maintains energy conservation and constraint manifold consistency over long durations.

```@example
# Construct with default settings
integrator = Zhong06()
solver = DynamicsSolver(integrator)
```

### Runge-Kutta (RK) Family
Supports various generic implicit Runge-Kutta tableaus, like Gauss-Legendre (GL), via the `RungeKutta.jl` package. It is excellent for high-accuracy requirements in smooth systems.

```@example
using RungeKutta
# Construct a 2-stage Gauss-Legendre (implicit) integrator
integrator = RKIntegrator(RungeKutta.TableauGauss(2))
solver = DynamicsSolver(integrator)
```

### Extra Integrators (Generalized-α, Moreau-Jean)
For specialized needs in comparison and evaluation, Rible provides additional integrators in the [RibleExtraIntegrators](@ref "Generalized Alpha and Moreau-Jean Integrators") package.

---

## Typical Simulation Workflow

Rible provides a two-level interface that balances out-of-the-box convenience with deep customization.

### 1. Convenient API: `solve!`

For most standard scenarios, use the `solve!` method directly. It handles memory allocation and runs the main loop internally.

The following example reuses the repository's spinning-top model. The top starts with high angular velocity, touches a planar ground, and produces a trajectory that we can inspect at the end of the page.

```@example
include(joinpath(pathof(Rible), "../../examples/robots/spinningtop.jl"))

origin_position = [0.0, 0.0, 0.5]
R = RotX(0.0)
origin_velocity = [1.0, 0.0, 0.0]
Ω = [0.0, 0.0, 200.0]

top = make_top(origin_position, R, origin_velocity, Ω, :NCF; μ = 0.95, e = 0.5, loadmesh = false)
top.structure.connectivity.num_of_bodies
```

The model already bundles its `Structure` and `ControlHub`, so the next step is to define the contact surface and interaction law.

```@example
planes = StaticContactSurfaces([
    HalfSpace([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
])

contact_model = RestitutionFrictionCombined(
    NewtonRestitution(),
    CoulombFriction(),
)

prob = DynamicsProblem(top; env=planes, contact_model)
```

With the setup complete, define a contact-aware solver and run a short simulation.

```@example
solver = DynamicsSolver(
    Zhong06(),
    InnerLayerContactSolver(InteriorPointMethod()),
)

sim_result = solve!(
    prob,
    solver;
    tspan = (0.0, 0.2),
    dt = 5e-4,
    ftol = 1e-14,
    maxiters = 50,
    exception = false,
)

tip_position = get_trajectory!(top, 1, 1).u[end]
tip_velocity = get_velocity!(top, 1, 1).u[end]
tip_position, tip_velocity
```

### 2. Simulation Container — Simulator

`Simulator` is the execution container responsible for maintaining runtime state.

When you need fine-grained control over the main loop (e.g., integrating with an external reinforcement learning framework, adding runtime controllers, or implementing hot-restart from exceptions), explicitly create a `Simulator` and step through manually:

```@example
top_restart = make_top(origin_position, R, origin_velocity, Ω, :NCF; μ = 0.95, e = 0.5, loadmesh = false)
prob_restart = DynamicsProblem(top_restart; env=planes, contact_model)

sim = Simulator(prob_restart, solver; tspan = (0.0, 0.02), dt = 5e-4, restart = true)
solve!(sim, solver; dt = 5e-4)
traj_len_sim = length(top_restart.traj)
t_end_sim = top_restart.traj.t[end]
traj_len_sim, t_end_sim
```

After integration, extract post-processing metrics directly from the top's trajectory containers.

```@example
E = mechanical_energy!(top)
E[begin], E[end]
```

To visualize the path directly, extract one point trajectory from the recorded simulation and draw it with `Makie.lines`.

```@example
tip_traj = get_trajectory!(top, 1, 1).u
tip_points = Point3.(tip_traj)

fig_traj = Figure()
ax_traj = Axis3(fig_traj[1, 1], aspect = :data)
Makie.lines!(ax_traj, tip_points, linewidth = 2)
fig_traj
```

Additionally, trajectory data can be fed into the [Visualization](vis.md) workflows for further plotting and animation.

---

## Core Mechanisms and Performance

Rible achieves high performance through zero-overhead abstractions enabled by multiple dispatch and the traits system.

1. **Trajectory memory pre-allocation (`prepare_traj!`)**
   Called automatically when initializing a `Simulator`. Pre-allocates a contiguous memory block for the state sequence over the entire time span.

2. **Automatic workspace generation (`generate_cache`)**
   Before simulation starts, generates a type-specific `Workspace` at compile time based on the concrete problem and solver types, to allow algorithm safely owns its cache matrices (e.g., Jacobians).

3. **Trait-based dispatch optimization**
    Uses static analysis of physical traits, like `has_constant_mass_matrix(bot)`, to check whether the mass matrix is constant. If so, the compiler automatically skips mass matrix recomputation.

## References

```@bibliography
Pages = ["dynamics.md"]
```

