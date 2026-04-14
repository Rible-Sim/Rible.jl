# Errors, Costs, and Objectives

```@meta
CurrentModule = Rible
ShareDefaultModule = true
```

```@setup
using Rible
import Rible as RB
using StaticArrays
using LinearAlgebra
import TypeSortedCollections as TSC
```

In trajectory optimization (e.g., iLQR), and reinforcement learning training, we need a quantitative scalar metric to evaluate a system's performance. Rible uses the `AbstractObjective` framework to aggregate scattered physical observations (Gauges) and execution costs (Actions) into a unified mathematical objective.

## 1. From Error to Cost Chain

In Rible, the **Error** is a local scalar defined at the gauge level, while the **Cost** is a weighted global indicator defined at the system level.

### 1.1 ErrorGauge: The Atomic Unit of Evaluation

As mentioned in the [Control Hub](../hub.md), `ErrorGauge` is responsible for computing the deviation in a single dimension. For example:

- "The Euclidean distance from the end effector to the target"
- "The angle of the rod deviating from the vertical"

Each `ErrorGauge` yields a differentiable scalar error.

### 1.2 Objective: The Baton of Weights

The core role of the `Objective` is to perform a weighted sum. It defines which errors are more important and how much control effort is required.

The total cost ``J`` typically consists of two parts:

```math
J = \underbrace{\sum_{k=0}^{N-1} L(x_k, u_k)}_{\text{Trajectory Cost}} + \underbrace{\Phi(x_N)}_{\text{Terminal Cost}} 
```

Where the instantaneous cost ``L`` is the weighted sum of all active `ErrorGauge`s and `Actuator`s.

## 2. Structure of the Objective Function (Objective)

A typical `Objective` object contains the following key pieces of information:

- **Trajectory Cost Weights**: weights computed at each time step.
- **Terminal Cost Weights**: weights computed only at the final time step. Typically used to emphasize achieving the “final goal” (e.g., must stop at a specific position).
- **Actuator Weights**: penalties on the magnitude of the control input ``u`` (usually related to energy consumption).
- **Time Scaling Factor (``\eta``)**: used to handle weight adjustments under variable step sizes or specific integration schemes.

### Core Interface: Cost Computation

Evaluations based on the `Objective` are implemented via the following functions:

- **`cost!(bot, objt, inst_state, u; mode=:trajectory)`**: computes a single-step instantaneous cost.
- **`cost!(bot, objt, solver, dt)`**: computes the total cost for the entire trajectory.

## 3. Sensitivity Analysis: Differentiation Interfaces

For trajectory optimization and gradient-based algorithms (such as adjoint sensitivity analysis), cost values alone are not enough. Rible provides efficient, cache-free differentiation interfaces:

### Gradient

**`cost_gradient!`**: computes the gradient of the cost with respect to the following variables:

- position vector ``q`` and momentum vector ``p``.
- control input vector ``u``.
- system parameters ``c`` or policy parameters ``\theta``.

### Hessian

**`cost_hessian!`**: computes the second-order partial derivatives of the cost with respect to state and control. This is used in algorithms like iLQR to build a local quadratic approximation.

> [!TIP]
> These differentiation interfaces internally leverage `measure_jacobian!`. Since the gauge-level computations are already highly efficient, differentiating at the Objective level adds virtually no extra overhead.

## 4. Example: Documenter-Ready Objective and cost! Workflow

The following examples all use `@example` and can be executed inline by Documenter.

### Step 1: Build a Minimal Runnable `Structure`

The goal of this step is to prepare a minimal physical system that can be evaluated by `cost!`. We will sequentially create `Locus`, `RigidBodyProperty`, `RigidBodyState`, and `NCF.NC1P3V`, then assemble them into a `Structure` and perform a single `update!`. This ensures subsequent gauge and objective function examples focus on the cost calculation logic, without being disturbed by detailed model specifics.

```@example
mass_locus = RB.Locus(SVector(0.0, 0.0, 0.0))
loci = [RB.Locus(SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0))]
inertia = SMatrix{3,3,Float64}(I)

prop = RB.RigidBodyProperty(1, true, 1.0, inertia, mass_locus, loci; visible = false)
r0 = SVector(0.0, 0.0, 0.0)
R0 = SMatrix{3,3,Float64}(I)
state = RB.RigidBodyState(prop, r0, R0, zero(r0), zero(r0))
coords = RB.NCF.NC1P3V(r0, r0, R0)
body = RB.RigidBody(prop, state, coords, nothing)
body.prop.id, length(body.prop.loci)
```

The previous block completes the single-body modeling; next we wrap it into a container and initialize the structure state.

```@example
bodies = TSC.TypeSortedCollection([body])
apparatuses = TSC.TypeSortedCollection(Int[])
st = RB.Structure(bodies, apparatuses, RB.Connectivity(bodies, apparatuses))
RB.update!(st)
n_q = RB.get_num_of_coords(st)
n_c = RB.get_num_of_cstr(st)
n_q, n_c
```

### Step 2: Define an Error Gauge and assemble `ControlHub`

`ErrorGauge` maps the “deviation of the current state from the reference” into a differentiable scalar; `Signifier` specifies which body and which locus to measure. Then, via `Coalition` and `ControlHub`, the indices of error and action signals in the system vectors are established, and a `Robot` is constructed as a single entry point.

```@example
error_gauge = RB.ErrorGauge(
    1,
    RB.Signifier(body, 1),
    RB.PositionCaptum(),
    [0.2, 0.0, 0.0],
)

capta_gauges = TSC.TypeSortedCollection(Int[])
error_gauges = TSC.TypeSortedCollection([error_gauge])
actuators = TSC.TypeSortedCollection(Int[])
num_error_gauges = length(error_gauges)
error_gauge.id, num_error_gauges
```

Next bind the gauge, actuators, and structure into a control hub and generate a `Robot` ready for `cost!`.

```@example
coalition = RB.Coalition(st, capta_gauges, error_gauges, actuators)
hub = RB.ControlHub(st, capta_gauges, error_gauges, actuators, coalition)
bot = RB.Robot(st, hub)
num_errors = length(error_gauges)
num_actions = RB.get_num_of_actions(bot)
num_errors, num_actions
```

### Step 3: Construct an `Objective` and Inspect Time Weights

The four weight vectors here correspond to: trajectory error terms, trajectory control terms, terminal error terms, and terminal control terms. Because the semantic meaning of these constructor arguments is not self-evident, the example keeps the minimal necessary comments to indicate the role of each vector.

```@example
objt = RB.Objective(
    [10.0],      # trajectory_error_gauges_weights
    Float64[],   # trajectory_actuators_weights
    [100.0],     # terminal_error_gauges_weights
    Float64[],   # terminal_actuators_weights
    t -> 1 / (1 + t),
)

ηs = RB.get_trajectory_cost_weights(objt, [0.0, 0.1, 0.2], 0.1)
ηs
```

### Step 4: Compute Single-Step Trajectory / Terminal Cost

We directly call the single-step interface `cost!(bot, objt, inst_state, u; mode=...)`. `mode=:trajectory` uses the running weights, and `mode=:terminal` uses the terminal weights; the two costs can be displayed side-by-side to visually compare the same state under the two weight sets.

```@example
inst = bot.structure.state.system
u = Float64[]

ϕ_traj = RB.cost!(bot, objt, inst, u; mode = :trajectory)
ϕ_term = RB.cost!(bot, objt, inst, u; mode = :terminal)
ϕ_traj, ϕ_term
```

## 5. Summary

By decoupling what we measure (Gauge) from how much weight it carries (Objective), Rible enables developers to:

1. **Rapid parameter tuning**: simply adjust the weight arrays in the `Objective` to change the robot's behavior style.
2. **Multi-target switching**: easily switch between "speed-first" and "smoothness-first" task goals on the same physical model.
3. **Support for advanced algorithms**: leverage built-in `cost_gradient!` and `cost_hessian!` to interface with high-performance optimizers.