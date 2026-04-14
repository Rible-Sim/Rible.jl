# Errors, Costs, and Sensitivity Analysis

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

> [!WARNING]
> Sensitivity analysis support is under development. API patterns described in this page may still change across solver families and future releases.

In trajectory optimization (e.g., iLQR) and reinforcement learning, we need a scalar metric to evaluate system performance. Rible uses the `AbstractObjective` framework to aggregate physical observations (Gauges) and execution costs (Actions) into a unified mathematical objective. Sensitivity analysis then computes gradients of this objective with respect to parameters, enabling gradient-based optimization.

## From Error to Cost

In Rible, the **Error** is a local scalar defined at the gauge level, while the **Cost** is a weighted global indicator defined at the system level.

### Error Gauge

As mentioned in the [Control Hub](../hub.md), `ErrorGauge` computes the deviation in a single dimension — e.g., the Euclidean distance from an end effector to a target, or the angle of a rod deviating from vertical. Each `ErrorGauge` yields a differentiable scalar error.

### Objective

`Objective` performs weighted summation over all active gauges and actuators. The total cost ``J`` consists of two parts:

```math
J = \underbrace{\sum_{k=0}^{N-1} L(x_k, u_k)}_{\text{Trajectory Cost}} + \underbrace{\Phi(x_N)}_{\text{Terminal Cost}}
```

A typical `Objective` object contains:

- **Trajectory error weights** — applied at every time step.
- **Terminal error weights** — applied only at the final step, emphasizing goal achievement.
- **Actuator weights** — penalties on control input magnitude ``u`` (energy cost).
- **Time scaling factor** (``\eta``) — handles weight adjustments under variable step sizes.

The core evaluation interface is:

- **`cost!(bot, objt, inst_state, u; mode=:trajectory)`** — single-step instantaneous cost.
- **`cost!(bot, objt, solver, dt)`** — total cost over the entire trajectory.

### Differentiation Interface

For gradient-based algorithms, Rible provides buffer-free differentiation:

- **`cost_gradient!`** — gradient of the cost with respect to ``q``, ``p``, ``u``, and parameters ``\theta``.
- **`cost_hessian!`** — second-order derivatives for local quadratic approximation (used by iLQR).

> [!TIP]
> These interfaces internally leverage `measure_jacobian!`. 

## Sensitivity Analysis Overview

Sensitivity analysis measures how system output responds to changes in inputs and parameters. Rible provides two methods, both compatible with the Zhong06 discrete integration framework:


| Feature | Adjoint | Direct |
| :--- | :--- | :--- |
| Computation direction | Reverse Mode | Forward Mode |
| Best for | Many parameters, few outputs | Few parameters, high state dimension |
| Complexity | ``O(N_{obj})``, weakly dependent on parameter count | ``O(N_{param})``, linear in parameter count |
| Memory | High (stores forward trajectory) | Low (online progression) |
| Typical use | Large-scale optimization, policy training | Real-time control, online monitoring |

- If parameter dimension is high (e.g., neural networks) but the objective is scalar, prefer the **Adjoint method**.
- If parameters are few and online sensitivity trajectories are needed, prefer the **Direct method**.

## Adjoint Sensitivity Analysis

The adjoint method constructs backward-propagating variables and obtains the gradient of the objective with respect to many parameters in a single backward sweep.

### Theory

Consider a discrete-time system:

```math
x_{k+1} = f(x_k, u_k, \theta), \quad k = 0, \dots, N-1
```

with objective:

```math
J = \sum_{k=0}^{N-1} L_k(x_k, u_k) + \Phi(x_N)
```

Introduce adjoint variable ``\lambda_k`` with backward recurrence:

```math
\lambda_k = \frac{\partial L_k}{\partial x_k} + \left(\frac{\partial f}{\partial x_k}\right)^T \lambda_{k+1}
```

The parameter gradient is then:

```math
\frac{dJ}{d\theta} = \sum_{k=0}^{N-1} \left(\frac{\partial L_k}{\partial \theta} + \lambda_{k+1}^T \frac{\partial f}{\partial \theta}\right)
```

## Direct Sensitivity Analysis

The direct method advances both the state and sensitivity equations forward in time.

### Theory

Define the sensitivity matrix ``S = \partial x / \partial \theta``. The sensitivity equation is:

```math
\dot{S} = \frac{\partial f}{\partial x} S + \frac{\partial f}{\partial \theta}
```

This integrates forward concurrently with the original equations, yielding full trajectories ``x(t)`` and ``S(t)``.

## Sensitivity Analysis API

This page focuses on the public workflow first: evaluate the objective, construct the sensitivity solver, run `solve!`, then inspect `sim.solver_cache`. Cache type names are solver-family internals and should only matter if you are extending the implementation.

The sensitivity solvers are thin wrappers around an existing dynamics solver and an objective.

- **Adjoint** — `adj_solver = DiscreteAdjointDynamicsSolver(dyn_solver, objt)` followed by `adj_sen_solver = AdjointDynamicsSensitivitySolver(dyn_solver, adj_solver)`.
- **Direct** — `drc_sen_solver = DirectDynamicsSensitivitySolver(dyn_solver, objt)`.
- **Execution** — `solve!(prob, sensitivity_solver; ...)` returns a simulation result whose sensitivity outputs live in `sim.solver_cache`.

This keeps the public entry point consistent with the rest of the dynamics API: the model stays in `DynamicsProblem`, the numerical method stays in `DynamicsSolver`, and the sensitivity method is added as an outer wrapper.

After the solve finishes, inspect `sim.solver_cache` instead of relying on solver-family-specific cache type names.

- **Adjoint gradients** — `∂J∂x₀ᵀ`, `∂J∂θᵀ`, `∂J∂cᵀ`
- **Direct Jacobians** — `Jac_state`, `Jac_action`, `Jac_control_params`

The adjoint fields summarize how the scalar objective changes with respect to the initial state and parameter groups. The direct fields describe how the propagated state changes with respect to initial state, actions, and control parameters.


## Example: Cost and Sensitivity Workflow

### Step 1: Build a Minimal System

Create a single-body structure with an error gauge and control hub.

```@example
using Rible
using StaticArrays
using LinearAlgebra
import TypeSortedCollections as TSC

mass_locus = Locus(SVector(0.0, 0.0, 0.0))
loci = [Locus(SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0))]
inertia = SMatrix{3,3,Float64}(I)

prop = RigidBodyProperty(1, true, 1.0, inertia, mass_locus, loci; visible = false)
r0 = SVector(0.0, 0.0, 0.0)
R0 = SMatrix{3,3,Float64}(I)
state = RigidBodyState(prop, r0, R0, zero(r0), zero(r0))
coords = NCF.NC1P3V(r0, r0, R0)
body = RigidBody(prop, state, coords, nothing)

bodies = TSC.TypeSortedCollection([body])
apparatuses = TSC.TypeSortedCollection(Int[])
st = Structure(bodies, apparatuses, Connectivity(bodies, apparatuses))
update!(st)
```

### Step 2: Define Gauge, Hub, and Objective

Wire up an error gauge, add an external force actuator, and construct the objective.

```@example
error_gauge = ErrorGauge(
    1,
    Signifier(body, 1),
    PositionCaptum(),
    [0.1, 0.0, 0.0],
)

capta_gauges = TSC.TypeSortedCollection(Int[])
error_gauges = TSC.TypeSortedCollection([error_gauge])
force_actuator = ExternalForceActuator(
    1,
    Signifier(body, 1),
    NaiveOperator(2),
    [1.0 0.0; 0.0 1.0; 0.0 0.0],
    [0.0, 0.0],
)
actuators = TSC.TypeSortedCollection([force_actuator])

coalition = Coalition(st, capta_gauges, error_gauges, actuators)
hub = ControlHub(st, capta_gauges, error_gauges, actuators, coalition)
bot = Robot(st, hub)
prob = DynamicsProblem(bot; env=GravityEnv())

objt = Objective(
    [10.0],      # trajectory_error_gauges_weights
    [1.0],       # trajectory_actuators_weights
    [100.0],     # terminal_error_gauges_weights
    [5.0],       # terminal_actuators_weights
    t -> 1 / (1 + t),
)
```

### Step 3: Compute Cost

```@example
inst = bot.structure.state.system
u = [0.2, -0.1]

ϕ_traj = cost!(bot, objt, inst, u; mode = :trajectory)
ϕ_term = cost!(bot, objt, inst, u; mode = :terminal)
ϕ_traj, ϕ_term
```

### Step 4: Adjoint Sensitivity Solve

```@example
dyn_solver = DynamicsSolver(Zhong06())
adj_solver = DiscreteAdjointDynamicsSolver(dyn_solver, objt)
adj_sen_solver = AdjointDynamicsSensitivitySolver(dyn_solver, adj_solver)

adj_sim = solve!(prob, adj_sen_solver; tspan = (0.0, 0.01), dt = 1e-3)
```

### Step 5: Direct Sensitivity Solve

Both methods reuse the same `dyn_solver` and `objt` — only the outer solver type changes:

```@example
drc_sen_solver = DirectDynamicsSensitivitySolver(dyn_solver, objt)
drc_sim = solve!(prob, drc_sen_solver; tspan = (0.0, 0.01), dt = 1e-3)
```


## Extract Gradients and Jacobians

The examples below show how to read results directly from `sim.solver_cache`.

### Adjoint Result

```@example
∂J∂x0 = adj_sim.solver_cache.∂J∂x₀ᵀ
∂J∂θ = adj_sim.solver_cache.∂J∂θᵀ
∂J∂c = adj_sim.solver_cache.∂J∂cᵀ

∂J∂x0, ∂J∂θ, ∂J∂c
```

`∂J∂x0` is the gradient with respect to the initial state, while `∂J∂θ` and `∂J∂c` collect the sensitivity trajectories with respect to control and structural parameters.

### Direct Result

The direct solver exposes the state Jacobian blocks directly from `drc_sim.solver_cache`:

```@example

Jac_state = drc_sim.solver_cache.Jac_state
Jac_action = drc_sim.solver_cache.Jac_action
Jac_control_params = drc_sim.solver_cache.Jac_control_params

state_sample = Jac_state[end][1:3, 1:3]
action_sample = Jac_action[end][1:3, :]
control_param_total = sum(sum(abs, block) for block in Jac_control_params)

(; state_sample, action_sample, control_param_total)
```

These Jacobians describe sensitivity of the propagated state with respect to the initial state, actions, and control parameters. In this example, `control_param_total` stays zero because the setup does not introduce a parameterized controller.
