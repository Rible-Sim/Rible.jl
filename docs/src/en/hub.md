```@meta
CurrentModule = Rible
ShareDefaultModule = true
```

# Control Hub

If `Structure` describes the robot's "body" — the musculoskeletal system — then `ControlHub` defines its "nervous system": it specifies what the robot can sense and what it can actuate.

1. **Sensing via Gauges**: Gauges define how to extract feature vectors from the complex physical state (coordinates, velocities) for feedback control policies or value functions.
2. **Actuate via Actuators**: define how to transform abstract control commands into generalized forces in the physical system.

This architecture allows the same physical model to seamlessly work with different control interfaces, facilitating control algorithm implementations.

## 2. Actuators

All actuators inherit from the abstract type `AbstractActuator`. An actuator's responsibility is to map the actions to forces applied in the physical system.

- **`ExternalForceActuator`**: the most commonly used concrete implementation, for applying forces or torques at a specific Locus on a body.

Actuator behavior is defined through the following functions:

- **`actuate!(bot, policy, state)`**: Top-level entry point. Computes the actions for the entire robot based on the current policy and triggers all sub-actuators.
- **`execute!(structure, actuator, u)`**: Executes a single actuator, converting action `u` to generalized forces and accumulating them into the system force vector.
- **`gen_force_actu_jacobian!(∂F∂u, structure, actuator, u)`**: Computes the Jacobian of the generated generalized forces with respect to actions `u`. This is the foundation for adjoint computations.

## 3. Gauges

Gauges define how the robot observes itself and its environment. All gauges inherit from `AbstractGauge`.

### Captum Gauges vs. Error Gauges

To serve both feedback control and optimization solving, Rible classifies gauges into two categories:

#### CaptumGauge

`CaptumGauge` retrieves **raw feedback signals**. It is typically fed to a policy function or used for state monitoring.

- **Composition**: a `Signifier` (pointing to a body locus) and a `Captum` (e.g., `PositionCaptum`).
- **Output**: typically a vector (e.g., a 3D position vector).

#### Error Gauge

`ErrorGauge` introduces a **reference value** on top of the observation. It directly computes the deviation between the current state and the reference.

- **Mathematical form**: typically computes a weighted squared error

```math
e = \frac{1}{2}\|a - a_{ref}\|^2
```

- **Output**: a scalar.
- **Usage**: `ErrorGauge` is the atomic building block for constructing Objectives (cost functions for trajectory optimization).

### Core Interface: Measurement

Gauge computation revolves around the following interface:

- **`measure!(out, st, sig, cap)`**: Writes physical attributes (e.g., position or velocity at a specific locus) directly into the output buffer `out`, based on the concrete `Signifier` and `Captum` types.
- **`measure_jacobian!(Jq, Jv, Js, st, sig, cap)`**: Computes the Jacobian of the measurement with respect to system generalized coordinates ``q``, velocities ``v``, and auxiliary variables ``s``. 

```@setup
using Rible
using StaticArrays
using LinearAlgebra
import TypeSortedCollections as TSC
```

## 4. Example: Building a Simple Cart-Pole

The following demonstrates how to build a `ControlHub` for a simplified cart-pole system.

### Step 1: Physical Structure Setup

First, create a physical structure (`Structure`) for the robot. In this example, we define a cart with a pole.

```@example cartpole
using Rible
using StaticArrays
using LinearAlgebra
import TypeSortedCollections as TSC

# Physical structure setup (cart, L=1.0)
l = 1.0
mass_locus = Locus(@SVector zeros(3))
loci = [Locus(@SVector zeros(3)), Locus(SVector(0.0, l, 0.0))]
prop_cart = RigidBodyProperty(1, true, 1.0, SMatrix{3,3,Float64}(I), mass_locus, loci)
prop_pole = RigidBodyProperty(2, true, 0.5, SMatrix{3,3,Float64}(I), mass_locus, loci)

cart = RigidBody(prop_cart, RigidBodyState(prop_cart, zeros(3), SMatrix{3,3,Float64}(I), zeros(3), zeros(3)), NC1P3V(zeros(3), zeros(3), SMatrix{3,3,Float64}(I)), nothing)
pole = RigidBody(prop_pole, RigidBodyState(prop_pole, SVector(0.0, l/2, 0.0), SMatrix{3,3,Float64}(I), zeros(3), zeros(3)), NC1P3V(SVector(0.0, l/2, 0.0), SVector(0.0, l/2, 0.0), SMatrix{3,3,Float64}(I)), nothing)

bodies = TSC.TypeSortedCollection([cart, pole])
apparatuses = TSC.TypeSortedCollection(Int[])
st_obj = Structure(bodies, apparatuses, Connectivity(bodies, apparatuses))
```

### Step 2: Define Gauges

Next, define two types of gauges: one for observing the cart's state, and one for measuring the pole's position error.

```@example cartpole
# Define gauges: observation (Captum) and error (Error)
cart_captum = CaptumGauge(1, Signifier(cart, 1), PosVelCaptum())
pole_pos_error = ErrorGauge(1, Signifier(pole, 2), PositionCaptum(), [0.0, l, 0.0])

capta_gauges = TSC.TypeSortedCollection([cart_captum])
error_gauges = TSC.TypeSortedCollection([pole_pos_error])
```

### Step 3: Define Actuators

Actuators define how action values affect the physical system. Here we define an actuator that applies a horizontal push force on the cart.

```@example cartpole
# Define actuator: horizontal push force
force_actuator = ExternalForceActuator(
    1, Signifier(cart, 1), NaiveOperator(1), [0; 1.0; 0;;], [0.0]
)
actuators = TSC.TypeSortedCollection([force_actuator])
```

### Step 4: Assemble the Hub and Measure

Finally, assemble all components into a `ControlHub` and demonstrate measurement.

```@example cartpole
# Assemble Hub and measure
coalition = Coalition(st_obj, capta_gauges, error_gauges, actuators)
hub = ControlHub(st_obj, capta_gauges, error_gauges, actuators, coalition)

# Execute measurements
c = measure(st_obj, cart_captum)
e = measure(st_obj, pole_pos_error)
(c, e)
```
