```@meta
CurrentModule = Rible
ShareDefaultModule = true
```

```@setup
using Rible
import Rible as RB
using StaticArrays
using LinearAlgebra
using FileIO, MeshIO
using Makie
import CairoMakie
CairoMakie.activate!()
Makie.inline!(true)
```

# Bodies and Apparatuses

## Locus

!!! note "The place where something is situated or occurs. — Cambridge Dictionary"

In multibody dynamics, accurately describing the position and orientation of material points on a body is fundamental for establishing connections, applying external forces, and handling contacts. Additionally, physical properties such as surface normal direction, friction coefficient, and restitution coefficient at these points are essential for contact and collision processing.

To manage this, Rible introduces the **`Locus`** and **`LocusState`** types.

**`Locus`** is an immutable type that defines the invariant properties of a material point in the body's **local coordinate frame**:

- **`position`** — the point's position vector in local coordinates.
- **`axes`** — the local reference axis system (`Axes`) at this point, typically used to define contact surface normals or joint rotation/translation axes. By default, the first axis (x-axis) is the contact normal.
- **`friction_coefficient`** and **`restitution_coefficient`** — define the friction and restitution coefficients at this point.

### Constructing a Locus

When building a rigid body, you provide a `mass_locus` for the center of mass and a vector of `loci` for subsequent connections or contacts.

```@example
# Define loci in local coordinates
mass_locus = Locus(SVector(-0.2, 0.0, 0.0)) # center of mass, default friction/restitution

# Define two loci. The flexible constructor dispatch lets you provide parameters as needed
loci = [
    # Locus 1: position and reference axis normal
    Locus(SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0)),
    # Locus 2: position, reference axis normal, friction (0.5) and restitution (0.1)
    Locus(SVector(0.5, 0.0, 0.0), SVector(1.0, 0.0, 0.0), 0.5, 0.1)
]
```

**`LocusState`** is a mutable type that records the point's state in the **global (world) coordinate frame** during simulation:

- **`frame`** — the complete kinematic frame (`CartesianFrame`) in global coordinates, including current global position, orientation (rotation matrix), linear velocity, and angular velocity.
- **`force`** and **`torque`** — the global force and torque currently acting at this point.
- **`contact_state`** — maintains the current contact state, such as whether the point is in active contact, the contact gap, and contact forces.

## Bodies

!!! note "The body and its functions share a single source; the apparent and the hidden are not separate. — Cheng Yi, *Preface to the Commentary on the I Ching*"

Bodies are the fundamental building blocks of a multibody system. In Rible, all bodies are subtypes of `AbstractBody{N,T}`, where `N` is the spatial dimension (2 or 3) and `T` is the numeric type (typically `Float64`). This parametric design ensures the system handles both planar and spatial problems while maintaining type safety and performance.

### Rigid Body

A `RigidBody` is one of the most basic components in a multibody system. In Rible, it is composed of `RigidBodyProperty`, `RigidBodyState`, a coordinate description (Coordinates), and an optional visualization mesh.

#### Properties of Rigid Body

`RigidBodyProperty` defines the invariant properties of a rigid body. It contains dynamics parameters and simulation metadata:

- **Dynamics parameters**: mass (`mass`) and inertia tensor (`inertia`), where the inertia tensor is typically stored as an `SMatrix` for optimized memory and computation speed.
- **Geometric features**: center of mass position (`mass_locus`) and attachment points (`loci`).
- **Metadata**: unique `id` and optional `type` (symbol tag).
- **Flags**: `contactable` controls whether the body participates in contact detection; `visible` controls whether it is rendered in visualization.

```@example
id = 1
mass = 2.5
inertia = SMatrix{3,3,Float64}(I)  # unit inertia matrix

prop = RigidBodyProperty(
    id,
    true,           # contactable
    mass,
    inertia,
    mass_locus,
    loci;
    visible=true
)
```

#### State of Rigid Body

`RigidBodyState` maintains the body's instantaneous state. This includes the position vector and orientation (rotation matrix) in the global frame, plus linear and angular velocities. These are the quantities updated by the simulator at each time step.

```@example
r0 = SVector(0.0, 0.0, 0.0)      # initial position
R = SMatrix{3,3}(I)              # initial orientation (no rotation)
r_dot = zero(r0)                 # initial linear velocity
omega = zero(r0)                 # initial angular velocity
state = RigidBodyState(prop, r0, R, r_dot, omega)
```

#### Coordinate Description

The system uses the **Natural Coordinate Formulation (NCF)** to describe motion. Depending on the body type (e.g., a freely floating block, a link connecting two points), choose the appropriate coordinate formula:

- **`NC1P3V`** (1-point 3-vector): for standard 6-DOF rigid bodies
- **`NC3D2P`** (3D 2-point): for links defined by two endpoints
- **`NC2D2P`** (2D 2-point): for planar links

```@example
# For a general spatial rigid body, use one basic point and three base vectors (NC1P3V)
coords = NCF.NC1P3V(r0, r0, R)
```

#### Visualization Mesh

An optional mesh can be loaded for visualization.

```@example
# Load STL file and scale/translate (example)
# mesh = load(assetpath("link.STL")) |> make_patch(; scale=0.1)
mesh = nothing # mesh omitted here
```

#### Constructor

Finally, a rigid body can be constructed by supplying all the above instances.

```@example
body = RigidBody(prop, state, coords, mesh)
```

Under the hood, a `RigidBodyCache` is automatically constructed to handle intermediate computational cache.

## Apparatuses

!!! note "A great apparatus requires no completion. — *Daodejing* (Mawangdui Silk Manuscript)"

Apparatuses are composed of **Joints** and **Forces**, representing the physical joints or forces (such as springs) in the system. They are the carriers of **constraint terms** and **control forces** in the multibody dynamics equations.

### Forces

Forces compute the generalized forces acting on bodies. All force models subtype `AbstractForce`.

- **`DistanceSpringDamper`** — the most fundamental force element. It connects two points (Loci) and maintains a `DistanceSpringDamperState` that tracks the current length, elongation ratio, and tension.

- **`TorsionalSpringDamper`** and **`RotationalSpringDamper`** — produce restoring torques, typically acting on a joint's rotation axis to resist relative rotation.

### Joints

Joints define the relative connection between bodies. All joint models inherit from `AbstractJoint`.

Joints typically, but not always, introduce **constraint forces** into the multibody dynamics equations. Standard joints and their variants are implemented via the `ProtoJoint` type:

| Joint Type | DOF | Description |
| :--- | :---: | :--- |
| `:Fixed` | 0 | **Fixed joint**. Locks two rigid bodies completely, or locks a point on a rigid body to a fixed position in space (`FixedPointJoint`). Commonly used for base fixation. |
| `:Revolute` | 1 | **Revolute joint**. Allows rotation about a shared axis. Requires specifying coincident points on both bodies and collinear axes. |
| `:Prismatic` | 1 | **Prismatic joint**. Allows translation along an axis while restricting rotation. |
| `:Cylindrical` | 2 | **Cylindrical joint**. Allows both rotation about and sliding along an axis. |
| `:Universal` | 2 | **Universal joint**. Composed of two orthogonal revolute axes, allowing rotation in two directions. |
| `:Spherical` | 3 | **Spherical joint**. Constrains only translation, allowing arbitrary rotation about the connection point. |

#### Defining Connections

The core of joint construction is precisely describing the "connection relationship." This is achieved through two types:

**`Anchor`** — precisely locates the connection point's geometric information (position and constraint axes).

- `Anchor(body, pid)` — reads both position and axes from the same Locus.
- `Anchor(body, position_pid, axes_pid)` — position and axes come from different Loci (useful for revolute joints and other cases requiring a specific axis direction).

**`Hen2Egg`** — defines the pairing between the "parent body" (Hen) and "child body" (Egg).

- `Hen2Egg(hen_anchor, egg_anchor)` — establishes the correspondence between two anchors.

When building a joint, first use `Anchor` to specify the geometry at both ends, then use `Hen2Egg` to describe the topological direction.

#### Joint Construction Examples

```@example
prop1 = RigidBodyProperty(1, true, 1.0, inertia, mass_locus, loci)
state1 = RigidBodyState(prop1, r0, R, r_dot, omega)
nmcs1 = NCF.NC1P3V(r0, r0, R)
body1 = RigidBody(prop1, state1, nmcs1, nothing)

prop2 = RigidBodyProperty(2, true, 1.0, inertia, mass_locus, loci)
r0_2 = r0 + SVector(1.0, 0.0, 0.0)
state2 = RigidBodyState(prop2, r0_2, R, r_dot, omega)
nmcs2 = NCF.NC1P3V(r0_2, r0_2, R)
body2 = RigidBody(prop2, state2, nmcs2, nothing)
```

```@example
# 1. Define connection anchors: body1's 2nd locus and body2's 1st locus
hen = Anchor(body1, 2)
egg = Anchor(body2, 1)

# 2. Establish the connection
conn = Hen2Egg(hen, egg)
```

If the joint requires a specific rotation or translation axis (e.g., revolute joint), use `Anchor(body, position_pid, axes_pid)` to specify position and constraint axis separately:

```@example
# Specify body1's 1st locus for position/axis, body2's 2nd locus for position/axis
hen_axis = Anchor(body1, 1, 1)
egg_axis = Anchor(body2, 2, 2)
conn_axis = Hen2Egg(hen_axis, egg_axis)
```

Create the joint object using `ProtoJoint` with a joint type symbol:

```@example
# Create a spherical joint
joint_spherical = ProtoJoint(conn, :Spherical)

# Create a revolute joint
joint_revolute = ProtoJoint(conn_axis, :Revolute)
```
