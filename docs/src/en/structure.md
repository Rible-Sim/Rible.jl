```@meta
CurrentModule = Rible
ShareDefaultModule = true
```

```@setup
using Rible
import Rible as RB
using StaticArrays
using LinearAlgebra
```

# Structure

`AbstractStructure`/`Structure` is the central container that manages all mechanical properties, topology, and kinematic state of a multibody system. It integrates bodies and apparatuses, uses a semi-auto-generated `Connectivity` for numbering and indexing, and maintains `StructureState` and `StructureCache` that exposes an internal interface for dynamics.

Rible's interface design relies on generic functions and multiple dispatch. The family of functions defined around `AbstractStructure` (e.g., `update!`, `cstr_function!`) establishes behavioral contracts — developers implement the corresponding methods for new subtypes in order to plug into the simulation framework seamlessly.

- `Structure` is the standard concrete implementation, suitable for most scenarios involving rigid bodies and complex connections.
- Users can extend functionality by subtyping `AbstractStructure` and implementing specific methods. For example, `TensegrityStructure` in `RibleTensegrity` provides specialized dispatches for tensegrity structures.

## Composition

A structure consists of three typical parts:

### Bodies and Apparatuses

All motion units in the system are stored in the `bodies` field. Because different bodies may use different mathematical descriptions, this is typically a heterogeneous collection. All connection elements (springs, dampers, constraint joints, etc.) are stored in the `apparatuses` field.

To maintain type stability and high performance when processing heterogeneous collections, Rible uses `TypeSortedCollection` for storage. This ensures that bodies and apparatuses of the same type trigger specialized method dispatch, avoiding the overhead of dynamic dispatch during large-scale system updates.

### Structure State

`StructureState` acts as the canonical structural state, managing variables through multiple views:

- **System view**: stores global coordinates `q`, velocities ``q̇``, accelerations ``q̈``, constraint forces `λ`, etc. in flat arrays under the `system` field.
- **Member view**: maps slices of the global arrays to individual bodies via `view` mechanics (stored in `members`), so that components operate directly on shared global memory without data duplication.

### Coordinates State

`AbstractCoordinatesState` is the abstract supertype for all coordinate-level state containers used by dynamics solvers. Its standard concrete implementation is `CoordinatesState`, which stores:

- **`q`**, **`q̇`**, **`q̈`** — generalized coordinates, velocities, and accelerations as flat vectors.
- **`p`** — conjugate momentum.
- **`F`** — the generalized force vector assembled from body forces, apparatus forces, and external fields.
- **`λ`** — Lagrange multipliers (constraint forces).
- **`s`** — auxiliary variables for apparatus internal states.
- **`c`** — structural parameters.

`StructureState.system` is a `CoordinatesState`, and solvers operate on `CoordinatesState` directly when calling structure-level functions such as `cstr_function!` and `assemble_forces!`. This separation lets the solver manipulate a lightweight state object while the structure maintains the full member-level mapping internally.

## Connectivity

`Connectivity` is typically auto-generated from the topological layout of bodies and apparatuses. It manages the complex index mapping logic.

### Auto-generation

Calling `Connectivity(bodies, apparatuses)` traverses all components and automatically derives the system's global degrees of freedom, index offsets, and connection topology.

### Internal Index Mapping

`Connectivity` maintains key index information:

- **Coordinate indices**: maps each body's coordinates to their position in the system-level coordinate vector.
- **Constraint indices**: distinguishes between a body's intrinsic constraints (e.g., body constraints from natural coordinates) and extrinsic constraints introduced by apparatuses (e.g., kinematic pairs).
- **Position mapping**: establishes the association between body attachment points (Loci) and the global system state.

## Core Operations

The structure provides core functions for dynamics simulation.

### Update Pipeline

`update!` is the top-level entry point for a full structure update. It executes the following steps in order:

```julia
function update!(st::AbstractStructure, field=NoField(); )
    clear_forces!(st)
    stretch!(st)
    update_bodies!(st)
    update_apparatuses!(st)
    apply_field!(st, field)
    assemble_forces!(st)
end
```

Each step has a specific role:

- **`clear_forces!(st)`** — zeroes out all force vectors in both the system-level state and individual body caches.
- **`stretch!(st)`** — propagates structual parameter values to each body's loci and each apparatus, preparing them for evaluation.
- **`update_bodies!(st)`** — populates each body's local state from system-level coordinates, recomputes inertia caches, transformations, and loci states.
- **`update_apparatuses!(st)`** — evaluates all apparatuses (joints, forces) using the updated body states.
- **`apply_field!(st, field)`** — applies external fields (e.g., gravity) to each body.
- **`assemble_forces!(st)`** — collects all body-level and apparatus-level forces into the system-level force vector `F`.

A lighter variant, `lazy_update!`, replaces `update_bodies!` with `lazy_update_bodies!`, which skips the inertia cache recomputation. This is useful when the mass matrix is not needed (e.g., during iterative constraint solving).

### Constraints and Jacobians

- **`cstr_function!(Φ, st, inst_state)`**: Computes the residual vector of all active constraints in the system.
- **`cstr_jacobian!(A, st, inst_state)`**: Computes the Jacobian matrix of system constraints with respect to generalized coordinates.

### Mass Matrix and Energy

- **`assemble_M(st)`**: Assembles and returns the system-level mass matrix.
- **`kinetic_energy(st)`**: Computes the system's total kinetic energy from the current state via dispatch.

Through these highly abstract generic interfaces, dynamics solvers can process complex physical systems with a unified logic.

## Example: Minimal Runnable Structure

The following example uses only the package API to build a minimal system.
This page enables `ShareDefaultModule = true`, so multiple `@example` blocks share variables sequentially and inline results at each step.

### Step 1: Define a Rigid Body

```@example
using Rible
using StaticArrays
using LinearAlgebra

mass_locus = Locus(SVector(0.0, 0.0, 0.0))
loci = [Locus(SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0))]
inertia = SMatrix{3,3,Float64}(I)

prop = RigidBodyProperty(1, true, 1.0, inertia, mass_locus, loci; visible = false)
r0 = SVector(0.0, 0.0, 0.0)
R0 = SMatrix{3,3,Float64}(I)
state = RigidBodyState(prop, r0, R0, zero(r0), zero(r0))
coords = NCF.NC1P3V(r0, r0, R0)
body = RigidBody(prop, state, coords, nothing)
```

### Step 2: Assemble Connectivity and Structure

```@example
bodies = [body]
apparatuses = Int[]
cnt = Connectivity(bodies, apparatuses)
st = Structure(bodies, apparatuses, cnt)
inst = st.state.system
(; n_bodies = length(st.bodies), n_apparatuses = length(st.apparatuses))
```

### Step 3: Access System State Vectors

```@example
q = inst.q
qdot = inst.q̇
(; q_len = length(q), qdot_len = length(qdot))
```

### Step 4: Run a Structure Update and Query Dimensions

```@example
update!(st)
n_q = get_num_of_coords(st)
n_c = get_num_of_cstr(st)
(; n_q, n_c, dof = n_q - n_c)
```

### Step 5: Compute Constraint Residual and Jacobian

```@example
Φ = zeros(n_c)
cstr_function!(Φ, st, inst)
A = cstr_jacobian(st, inst)
(; residual_norm = norm(Φ), jacobian_size = size(A))
```

### Step 6: Assemble the Mass Matrix

```@example
M = assemble_M(st)
size(M)
```

### Step 7: Build a Robot and Run a Short Integration

```@example
bot = Robot(st)
prob = DynamicsProblem(bot; env=GravityEnv())
solver = DynamicsSolver(Zhong06())
solve!(prob, solver; dt = 1e-3, tspan = (0.0, 0.01))
length(bot.traj)
```

### Step 8: Post-processing — Mechanical Energy and Mid-point Velocity

```@example
E_total = mechanical_energy!(bot).E[end]
v_mid = first(get_mid_velocity!(bot, 1, 1))
(; E_total, v_mid)
```
