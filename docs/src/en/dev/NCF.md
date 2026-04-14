# Natural Coordinates

Natural Coordinates is the builtin coordinate formulations used in Rible for modeling rigid bodies. It is implemented as a submodule `NCF` living in `src/coordinates/natural_coordinates`, and is wired through multiple dispatch to the mass matrix, constraint equations, and joint solving workflows.

> Core features: This approach describes pose directly using point coordinates and direction vectors. It features a non-singular orientation description and a constant mass matrix.

## Basic Representation

In the natural coordinates framework, the state of a rigid body is described by a basic point and direction vectors:

- 2D bar representation: 1 base point + 1 base vector
- 2D rigid body representation: 1 base point + 2 base vectors
- 3D rigid body representation: 1 base point + 3 base vectors

Variants, such as 2 base points + 2 base vectors, are also fully supported to meet special needs.

Correspondingly, `NC` is a unified coordinate structure; the specific dimension and degrees of freedom are determined by the constructor.

## Module Structure and Code Entry Points

The main files under `src/coordinates/natural_coordinates` are:

- `types.jl`: defines `NC` and underlying data structures.
- `constructors.jl`: defines common constructors, such as `NC2D1P1V`, `NC1P2V`, `NC1P3V`.
- `functions.jl`: coordinate mappings and Jacobian calculations, e.g., `to_position`, `to_position_jacobian`.
- `constraints.jl`: intrinsic constraints and Jacobian implementations, core interfaces `cstr_function!`, `cstr_jacobian!`.
- `mass_matrix.jl`: mass matrix construction, e.g., `make_M`.
- `joints.jl`: joint caches, constraint violations, and Jacobian assembly.

## Common Usage

### 1 Construct Coordinate Objects

Constructors of natural coordinates expect at least a basic point `ri` and optionally other basic points and base vectors like `u`, all expressed in the global frame.

```julia
# 2D bar: 1 point + 1 vector
nc_bar_2d = NC2D1P1V(ri, u)

# 2D rigid body: 1 point + 2 vectors
nc_body_2d = NC1P2V(ri)

# 3D rigid body: 1 point + 3 vectors
nc_body_3d = NC1P3V(ri)
```

### 2 Position Mapping and Jacobian

```julia
# Position of local point c in global frame
pos = to_position(nc_body_3d, q, c)

# Jacobian with respect to coordinates
Jpos = to_position_jacobian(nc_body_3d, q, c)
```

### 3 Constraint Values and Constraint Jacobians

`constraints.jl` core interface uses in-place update variants.

```julia
ret_c = zeros(nc_num_constraints)
ret_J = zeros(nc_num_constraints, length(q))

cstr_function!(ret_c, nc_body_3d, q)
cstr_jacobian!(ret_J, nc_body_3d, similar(ret_J), q)
```

### 4 Mass Matrix

Natural coordinates produce a constant mass matrix for a rigid body.

```julia
M = make_M(nc_body_3d, mass, inertia, mass_center)
```

## Joint Development Guide

This section targets developers who need to add new joint types. Rible's joint implementation uses a "base constraint composition" strategy: first define the basic constraint types, then in `src/joint/joints.jl`’s `get_joint_info` combine masks to construct the specific joint.

### Theoretical Foundation: Basic Constraints

In the joint module, all constraints are written in a unified "joint coordinates" form. Let

- `q_h`: parent (hen) body coordinates
- `q_e`: child (egg) body coordinates
- `q_j = [q_h; q_e]`: concatenated joint coordinates

`build_joint_cache` aims to precompile geometric constraints into two matrix objects:

- `transformations`: first-order (linear) term matrices
- `halves`: list of quadratic form matrices (one per constraint is `H_i`)

Thus each constraint can be written as:

```math
C_i(q_j) = q_j^T H_i q_j + T_i q_j - c_i
```

Where `T_i` is the i-th row of `transformations`, and `H_i` is `halves[i]`. `violations` stores constant offsets `c_i` for the reference configuration, which are subtracted at runtime to center the baseline at zero violation.

> Linear vs quadratic meaning:
> - Linear constraint: `C_i(q_j)` is linear in `q_j`; the Jacobian is constant.
> - Quadratic constraint: `C_i(q_j)` is a quadratic polynomial in `q_j`; the Jacobian is linear in `q_j`, but the Hessian is constant.

### Mask 1st Coincidence Constraint Type 1

Constraint that two points coincide:

```math
r_i - r_j = 0
```

In code this constraint mainly comes from the linear part of the point-mapping Jacobian concatenation `J = [-C_h, C_e]`:

- Code index: `mask_1st`
- Matrix contribution: added to `transformations = J[mask_1st, :]`
- Constraint order: Linear (corresponding to `H_i = 0`, only `T_i q_j` remains)

### Mask 2nd Distance Constraint Type 4

Constraint that the distance between two points is a constant `d`:

```math
\|r_i - r_j\|^2 - d^2 = 0
```

These constraints are constructed in code via a quadratic form:

- Code index: `mask_2nd`
- Matrix contribution: `half_2nd[1] = J' * J`
- Constraint order: Quadratic (only `q_j^T H_i q_j`; the constant term is provided by `violations`)

### Mask 3rd Orthogonal Constraint I Type 3

Constraint that the direction vector is orthogonal to the relative displacement:

```math
 v_h^T (r_i - r_j) = 0
```

Its quadratic form is obtained by combining the direction-selection matrix with a Kronecker product:

- Code index: `mask_3rd`
- Matrix construction: generated via `kron(...)` to create `half_3rd`
- Constraint order: Quadratic
- Common uses: planar motion constraints, slider direction constraints

### Mask 4th Orthogonal Constraint II Type 2

Constraint that two direction vectors are orthogonal:

```math
 v_h^T v_l = 0
```

This constraint is entirely formed from the orientation coordinates, i.e., a pure pose-coupling term:

- Code index: `mask_4th`
- Matrix construction: built by combining orientation bases with `kron(...)` to form `half_4th`
- Constraint order: Quadratic
- Common uses: shaft alignment, universal joint constraint composition

## Implementation Details and Performance

### Cache Mechanism: ApparatusCache

`build_joint_cache` returns `cache` and `violations` at initialization. The `cache` contains all runtime sparse structures:

- `transformations`: linear constraint matrices
- `halves`: list of quadratic form matrices
- `hessians`: each constraint’s symmetric Hessian `H_i + H_i^T`
- `joint_q`, `joint_work`, `trf_work`: in-place computation buffers to avoid allocations

This design ensures the runtime only performs value substitutions; the structure matrices are not rebuilt.

### Key Functions in joints.jl

- `build_joint_cache`: precompute `transformations`, `halves`, `hessians` according to the masks.
- `get_joint_violations!`: first compute `q_j^T H_i q_j`, then add the linear term `T_i q_j`, and finally subtract `violations`.
- `get_joint_jacobian!`: assemble `C_q` using the precomputed `hessians` with the current `q_j`.

## Extending to New Joints

Extending to new joints typically does not require new underlying mathematical kernels. A recommended workflow:

1. In `src/joint/joints.jl`, in `get_joint_info` define the target joint’s mask combination.
2. Use the existing `mask_1st` through `mask_4th` to compose the desired degrees-of-freedom constraint.
3. Reuse `build_joint_cache`, `get_joint_violations!`, and `get_joint_jacobian!` to hook into the solver.

This approach keeps new joints fully compatible with the existing solving pipeline and leverages the caching performance benefits.
