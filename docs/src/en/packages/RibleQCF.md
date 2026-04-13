# Quaternion Coordinates — RibleQCF

```@meta
CurrentModule = Rible
ShareDefaultModule = true
```

```@setup qcf
using Rible
import Rible as RB
using RibleQCF
import RibleQCF as QCF
using StaticArrays
using LinearAlgebra
using Rotations
using TypeSortedCollections
using Makie
import CairoMakie
using FileIO
using MeshIO
import GeometryBasics as GB

CairoMakie.activate!()
Makie.inline!(true)
```

RibleQCF is an extension package that adds quaternion-based coordinates (`QC`) as an alternative to the default Natural Coordinate Formulation (NCF). A single rigid body is described by **7 generalized coordinates** — 3 for translation and 4 for a unit quaternion — yielding 6 degrees of freedom after the unit-norm constraint is enforced. The inertia representation and its numerical influence follow the framework established in [xuNumericalInfluencesInertia2016](@cite) and [xuNumericalInfluenceAdditional2020](@cite).

!!! tip 
    **QCF** produces a non-constant mass matrix, which requires recomputations during updates and has involved jacobian formulations.

## Coordinate Type `QC`

The core type `QCF.QC` stores mass and a diagonalized inertia tensor along with several precomputed scalars derived from them:

```julia
struct QC{T,JT} <: AbstractCoordinates{3,T}
    m      # mass
    m⁻¹    # 1/mass
    γ      # max diagonal of inertia (regularization parameter)
    Jγ     # J - γ·I
    γ⁻¹    # 1/γ
    J⁻¹γ   # J⁻¹ - (1/γ)·I
end
```

Construct a `QC` instance by passing mass and a 3×3 inertia matrix:

```julia
m = 0.5
J = Diagonal(SA[0.001, 0.001, 0.002])

qc = QCF.QC(m, J)
```

Only the diagonal of the inertia matrix is used; off-diagonal terms are discarded internally. The optional keyword `γ` defaults to `maximum(diag(J))` and acts as a regularization parameter that keeps the rotational mass matrix well-conditioned.

!!! tip "Checking the coordinate trait"
    You can query whether a body uses constant-mass coordinates at any level:
    ```julia
    has_constant_mass_matrix(body.coords)   # Val(false) for QC
    has_constant_mass_matrix(body)          # propagates from coords
    has_constant_mass_matrix(structure)     # propagates from bodies
    ```

## Nonconstant Mass Matrix

Unlike NCF (where `has_constant_mass_matrix` returns `Val(true)`), the QCF generalized mass matrix **depends on the current quaternion state**:

```julia
has_constant_mass_matrix(::QCF.QC) = Val(false)
```

The 7×7 mass matrix is block-diagonal:

```math
M(q) = \begin{bmatrix} mI_3 & 0 \\ 0 & M_\gamma(q) \end{bmatrix}, \quad M_\gamma(q) = 4\bigl(L^\top(q)\,J_\gamma\,L(q) + \gamma I_4\bigr)
```

where ``L(q)`` is the 4×4 left-quaternion-multiplication matrix. Because ``M`` changes every step, the solver must recompute and refactorize it at each iteration. Rible handles this automatically through its solver dispatch: when it detects `Val(false)`, it selects the `Zhong06_Nonconstant_Mass` solver family instead of the constant-mass family.

### Practical implications

- **Cost per step** is higher than NCF due to per-step mass-matrix assembly and factorization. For single-body simulations the difference is negligible; for large multi-body systems with many QCF bodies, consider whether NCF suffices.
- **Tight tolerances help**. The quaternion constraint ``\|q\|^2 = 1`` drifts if the solver tolerances are loose. Using `ftol = 1e-14` and `maxiters = 50` is a good default.

## Intrinsic Constraint

`QC` carries one intrinsic algebraic constraint — the unit quaternion norm:

```math
c(q) = \frac{1}{2}(q^\top q - 1) = 0
```

This constraint is enforced automatically by the solver alongside any joint constraints. The framework provides the constraint function, its Jacobian, and the force-Jacobian in `constraints.jl`.

## Complete Example: Spinning Top with QCF

The bundled spinning-top model supports both NCF and QCF via the `cT` argument. Below we construct a QCF top, verify its nonconstant-mass trait, and run a full contact simulation.

### Construction

```@example qcf
include(joinpath(pathof(Rible), "../../examples/robots/spinningtop.jl"))

origin_position = [0.0, 0.0, 0.5]
R = RotX(0.0)
origin_velocity = [1.0, 0.0, 0.0]
Ω = [0.0, 0.0, 200.0]

top = make_top(origin_position, R, origin_velocity, Ω, :QCF; μ = 0.95, e = 0.5, loadmesh = true)
```

### Verify the nonconstant-mass trait

```@example qcf
body1 = get_bodies(top.structure)[1]
has_constant_mass_matrix(body1.coords) === Val(false)
```

### Contact setup and simulation

The solver object is the same as for NCF — a `Zhong06` integrator with an inner contact solver. Internally, the framework routes QCF bodies to the nonconstant-mass code path.

```@example qcf
planes = StaticContactSurfaces([
    HalfSpace([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
])

contact_model = RestitutionFrictionCombined(
    NewtonRestitution(),
    CoulombFriction(),
)

prob = DynamicsProblem(top; env=planes, contact_model)

solver = DynamicsSolver(
    Zhong06(),
    InnerLayerContactSolver(InteriorPointMethod()),
)

sim_result = solve!(
    prob,
    solver;
    tspan = (0.0, 0.2),
    dt = 1e-3,
    ftol = 1e-14,
    maxiters = 50,
    exception = false,
)

tip_pos = get_trajectory!(top, 1, 1).u[end]
tip_vel = get_velocity!(top, 1, 1).u[end]
tip_pos, tip_vel
```

### Visualization

Because the top was constructed with `loadmesh = true`, `vis!` renders the STL geometry directly:

```@example qcf
key_steps = round.(Int, range(1, length(top.traj.t), length = 4))

fig_traj = plot_traj!(
    top;
    do_slide = false,
    show_info = false,
    show_loci = false,
    show_background = false,
    show_ground = false,
    gridsize = (2, 2),
    at_steps = key_steps,
)
```

### Comparing QCF and NCF results

You can reuse the same initial conditions with `:NCF` to compare:

```@example qcf
top_ncf = make_top(origin_position, R, origin_velocity, Ω, :NCF; μ = 0.95, e = 0.5, loadmesh = false)

has_constant_mass_matrix(get_bodies(top_ncf.structure)[1].coords) === Val(true)
```

## References

```@bibliography
Pages = ["RibleQCF.md"]
```
