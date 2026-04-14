# RibleTensegrity

```@meta
ShareDefaultModule = true
```

```@setup tensegrity
using Rible
import Rible as RB
using RibleTensegrity
import RibleTensegrity as RT
using StaticArrays
using LinearAlgebra
using SparseArrays
using TypeSortedCollections
using CircularArrays
import CircularArrays as CA
using Rotations
using Unitful
using ElasticArrays
using DataStructures
using EponymTuples
using StructArrays
using FileIO, MeshIO
import GeometryBasics as GB
using LoggingExtras
using Polyhedra
import CDDLib
using Meshes
using Makie
import CairoMakie

CairoMakie.activate!()
Makie.inline!(true)
global_logger(ConsoleLogger(stdout, Logging.Warn; show_limited=false))

# Body definitions (from main Rible repo)
include(joinpath(pathof(Rible), "../../examples/bodies/rigidbar.jl"))
include(joinpath(pathof(Rible), "../../examples/bodies/rigidbar_nonsmooth_repro.jl"))
include(joinpath(pathof(Rible), "../../examples/bodies/make_3d_bar.jl"))
include(joinpath(pathof(Rible), "../../examples/bodies/make_3d_plate.jl"))
include(joinpath(pathof(Rible), "../../examples/bodies/make_3d_tri.jl"))
include(joinpath(pathof(Rible), "../../examples/bodies/new_deck.jl"))

# Robot definitions
include(joinpath(pathof(Rible), "../../examples/robots/superball.jl"))
include(joinpath(pathof(RibleTensegrity), "../../examples/robots/bridge3d.jl"))

const nofield = RB.NoField()
const gravity = RB.Gravity()
nothing
```

RibleTensegrity is an extension package specialized for modeling and analyzing tensegrity structures [luoUnifiedApproachDynamic2024a](@cite), which originate from architecture innovations by Kenneth Snelson and Buckminster Fuller. Their core feature is the idea of tensional integrity: discrete compression members (struts) suspended within a continuous network of cables under tension.

In robotics, tensegrity structures are prized for their high strength-to-weight ratio, excellent energy absorption, large deformation capabilities, and inherent collision resilience. They are considered an ideal framework for rigid-flexible robots. This package provides an toolchain from morphology design and topology expression, through inverse static form finding, to dynamic simulation and stiffness optimization.

## Core Types

### TensegrityStructure

`TensegrityStructure` is the central container that manages all rigid bodies and cable apparatuses. It is a subtype of `AbstractStructure` and holds the system state, connectivity, and contact-related data:

```julia
struct TensegrityStructure{BodyType,TenType,CntType,StateType,CRType,CacheType} <:
    AbstractStructure{BodyType,TenType,CntType}
```

Constructing a `TensegrityStructure` requires three arguments — bodies, apparatuses, and a connectivity object:

```julia
st = RT.TensegrityStructure(bodies, apparatuses, connectivity)
```

Unlike vanilla structures, `TensegrityStructure` treats force assembly with special consideration for the **unilateral loading** characteristics of cables: cables only pull, never push.

### CableJoint

`CableJoint` is the key type implementing cable logic. Unlike a standard hinge, it does **not** introduce kinematic hard constraints — it only defines the direction vector of the applied force:

```julia
struct CableJoint{hen2eggType} <: AbstractJoint
```

Each `CableJoint` connects a pair of loci on two bodies via a `Hen2Egg` pairing. The unilateral elastic rule that governs cable force is:

```math
f_{\text{tension}} = \max(0, k(l - l_0) + c \dot{l})
```

where ``k`` is stiffness, ``l`` the current length, ``l_0`` the rest length, and ``c`` the damping coefficient. This models the slack property: when the cable is shorter than its rest length, it exerts zero force.

## Structure Construction

The package provides helper functions to streamline construction. The typical workflow is:

1. **Define rigid bodies** (e.g., bars, plates)
2. **Specify a connecting matrix** indicating which body loci are linked by cables
3. **Convert format** and **create cables** using `change_connecting_format` + `connect_spring`
4. **Contructing** a `TensegrityStructure`

### Example: Superball (6-bar Tensegrity)

The Superball is a classic 6-bar tensegrity robot with 24 cables. The bundled `superball()` constructor demonstrates the full construction pipeline:

```@example tensegrity
# Construct a superball with default parameters
bot = superball(0.0;
    l = 1.7 / 2,        # half bar length
    d = 1.7 / 4,        # half bar distance
    z0 = 2.0,
    visible = false,
    loadmesh = false,
    constrained = false,
)

st = bot.structure

println("Number of bodies: $(length(st.bodies))")
println("Number of cables: $(length(RT.get_cables(st)))")
println("Full coordinates: $(st.connectivity.num_of_full_coords)")
```

Internally, `superball()` follows this pattern:

```julia
# 1. Define rigid bodies using the rigidbar helper
rbs = [rigidbar(i, p₁, p₂; m=5.0, μ, e) for i in 1:6]
rigidbodies = TypeSortedCollection(rbs)

# 2. Define spring-damper devices for each cable
spring_dampers = [RB.DistanceSpringDamper3D(restlen, k, c; slack=false) for i in 1:24]

# 3. Build the connecting matrix and convert format
#    Format: each row is [bid₁, pid₁, bid₂, pid₂]
#    positive entry = Hen endpoint, negative = Egg endpoint
cm = RT.change_connecting_format(rigidbodies, connecting_matrix)

# 4. Create CableJoint + DistanceSpringDamper pairs
cables = RT.connect_spring(rigidbodies, spring_dampers; cm, istart=0)

# 5. Choose connectivity type and assemble
cnt = RB.Connectivity(rigidbodies, cables)           # free structure
# or: cnt = RB.PresFreeConnectivity(rigidbodies, cables)  # prescribed-free
st = RT.TensegrityStructure(rigidbodies, cables, cnt)
```

### Cable Inspection API

After construction, use these query functions to inspect cable properties:

```@example tensegrity
# Cable tensions, lengths, stiffness at the current configuration
tensions  = RT.get_cables_tension(st)
lengths   = RT.get_cables_len(st)
stiffness = RT.get_cables_stiffness(st)
restlens  = RT.get_cables_restlen(st)

println("Tension range: [$(minimum(tensions)), $(maximum(tensions))]")
println("Length range:  [$(minimum(lengths)), $(maximum(lengths))]")
```

Other inspection functions include `get_cables_deform`, `get_cables_len_dot`, `get_cables_force_density`.

## Inverse Statics and Form Finding

Given a target configuration ``q``, the static equilibrium condition balances cable forces against external loads:

```math
B \gamma = \tilde{F}_{ext}
```

where ``B \in \mathbb{R}^{n \times m}`` is the equilibrium matrix mapping ``m`` cable force densities ``\gamma`` to ``n`` system DOFs, and ``\tilde{F}_{ext}`` is the generalized external force (e.g., gravity).

### Inverse for Rest Length

Given a target configuration, `inverse_for_restlength` solves for the cable rest lengths that achieve static equilibrium under a specified external field:

```julia
μ = RT.inverse_for_restlength(bot, botref, field; fmin=0.0, eps_rel=1e-6, verbose=false)
```

- `bot`: the robot/structure at the target configuration
- `botref`: reference structure (usually same as `bot`)
- `field`: external field (e.g., `RB.NoField()`, `RB.Gravity()`)
- `fmin`: minimum pretension force to ensure no cable goes slack

For determinate systems, it solves directly. For indeterminate or unilateral-constrained systems (``γ ≥ 0``), it automatically uses Quadratic Programming via COSMO:

```math
\min_{\mu} \frac{1}{2} \mu^T \text{diag}(\kappa)\,\mu + h^T\mu \quad \text{s.t.} \quad B\mu = \tilde{F}_{\text{ext}},\;\mu \geq 0,\;\kappa \cdot (\ell - \mu) \geq f_{\min}
```

- **`inverse_for_actuation(bot, target_q, field)`**: for structures with actuators (winches, motors), solves for the actuation variables needed to reach the target state.

### Example: Tensegrity Bridge

```@example tensegrity
# Build a 2-module tensegrity bridge
bridge = bridge3d(; n=2)

# Solve for rest lengths with minimum pretension
μ = RT.inverse_for_restlength(bridge, bridge, nofield; fmin=1e4, eps_rel=1e-10)

# Apply the solved rest lengths
RT.set_restlen!(bridge, μ)
RB.update!(bridge.structure, nofield)

# Verify tensions are within expected range
tension_min, tension_max = RT.get_cables_tension(bridge.structure) |> extrema
println("Tension range: [$tension_min, $tension_max]")
```

### GDR: Generalized Dynamic Relaxation

For complex structures where analytical inversion is difficult, `GDR!` (provided by the core `Rible` package) performs numerical relaxation to find equilibrium. It introduces fictitious kinetic energy dissipation to drive the structure toward minimum-energy equilibrium:

```@example tensegrity
# Run GDR under gravity to find the gravity-balanced equilibrium
bridge_gdr = deepcopy(bridge)
RB.GDR!(bridge_gdr, gravity; β=1e-4, res=1e-9)

# Extract the equilibrium configuration
q_eq = bridge_gdr.traj.q[end]
RB.set_initial!(bridge_gdr, q_eq, zero(q_eq))

# Verify static equilibrium
isequilibrium, _ = RT.check_static_equilibrium_output_multipliers(
    bridge_gdr.structure, gravity
)
println("At equilibrium: $isequilibrium")
```

### Eigenvalue Analysis

Once at equilibrium, you can compute the structure's natural frequencies:

```@example tensegrity
eigenvals, _ = RB.undamped_eigen(bridge_gdr.structure, gravity)
println("Natural frequencies (first 5): $(eigenvals[1:min(5, length(eigenvals))])")
```

## Dynamics Simulation

Tensegrity dynamics are handled by the standard Rible solver pipeline, with special consideration for contact modeling since tensegrity robots often roll on or interact with ground surfaces.

For a complete walkthrough of dynamics simulation with a rolling superball tensegrity robot, see [Quick Start](@ref "Rible.jl Quick Start").

## Stiffness and Stability Analysis

### Static-Kinematic Determinacy

`static_kinematic_determine` performs SVD on the transpose of the equilibrium matrix to decompose the structure into **self-stress states** and **mechanism modes**:

```@example tensegrity
# Build a superball for stiffness analysis
ballbot_s = superball(0.0;
    θ = 0.0,
    l = 2.0 / 2,
    d = 2.0 / 4,
    z0 = 2.0 / 2,
    visible = false,
    loadmesh = false,
    constrained = false,
)
st_s = ballbot_s.structure

# Verify static equilibrium and get cable tensions
RT.check_static_equilibrium_output_multipliers(st_s, nofield)
RB.update!(st_s)
f = RT.get_cables_tension(ballbot_s)

# Build the equilibrium matrix Bᵀ
Q̃ = RT.build_Q(st_s)
L̂ = RT.build_L̂(st_s)
Bᵀ = -Q̃ * L̂

# Project into the free-DOF nullspace
Ǎ = RB.cstr_jacobian(st_s, st_s.state.system)
Ň = RB.nullspace(Ǎ)
ℬᵀ = transpose(Ň) * Bᵀ

# Decompose into self-stress and mechanism modes
S, D = RT.static_kinematic_determine(ℬᵀ)
println("Self-stress states: $(size(S, 2))")
println("Mechanism modes:    $(size(D, 2))")
```

- **Self-stress states** (columns of `S`): internal force vectors that balance without external forces.
- **Mechanism modes** (columns of `D`): directions in which the structure can displace without stretching cables.

### Material and Geometric Stiffness

The tangent stiffness matrix ``K_T`` is key to assessing structural stability:

```math
K_T = K_{\text{material}} + K_{\text{geometric}} + K_{\text{constraint}}
```

Material stiffness ``K_{\text{material}}`` reflects elastic resistance to cable stretching. Both geometric stiffness ``K_{\text{geometric}}`` and constraint stiffness ``K_{\text{constraint}}`` stem from geometric pertubation and are affected by cables' prestress level, while only the ``K_{\text{constraint}}`` is found to be the sole source of instability.

```@example tensegrity
q = RB.get_coords(st_s)
k = RT.get_cables_stiffness(st_s)

# Material stiffness matrix
Ǩm = RT.build_material_stiffness_matrix!(st_s, q, k)

# Geometric stiffness matrix
Ǩg = RT.build_geometric_stiffness_matrix!(st_s, q, f)

# Project to free-coordinate space
𝒦m = transpose(Ň) * Ǩm * Ň |> Symmetric
𝒦g = transpose(Ň) * Ǩg * Ň |> Symmetric

vals_m = sort(eigvals(𝒦m))
println("Material stiffness eigenvalues (min): $(vals_m[1])")
```

A sufficient condition for stability is that the total tangent stiffness ``K_T`` is positive definite in the mechanism directions.

## Prestress Optimization

With convex optimization solvers (COSMO, Clarabel), RibleTensegrity can search for optimal performance points in the self-stress space [luoStabilityConditionsStiffness2024](@cite):

- **`optimize_maximum_stiffness`**: maximize the minimum eigenvalue of the stiffness matrix, improving resistance to external loads.
- **`optimize_zero_stiffness`**: locate the critical prestress level where the minimum eigenvalue crosses zero, identifying the onset of instability.

### Example: Maximizing Stiffness

```@example tensegrity
# Build prestress stiffness contributions for each self-stress state
ns = size(S, 2)
vec𝒦ps = map(1:ns) do i
    si = S[:, i]
    λi = inv(Ǎ * transpose(Ǎ)) * Ǎ * Bᵀ * si
    Ǩai = -RB.cstr_forces_jacobian(st_s, q, λi)
    Ǩgi = RT.build_geometric_stiffness_matrix!(st_s, q, si)
    𝒦pi = transpose(Ň) * (Ǩgi .+ Ǩai) * Ň |> Symmetric
    vec(𝒦pi)
end
mat𝒦ps = reduce(hcat, vec𝒦ps)

# Set up the optimization problem
vec𝒦m = vec(𝒦m)
vecI = vec(Matrix(1.0I, size(𝒦m)))
ᾱ = ones(ns)                         # self-stress combination weights
nx = ns + 2                           # decision variables: [α..., σ, ρ]

A = hcat(-Matrix(1.0I, ns, ns), ᾱ, zero(ᾱ))
b = zeros(ns)

# Maximize minimum stiffness (ρ)
result = RT.optimize_maximum_stiffness(mat𝒦ps, vec𝒦m, vecI, A, b, nx)
σ_opt = result.x[end-1]               # optimal prestress scale
ρ_opt = result.x[end]                 # maximum minimum eigenvalue

println("Optimal prestress scale: $σ_opt")
println("Maximum min eigenvalue:  $ρ_opt")
```

### Finding Zero-Stiffness Points

```@example tensegrity
# Search for the critical prestress where minimum stiffness → 0
result_zero = RT.optimize_zero_stiffness(
    mat𝒦ps, vec𝒦m, vecI,
    hcat(-Matrix(1.0I, ns, ns), ᾱ),   # equality constraints
    zeros(ns),                          # constraint RHS
    ns + 1,                             # number of decision variables
    result.x[1:end-1],                  # warm-start from max-stiffness result
)

σ_zero = result_zero.x[end]
println("Zero-stiffness prestress scale: $σ_zero")

# Verify: the stiffness matrix at this scale should have a near-zero eigenvalue
𝒦_zero = 𝒦m + σ_zero * reshape(mat𝒦ps * ᾱ, size(𝒦m))
ρ_zero = minimum(eigvals(Symmetric(𝒦_zero)))
println("Min eigenvalue at zero-stiffness point: $ρ_zero")
```

!!! tip "Solvers and warm-starting"
    `optimize_maximum_stiffness` and `optimize_zero_stiffness` use COSMO internally. The zero-stiffness search accepts a warm-start vector (`x_0`) — using the solution from `optimize_maximum_stiffness` as the starting point significantly improves convergence.

## References

```@bibliography
Pages = ["RibleTensegrity.md"]
```
