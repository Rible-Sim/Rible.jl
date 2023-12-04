# Body & Coordinates
Immutable properties, and coords type informations.
passed around with no allocations.

Mutable states, and coords cache.

Typical calling methods: use the properties and coords type to mutate state and cache.

Question: what about joints? can it be fully immutable?
Or immutable/mutable pairs?

# Heterogeneous collection 

type-stable code unavoidably requires `generated` functions.

# Modelling

aim for closed-loop kinematics, not open loop

performance
  no recursive algorithms
  parallel algorithms

aim for parametric optimization

excessive data: 
more than necessary, but not redundant

# Docs
run `make.jl` in the `yard` environment, where developmental/experimental codes reside.

run Literate to transform dev codes to `@example` code

run `makedocs` to actually run `@example` codes and generates docs.

Examples are reused for dev/docs/test/benchmark.
To do this, use `include`, docs/test/benchmark; use `includet` dev;
- [ ] where to put assets?
  - [ ] in robots too?
  - [ ] already independent of `.jl` files
- [ ] where to put `vis.jl`?
  - [ ] something like `deps.jl`

# Todo

- [ ] rename fields
  - [ ] less symbols, more words
  - [ ] use snake case
- [ ] connectivity, mask, use `BitVectors` isstead of `Vector{Int}` for free/pres, indice, 
- [ ] more preallocations/cache
    - [ ] free-flights
      - [ ] rigid bodies
        - [ ] natural coords
          - [x] get rid of functions fields
          - [x] where to put `pres_idx`, `free_idx`, `cstr_idx`? To `NonminimalCoordinates` for now
          - [x] how to cache hessians? to get both `∂Aq̇∂q` and `cstr_forces_jacobian`
        - [ ] quaternion-based
          - [x] differentiate position cache and velocity cache (acc cache not need for now) to_velocity_transformation is always at velocity level
          - [ ] similarly `∂M⁻¹p∂q` and `∂Mq̇∂q`
          - [ ] merge mass matrices cache and `∂T∂qᵀ`
          - [ ] use `∂T∂xᵀ∂x` and  `∂T∂xᵀ∂ẋ`
          - [ ] split `F!` into inertia (merged with mass matrices cache) and body internal strain and disspation and intra-bodies tensional forces and joints forces.
            - [ ] order matter! 
            - [ ] `Zhong06` familiy all depend on position coordinates
      - [ ] flexible bodies
        - [ ] absolute nodal coords
  - [x] integrator cache, used for all timesteps
    - [x] timestep cache, contacts related, because the num of contacts no known in advance,
      - [x] use solver / contact model / for dispatch
    - [x] Zhong06 family
    - [ ] Alpha family
    - [x] friction model: frictionless, Coulomb, PolyhedralCoulomb
    - [x] resitution: Inelastic, Newton, Poisson, Strange
    - [x] variants: Unclassified/Classified CCP, Mono/Two-layer,
    - [ ] Moreau-$\theta$ scheme
- [x] reduce functions fields
- [x] preallocations for autodiff: `QCF`  cstr_forces_jacobian
- [ ] joint cstr
  - [x] quadratic forms
  - [x] put joint into `NCF`, `QCF`
  - [x] reorder basic joint numbering
  - [x] reorder normal, tangent, bitangent
  - [x] `QCF` joints constraints hessians -> cstr_forces_jacobian
- [x] put contacts into loci
- [ ] local distribution law
- [ ] try SNAKE/RATTLE scheme
- [ ] visualize Locus/LocusState
- [ ] investigate scaling, may be scale relative velocity to position?
- [ ] A (moving) frame type for (origin_position,R,origin_velocity,ω) ? quaternion?
- [ ] ReinforcementLearning
  - [ ] value/reward will always depends on position variables: final state or intermediate states
  - [ ] uppon convergence, an additional Newton iteration is performed to compute all the derivatives/gradients/jacobians/adjoint and cached
  - [ ] cache differentials, w.r.t action(forces/torques), structual parameters (stiffness/geometry), joints(location, axes), initial conditions, 

# Show case
- [ ] pecking wood bird examples
- [ ] 走钢丝