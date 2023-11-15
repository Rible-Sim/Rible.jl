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
          - [ ] how to cache hessians? to get both `∂Aq̇∂q` and `cstr_forces_jacobian`
        - [ ] quaternion-based
          - [ ] differentiate position cache and velocity cache (acc cache not need for now) to_velocity_transformation is always at velocity level
          - [ ] similarly `∂M⁻¹p∂q` and `∂Mq̇∂q`
      - [ ] flexible bodies
        - [ ] absolute nodal coords
  - [ ] integrator cache, used for all timesteps
    - [ ] timestep cache, contacts related, because the num of contacts no known in advance,
    - [ ] 
    - [ ] Zhong06
    - [ ] ZhongQ06
    - [ ] ZhongCCP
    - [ ] ZhongQCCP
- [x] reduce functions fields
- [ ] preallocations for autodiff
- [ ] joint cstr
  - [x] quadratic forms
  - [ ] put joint into `NCF`, `QCF`
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
  走钢丝