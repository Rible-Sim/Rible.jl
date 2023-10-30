# Body & Coordinates
Immutable properties, and coordinates type informations.
passed around with no allocations.

Mutable states, and coordinates cache.

Typical calling methods: use the properties and coordinates type to mutate state and cache.

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

# Todo

- [ ] rename fields
  - [ ] less symbols, more words
  - [ ] use snake case
- [ ] more preallocations
    - [ ] free-flights
      - [ ] rigid bodies
        - [ ] natural coordinates
          - [ ] get rid of functions fields
          - [ ] where to put `pres_idx`, `free_idx`, `constraints_indices`?
          - [ ] how to cache hessians? to get both `∂Aq̇∂q` and `constraint_forces_jacobian`
          - [ ] similarly `∂M⁻¹p∂q` and `∂Mq̇∂q`
        - [ ] quaternion-based
      - [ ] flexible bodies
        - [ ] absolute nodal coordinates
    - [ ] ZhongCCP
    - [ ] ZhongQCCP
- [ ] reduce functions fields
- [ ] preallocations for autodiff
- [ ] joint constraints, quadratic forms
  - [ ] put joint into `NCF`, `QBF`
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