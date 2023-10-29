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
        - [ ] quaternion-based
      - [ ] flexible bodies
        - [ ] absolute nodal coordinates
    - [ ] ZhongCCP
    - [ ] ZhongQCCP
- [ ] reduce functions fields
- [ ] preallocations for autodiff
- [ ] joint constraints, quadratic forms
- [x] put contacts into loci
- [ ] local distribution law
- [ ] try SNAKE/RATTLE scheme
- [ ] visualize Locus/LocusState
- [ ] investigate scaling, may be scale relative velocity to position?