function new_pointmass(;
        origin_position = [0.0,0.0,1.0],
        origin_velocity = zero(origin_position),
        m = 1.0,
        μ = 0.3,
        e = 0.9
    )
    contactable = true
    visible = true
    Ia = SMatrix{3,3}(Matrix(m*I,3,3))
    mass_locus  = SVector{3}([ 0.0, 0.0, 0.0])
    r̄p1 = SVector{3}([ 0.0, 0.0, 0.0])
    loci = [r̄p1]
    axes = [SVector{3}([ 1.0, 0.0, 0.0])]
    friction_coefficients = [μ]
    restitution_coefficients = [e]
    prop = RB.RigidBodyProperty(
        1,contactable,m,Ia,
        mass_locus,loci,axes,
        friction_coefficients,restitution_coefficients,
        ;visible=visible
    )
    ω = zero(origin_position)
    R = RotX(0.0)
    loci = Ref(origin_position) .+ Ref(R).*loci
    nmcs = RB.NCF.NC3D1P(loci[1],)
    ci = Int[]
    cstr_idx = Int[]
    state = RB.RigidBodyState(prop,origin_position,R,origin_velocity,ω)
    coords = RB.NonminimalCoordinates(nmcs,ci,cstr_idx)
    rb1 = RB.RigidBody(prop,state,coords)

    rbs = TypeSortedCollection([rb1,])
    apparatuses = Int[]
    indexed = RB.index(rbs,apparatuses)
    numbered = RB.number(rbs,apparatuses)
    cnt = RB.Connectivity(indexed,numbered,)
    st = RB.Structure(rbs,apparatuses,cnt,)


    gauges = Int[]
    actuators = Int[]
    ## actuators = [
    ##     RB.RegisterActuator(
    ##         1,
    ##         collect(1:ncables_prism),
    ##         zeros(ncables_prism)
    ##     ),
    ## ]
    hub = RB.ControlHub(
        st,
        gauges,
        actuators,
        RB.Coalition(st,gauges,actuators)
    )
    bot = RB.Robot(st,hub)
end
