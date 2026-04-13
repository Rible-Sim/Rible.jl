

function pointmass3d(id = 1; 
        origin_position = [0.0,0.0,1.0],
        origin_velocity = zero(origin_position),
        m = 1.0,
        μ = 0.3,
        e = 0.9,
        contactable = true,
        visible = true,
    )
    Ia = SMatrix{3,3}(Matrix(m*I,3,3))
    mass_locus  = SVector{3}([ 0.0, 0.0, 0.0])
    r̄p1 = SVector{3}([ 0.0, 0.0, 0.0])
    loci = [r̄p1]
    axes = [SVector{3}([ 1.0, 0.0, 0.0])]
    friction_coefficients = [μ]
    restitution_coefficients = [e]
    prop = RB.RigidBodyProperty(id, contactable, m, Ia, RB.Locus(mass_locus), [RB.Locus(pos, ax, mu, e) for (pos,ax,mu,e) in zip(loci, axes, friction_coefficients, restitution_coefficients)])
    ω = zero(origin_position)
    R = RotX(0.0)
    coords = RB.NCF.NC3D1P(loci[1],)
    state = RB.RigidBodyState(prop,origin_position,R,origin_velocity,ω)
    RB.RigidBody(prop,state,coords)
end

function new_pointmass(;
        origin_position = [0.0,0.0,1.0],
        origin_velocity = zero(origin_position),
        term_position = [2.5,0.0,1.0],
        m = 1.0,
        μ = 0.3,
        e = 0.9,
        external_force = [1 0 ; 0 1 ; 0 0 ],
    )
    rb1 = pointmass3d(
        1;
        origin_position,
        origin_velocity,
        m,
        μ,
        e,
    )

    rbs = [rb1,]
    apparatuses = Int[]
    cnt = RB.Connectivity(rbs,apparatuses,)
    st = RB.Structure(rbs,apparatuses,cnt,)
    
    actuators = [
        RB.ExternalForceActuator(
            1,
            RB.Signifier(rb1,1),
            RB.NaiveOperator(size(external_force,2)),
            # External force, may not be needed
            external_force,
            # Control input values
            zeros(size(external_force,2)),
        )
    ]

    pos_vel_errors = RB.ErrorGauge(
        1,
        RB.Signifier(rb1,1),
        RB.PositionCaptum(),
        term_position,
        # RB.PosVelCaptum(),
        # zeros(6),
        # RB.VelocityCaptum(),
        # zeros(3)
    )

    capta_gauges = Int[]
    error_gauges = [
        pos_vel_errors,
    ]

    hub = RB.ControlHub(
        st,
        capta_gauges,
        error_gauges,
        actuators,
        RB.Coalition(st,capta_gauges,error_gauges,actuators)
    )
    bot = RB.Robot(st,hub)
end
