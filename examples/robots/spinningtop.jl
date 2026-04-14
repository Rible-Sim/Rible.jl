
function make_top(origin_position = [0.0,0.0,0.0],
        R = one(RotMatrix{3}),
        origin_velocity = [0.0,0.0,0.0],
        Ω = [0.0,0.0,5.0],
        cT = :NCF;
        μ = 0.5,
        e = 0.9,
        color = :slategray,
        visible = true,
        loadmesh=false,
    )
    ω = R*Ω
    contactable = true

    m =  0.58387070
    mass_locus = @SVector zeros(3)
    # Ī = SMatrix{3,3}(Matrix(1.0I,3,3))
    Ī = SMatrix{3,3}(
        Diagonal(SA[
            0.00022129,
            0.00022129,
            0.00030207
            ])
        )

    # h = 0.02292582
    radius = 0.044/√2
    h = 2*0.01897941
    loci = [radius.*SVector(1.0,1.0,0.0) for i = 1:4]
    push!(loci,SVector(0.0,0.0,-h))
    push!(loci,SVector(0.0,0.0,0.0))
    axes = [SVector(1.0,0,0) for i = 1:6]
    friction_coefficients = [μ for i = 1:6]
    restitution_coefficients = [e for i = 1:6]
    if loadmesh
        topmesh = load(
            RB.assetpath("Toupise2.STL")
        ) |> RB.make_patch(;
            scale=1/1000,
            color,
        )
    else
        pts = Point3.(loci[1:5])
        fcs = GB.TriangleFace.([
            [5,1,2],
            [5,4,3],
            [5,3,1],
            [5,2,4],
            [1,4,2],
            [4,1,3],
            [3,2,1],
            [2,3,4]
        ])
        nls = GB.normals(pts,fcs)
        topmesh = GB.Mesh(fcs; position=pts,normal=nls,)
    end
    prop = RB.RigidBodyProperty(1, contactable, m, Ī, RB.Locus(mass_locus), [RB.Locus(pos, ax, mu, e) for (pos,ax,mu,e) in zip(loci, axes, friction_coefficients, restitution_coefficients)])
    ri = origin_position#+R*loci[5]
    @debug "top" ri
    # Test if RibleQCF module is loaded
    if cT == :QCF
        nmcs = QCF.QC(m,Ī)
    elseif cT == :NCF
        nmcs = RB.NCF.NC1P3V(ri,origin_position,R)
    end
    state = RB.RigidBodyState(prop,origin_position,R,origin_velocity,ω)
    coords = nmcs
    rb1 = RB.RigidBody(prop,state,coords,topmesh)
    bodies = TypeSortedCollection([rb1,])
    apparatuses = Int[]
    cnt = RB.Connectivity(bodies,apparatuses)

    structure = RB.Structure(bodies,apparatuses,cnt,)

    pos_vel_errors = RB.ErrorGauge(
        1,
        RB.Signifier(rb1,6),
        RB.PositionCaptum(),
        [0.5,0,0.2],
        # RB.PosVelCaptum(),
        # [1.0,0,0,0,0,0],
        # RB.VelocityCaptum(),
        # zeros(3)
    )
    ## angle_err = RB.ErrorGauge(
    ##     RB.Gauge(
    ##         RB.Signifier(j3,3),
    ##         RB.AngularPositionCaptum()
    ##     ),
    ##     π/2
    ## )
    capta_gauges = Int[]
    error_gauges = TypeSortedCollection([pos_vel_errors,])
    force_actuator = RB.ExternalForceActuator(
        1,
        RB.Signifier(rb1,5),
        RB.NaiveOperator(1),
        [0,1.0,0],
        [0.0],
    )
    actuators = TypeSortedCollection([force_actuator,])
    hub = RB.ControlHub(
        structure,
        capta_gauges,
        error_gauges,
        actuators,
        RB.Coalition(
            structure,capta_gauges,error_gauges,actuators,
        ),
        # control,
    )

    bot = RB.Robot(structure, hub)
end