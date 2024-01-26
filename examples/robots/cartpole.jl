function cart_pole(;
        y0 =  0.0,
        ẏ0 =  0.0,
        θ0 =  π/2,
        ω0 =  0.0,
        coordsType = RB.NCF.NC,
        f = (t)->[1.0],
    )
    SVo3 = SVector{3}([0.0,0.0,0.0])
    b = 0.02
    a = 0.01
    l = 1.0
    mass = [
        0.038,
        0.038,
        0.076
    ]
    J = [
        7.4e-5,
        5.9e-4,
        2.7e-6
    ]
    function make_base(i,
            R = RotX(0.0),
            ω = SVo3,
            ro = SVo3,
            ṙo = zero(ro)
        )
        contactable = false
        visible = true
        r̄g = SVo3
        r̄p1 = SVector{3}([0.0, 0,  0])
        loci_positions = [r̄p1,]
        axes_normals = [
            SVector{3}([0.0,1.0,0.0]),
        ]
        m = mass[1]
        Ī = SMatrix{3,3}([
            J[1] 0 0;
            0 J[1] 0;
            0 0 J[1];
        ])
        ri = ro

        prop = RB.RigidBodyProperty(
            i,contactable,m,
            Ī,
            r̄g,
            loci_positions,
            axes_normals;
            visible
        )

        state = RB.RigidBodyState(prop, ri, R, ṙo, ω)
        if coordsType isa Type{RB.NCF.NC}
            nmcs = RB.NCF.NC1P3V(ri, ro, R)
            pres_idx = Int[]
            cstr_idx = collect(1:6)
            coords = RB.NonminimalCoordinates(nmcs, pres_idx, cstr_idx)
        else
            qcs = RB.QCF.QC(m,Ī)
            pres_idx = Int[]
            cstr_idx = [1]
            coords = RB.NonminimalCoordinates(qcs, pres_idx, cstr_idx)
        end

        basemesh = load(RB.assetpath("crank_slider/base.STL")) |> RB.make_patch(;
            # trans=[-1.0,0,0],
            # rot = RotZ(π),
            scale=1/1000,
            color = :silver,
        )

        RB.RigidBody(prop,state,coords,basemesh)
    end
    function make_cart(i;
            ro=SVo3,
            R = RotX(0.0),
            ṙo = zero(ro),
            ω = zero(ro),
        )
        contactable = true
        visible = true
        r̄g = SVo3
        r̄p1 = SVector(0.0,-b,-a)
        r̄p2 = SVector(0.0, b,-a)
        r̄p3 = SVector(0.0, b, a)
        r̄p4 = SVector(0.0,-b, a)
        r̄p5 = SVector(0.0, 0, 0)
        loci_positions = [r̄p1,r̄p2,r̄p3,r̄p4,r̄p5]
        axes_normals = [
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
        ]
        m = mass[i-1]
        Ī = SMatrix{3,3}([
            J[i-1] 0 0;
            0 J[i-1] 0;
            0 0 J[i-1];
        ])
        
        prop = RB.RigidBodyProperty(
            i,contactable,m,
            Ī,
            r̄g,
            loci_positions,
            axes_normals,
            fill(0.5,length(loci_positions)),
            fill(0.4,length(loci_positions));
            visible,
        )

        ri = ro
        # @show q[1:3]
        # @show q[4:6]
        # @show q[7:9]
        # @show q[10:12]
        state = RB.RigidBodyState(prop, ri, R, ṙo, ω)
        if coordsType isa Type{RB.NCF.NC}
            nmcs = RB.NCF.NC1P3V(ri, ro, R)
            coords = RB.NonminimalCoordinates(nmcs,)
        else
            qcs = RB.QCF.QC(m,Ī)
            coords = RB.NonminimalCoordinates(qcs,)
        end
        slidermesh = load(RB.assetpath("crank_slider/slider.STL")) |> RB.make_patch(;
            trans=[0.0,0,0],
            scale=1/1000,
            color=:mediumpurple4
        )
        RB.RigidBody(prop,state,coords,slidermesh)
    end
    function make_pole(i;
            ro = SVo3,
            R=RotX(θ0),
            ω = SVector{3}(
                ω0,0,0,
            ),
            ṙo = [0,0,ω0*l/2]
        )
        contactable = false
        visible = true
        r̄p1 = SVector{3}([0.0, -l/2, 0])
        r̄p2 = SVector{3}([0.0,  l/2, 0])
        r̄g  = SVector{3}([0.0,  0.0, 0])
        loci_positions = [r̄p1,r̄p2,]
        axes_normals = [
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
        ]
        m = mass[i-1]
        Ī = SMatrix{3,3}([
            J[i-1] 0 0;
            0 J[i-1] 0;
            0 0 J[i-1];
        ])
        ri = ro
        # @show ro
        prop = RB.RigidBodyProperty(
            i,contactable,m,
            Ī,
            r̄g,
            loci_positions,
            axes_normals;
            visible
        )

        state = RB.RigidBodyState(prop, ro, R, ṙo, ω)
        if coordsType isa Type{RB.NCF.NC}
            nmcs = RB.NCF.NC1P3V(ri, ro, R)
            pres_idx = Int[] 
            cstr_idx = collect(1:6)    
            coords = RB.NonminimalCoordinates(nmcs, pres_idx, cstr_idx)
        else
            qcs = RB.QCF.QC(m,Ī)
            pres_idx = Int[]
            cstr_idx = [1]
            coords = RB.NonminimalCoordinates(qcs, pres_idx, cstr_idx)
        end

        linkmesh = load(RB.assetpath("crank_slider/crank$(i-1).STL")) |> RB.make_patch(;
            trans=[
                0.01,
                -l/2,
                0
            ],
            # rot = RotZ(π),
            scale=1/1000,
            color = :slategrey,
        )

        RB.RigidBody(prop,state,coords,linkmesh)
    end
    base = make_base(1)
    p1 = [0.0,y0,0.0]
    p2 = [0.0,y0,0.0]
    cart = make_cart(2;
        ro=p1,
        ṙo = [0,0,ω0*l/2]
    )
    pole = make_pole(3;
        ro=p2+RotX(θ0)*[0.0, l/2, 0],
        R=RotX(θ0),
        ω = SVector{3}(
            ω0,0,0,
        ),
        ṙo = [0,0,ω0*l/2]
    )
    rbs = [base,cart,pole]
    rigdibodies = TypeSortedCollection(rbs)

    j1 = RB.FixedBodyConstraint(1,base)
    j2 = RB.PrismaticJoint(2,RB.Hen2Egg(RB.Signifier(base,1,1),RB.Signifier(cart,5,5)))
    j3 = RB.RevoluteJoint(3,RB.Hen2Egg(RB.Signifier(cart,5,5),RB.Signifier(pole,1,1)))

    apparatuses = TypeSortedCollection([j1,j2,j3])

    numbered = RB.number(rigdibodies,apparatuses)
    indexed = RB.index(rigdibodies,apparatuses)
    cnt = RB.Connectivity(numbered,indexed)
    structure = RB.Structure(rigdibodies,apparatuses,cnt)

    pos_vel_errors = RB.ErrorGauge(
        1,
        RB.Signifier(cart,5,5),
        RB.PositionCapta(),
        zeros(3)
    )
    ## angle_err = RB.ErrorGauge(
    ##     RB.Gauge(
    ##         RB.Signifier(j3,3),
    ##         RB.AngleCapta()
    ##     ),
    ##     π/2
    ## )
    gauges = TypeSortedCollection([pos_vel_errors,])
    force_actuator = RB.ExternalForceActuator(
        1,
        RB.Signifier(cart,5,5),
        RB.TimeOperator(f),
        [0,1.0,0],
        [0.0],
    )
    actuators = TypeSortedCollection([force_actuator,])
    hub = RB.ControlHub(
        structure,
        gauges,
        actuators,
        RB.Coalition(
            structure,gauges,actuators
        )
    )
    RB.Robot(structure,hub)
end