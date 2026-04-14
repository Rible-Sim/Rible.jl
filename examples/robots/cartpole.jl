
const SVo3 = SVector{3}([0.0,0.0,0.0])
function make_base(i;
        R = RotX(0.0),
        ω = SVo3,
        ro = SVo3,
        ṙo = zero(ro),
        m, J,
        coordsType
    )
    contactable = false
    visible = false
    r̄g = SVo3
    r̄p1 = SVector{3}([0.0, 0,  0])
    loci_positions = [r̄p1,]
    axes_normals = [
        SVector{3}([0.0,1.0,0.0]),
    ]
    Ī = SMatrix{3,3}([
        J 0 0;
        0 J 0;
        0 0 J;
    ])
    ri = ro

    prop = RB.RigidBodyProperty(i, contactable, m, Ī, 
        RB.Locus(r̄g), 
        [RB.Locus(pos, ax, mu, e) for (pos,ax,mu,e) in zip(
            loci_positions, 
            axes_normals, 
            zeros(length(loci_positions)), 
            zeros(length(loci_positions)))
        ]
    )

    state = RB.RigidBodyState(prop, ri, R, ṙo, ω)
    if coordsType == :NCF
        nmcs = RB.NCF.NC1P3V(ri, ro, R)
        coords = nmcs
    else
        qcs = QCF.QC(m,Ī)
        coords = qcs
    end

    basemesh = RB.Meshes.Box(
        RB.Meshes.Point(-0.1, -1.0, -0.05),
        RB.Meshes.Point( 0.1,  2.0,  0.0)
    ) |> RB.Meshes.boundary |> RB.simple2mesh |> RB.make_patch(;
        color = RB.CT.RGBA(0.7,0.7,0.7,0.0),
    )

    RB.RigidBody(prop,state,coords,basemesh)
end
function make_cart(i;
        ro=SVo3,
        R = RotX(0.0),
        ṙo = zero(ro),
        ω = zero(ro),
        b, a, 
        m, J,
        coordsType,
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
    Ī = SMatrix{3,3}([
        J 0 0;
        0 J 0;
        0 0 J;
    ])

    prop = RB.RigidBodyProperty(i, contactable, m, Ī, RB.Locus(r̄g), [RB.Locus(pos, ax, mu, e) for (pos,ax,mu,e) in zip(loci_positions, axes_normals, fill(0.5,length(loci_positions)), fill(0.4,length(loci_positions)))])

    ri = ro
    # @show q[1:3]
    # @show q[4:6]
    # @show q[7:9]
    # @show q[10:12]
    state = RB.RigidBodyState(prop, ri, R, ṙo, ω)
    if coordsType == :NCF
        nmcs = RB.NCF.NC1P3V(ri, ro, R)
        coords = nmcs
    else
        qcs = QCF.QC(m,Ī)
        coords = qcs
    end
    slidermesh = RB.Meshes.Box(
        RB.Meshes.Point(-a, -b, -a),
        RB.Meshes.Point( a,  b,  a)
    ) |> RB.Meshes.boundary |> RB.simple2mesh |> RB.make_patch(;
        color=:royalblue
    )
    RB.RigidBody(prop,state,coords,slidermesh)
end
function make_pole(i;
        ro = SVo3,
        R=RotX(θ0),
        ω = SVector{3}(
            ω0,0,0,
        ),
        ṙo = [0,0,ω0*l/2],
        m, J,
        l,
        coordsType,
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
    Ī = SMatrix{3,3}([
        J 0 0;
        0 J 0;
        0 0 J;
    ])
    ri = ro
    # @show ro
    prop = RB.RigidBodyProperty(i, contactable, m, Ī, 
            RB.Locus(r̄g), 
            [RB.Locus(pos, ax, mu, e) for (pos,ax,mu,e) in zip(
            loci_positions, 
            axes_normals, 
            zeros(length(loci_positions)), 
            zeros(length(loci_positions)))]
    )

    state = RB.RigidBodyState(prop, ro, R, ṙo, ω)
    if coordsType == :NCF
        nmcs = RB.NCF.NC1P3V(ri, ro, R)
        coords = nmcs
    else
        qcs = QCF.QC(m,Ī)
        coords = qcs
    end

    linkmesh = RB.endpoints2mesh(
        r̄p1,
        r̄p2;
        radius=l/20,
        color = :crimson,
    )

    RB.RigidBody(prop,state,coords,linkmesh)
end

function cart_pole(;
        y0 =  0.0,
        ẏ0 =  0.0,
        θ0 =  π/2,
        ω0 =  0.0,
        b = 0.2,
        a = 0.1,
        l = 1.0,
        coordsType = :NCF,
        f = (t)->[1.0],
        mass = [
            1.0,
            0.1,
            1.0
        ],
        J = [
            1.0,
            1.0,
            1.0
        ],
    )
    
    base = make_base(1; m = mass[1], J = J[1], coordsType)
    # cart origin
    p1 = SVector{3}(0.0, y0, 0.0)
    cart_vel = SVector{3}(0.0, ẏ0, 0.0)
    pole_rot = RotX(θ0)
    pole_offset = pole_rot * SVector{3}(0.0, l/2, 0.0)
    pole_vel = cart_vel + SVector{3}(ω0, 0.0, 0.0) × pole_offset
    cart = make_cart(2;
        ro = p1,
        ṙo = cart_vel,
        m = mass[2],
        J = J[2],
        b = b,
        a = a,
        coordsType,
    )
    # pole end?
    pole = make_pole(3;
        ro = p1 + pole_offset,
        R = pole_rot,
        ω = SVector{3}(
            ω0,0,0,
        ),
        ṙo = pole_vel,
        m = mass[3],
        J = J[3],
        l = l,
        coordsType,
    )
    rbs = [base,cart,pole]
    rigdibodies = TypeSortedCollection(rbs)

    j1 = RB.FixedBodyApparatus(1,base)
    j2 = RB.PrismaticJoint(2,RB.Hen2Egg(RB.Anchor(base,1,1),RB.Anchor(cart,5,5)))
    j3 = RB.RevoluteJoint(3,RB.Hen2Egg(RB.Anchor(cart,5,5),RB.Anchor(pole,1,1)))

    apparatuses = TypeSortedCollection([j1,j2,j3])

    
    cnt = RB.Connectivity(rigdibodies,apparatuses)
    structure = RB.Structure(rigdibodies,apparatuses,cnt)

    # fullstate_captum = RB.CaptumGauge(
    #     1, structure, RB.FullStateCaptum()
    # )
    cart_captum = RB.CaptumGauge(
        1, RB.Signifier(cart,5), RB.PosVelCaptum()
    )
    pole_captum = RB.CaptumGauge(
        2, RB.Signifier(pole,2), RB.PosVelCaptum()
    )
    capta_gauges = TypeSortedCollection([cart_captum, pole_captum])
    
    cart_pos_errors = RB.ErrorGauge(
        1,
        RB.Signifier(cart,5),
        RB.PositionCaptum(),
        zeros(3),
    )
    cart_vel_errors = RB.ErrorGauge(
        2,
        RB.Signifier(cart,5),
        RB.VelocityCaptum(),
        zeros(3),
    )
    pole_pos_errors = RB.ErrorGauge(
        3,
        RB.Signifier(pole,2),
        RB.PositionCaptum(),
        [0.0,0.0,l],
    )
    pole_vel_errors = RB.ErrorGauge(
        4,
        RB.Signifier(pole,2),
        RB.VelocityCaptum(),
        zeros(3),
    )
    cart_pos_vel_errors = RB.ErrorGauge(
        5,
        RB.Signifier(cart,5),
        RB.PosVelCaptum(),
        zeros(6),
    )
    pole_pos_vel_errors = RB.ErrorGauge(
        6,
        RB.Signifier(pole,2),
        RB.PosVelCaptum(),
        [0.0,0.0,l,0.0,0.0,0.0],
    )
    error_gauges = TypeSortedCollection([cart_pos_errors,cart_vel_errors,pole_pos_errors,pole_vel_errors,cart_pos_vel_errors,pole_pos_vel_errors])
    force_actuator = RB.ExternalForceActuator(
        1,
        RB.Signifier(cart,5),
        RB.NaiveOperator(1),
        [0;1.0;0;;],
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
    )
    RB.Robot(structure,hub)
end

"""
Compute the target frames for cart and pole based on independent states.
Returns a Dict{Int, CartesianFrame}.
"""
function set_cartpole_state!(bot, y, θ, ẏ, θ̇; 
        l #pole length
    )
    T = RB.get_numbertype(bot)
    base_frame = RB.CartesianFrame(
        (@MVector zeros(T,3)),
        (@MVector zeros(T,3)),
        RB.Axes(SMatrix{3,3,T}(I)),
        (@MVector zeros(T,3))
    )
    # Calculate desired states
    # Cart
    cart_pos = SVector{3}(0.0, y, 0.0)
    cart_vel = SVector{3}(0.0, ẏ, 0.0)
    cart_rot = RotX(0.0) # Cart doesn't rotate
    cart_ang_vel = SVector{3}(0.0, 0.0, 0.0)

    # Pole
    # Pole origin (center) relative to cart
    pole_rot = RotX(θ)
    pole_offset = pole_rot * SVector{3}(0.0, l/2, 0.0)
    
    pole_pos = cart_pos + pole_offset
    pole_ang_vel = SVector{3}(θ̇, 0.0, 0.0)
    pole_vel = cart_vel + pole_ang_vel × pole_offset

    # Construct Frames
    # Cart Frame
    cart_axes = RB.Axes(SMatrix{3,3}(cart_rot))
    cart_frame = RB.CartesianFrame(
        MVector(cart_pos),
        MVector(cart_vel),
        cart_axes,
        MVector(cart_ang_vel)
    )

    # Pole Frame
    pole_axes = RB.Axes(SMatrix{3,3}(pole_rot))
    pole_frame = RB.CartesianFrame(
        MVector(pole_pos),
        MVector(pole_vel),
        pole_axes,
        MVector(pole_ang_vel)
    )
    
    frames = [base_frame, cart_frame,pole_frame]
    RB.set_robot_state!(bot, frames)
end

