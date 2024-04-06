
function planar_parallel()
    SVo3 = SVector{3}([0.0,0.0,0.0])
    a = 1.0
    b = 0.5
    θ = 2π/3
    function make_base(i)
        contactable = true
        visible = true
        mass_center_position = SVo3
        r̄p1 = SVector{3}([0.0, a,0.0])
        r̄p2 = RotZ( θ)*r̄p1
        r̄p3 = RotZ(-θ)*r̄p1
        nodes_positions = [r̄p1,r̄p2,r̄p3,]
        axes = [
            SVector{3}([1.0,0.0,0.0])
            for i in eachindex(nodes_positions)
        ]
        m = 1.0
        Ī = SMatrix{3,3}(Matrix(1.0I,3,3))
        ci = collect(1:12)
        cstr_idx = Int[]    
        R = RotZ(0.0)    
        ω = SVo3
        ri = SVo3
        ro = ri
        ṙo = zero(ro)

        prop = RB.RigidBodyProperty(
                    i,contactable,m,
                    Ī,
                    mass_center_position,
                    nodes_positions,
                    axes;
                    visible
                    )

        nmcs = RB.NCF.NC1P3V(ri, ro, R)

        state = RB.RigidBodyState(prop, nmcs, ri, R, ṙo, ω, ci, cstr_idx)

        RB.RigidBody(prop,state)
    end
    function make_platform(i;ri=SVo3,R=RotZ(0.0),)
        contactable = true
        visible = true
        mass_center_position = SVo3
        r̄p1 = SVector{3}([0.0, b,0.0])
        r̄p2 = RotZ( θ)*r̄p1
        r̄p3 = RotZ(-θ)*r̄p1
        nodes_positions = [r̄p1,r̄p2,r̄p3]
        axes = [SVector{3}([1.0,0.0,0.0]),]
        Ī = SMatrix{3,3}(Matrix(1.0I,3,3))
        ω = SVo3
        ro = ri
        ṙo = zero(ro)
        m = 1.0

        prop = RB.RigidBodyProperty(
                    i,contactable,m,
                    Ī,
                    mass_center_position,
                    nodes_positions,
                    axes;
                    visible
                    )

        nmcs = RB.NCF.NC1P3V(ri, ro, R)
        state = RB.RigidBodyState(prop, nmcs, ri, R)

        RB.RigidBody(prop,state)
    end
    base = make_base(1)
    platform = make_platform(2;ri = SVector(0.0,0.0,0.0))
    # slider2 = make_platform(3;ri = SVo3,)
    # bar = make_3d_bar(
    #     4,
    #     SVector(-a,0.0,0.0),
    #     SVo3;
    # )
    rbs = [base,platform]
    rigdibodies = TypeSortedCollection(rbs)
    numbered = RB.number(rigdibodies)
    # sm = [
    #     0 1 0 1;
    #     0 2 0 2;
    #     0 3 0 3;
    #     0 0 1 4;
    #     0 0 2 5;
    #     0 0 3 6;
    # ]
    # indexed = RB.index(rigdibodies,sm)
    indexed = RB.index(rigdibodies,)

    ncables = 9
    original_restlens = zeros(ncables)
    original_restlens = zeros(ncables)
    ks = fill(100.0,ncables)
    cs = zeros(ncables)
    ss = [RB.DistanceSpringDamper3D( original_restlens[i],ks[i],cs[i];slack=false) for i = 1:ncables]
    apparatuses = (cables=ss,)

    cm = [
        1 -1 ;
        2 -2 ;
        3 -3 ;
        1 -2 ;
        1 -3 ;
        2 -1 ;
        2 -3 ;
        3 -1 ;
        3 -2 ;
    ]

    connected = RB.connect(rigdibodies, cm)
    tensioned = @eponymtuple(connected,)

    # j1 = RB.PrismaticJoint(1,RB.Hen2Egg(RB.Signifier(base,5,1),RB.Signifier(slider1,1,1)))
    # j2 = RB.PrismaticJoint(2,RB.Hen2Egg(RB.Signifier(base,5,2),RB.Signifier(slider2,1,1)))

    # js = [
    #     j1,j2
    # ]
    # jointed = RB.join(js,indexed)
    cnt = RB.Connectivity(numbered,indexed,tensioned,)
    st = RB.Structure(rigdibodies,apparatuses,cnt)
    RB.Robot(st,)
end