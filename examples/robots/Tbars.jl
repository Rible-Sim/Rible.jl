# includet("../bodies/make_3d_bar.jl")
function Tbars(;θ = 0, coordsType = RB.NCF.NC)
    SVo3 = SVector{3}([0.0,0.0,0.0])
    a = 1.0
    b = 1.0
    c = 0.05
    function make_base(i)
        contactable = true
        visible = true
        mass_center_position = SVo3
        r̄p1 = SVector{3}([ -a, b,c])
        r̄p2 = SVector{3}([ -a,-b,c])
        r̄p3 = SVector{3}([0.0, b,c])
        r̄p4 = SVector{3}([0.0,-b,c])
        r̄p5 = SVector{3}([0.0,0.0,c])
        r̄p6 = SVector{3}([-2a,0.0,c])
        r̄p7 = SVector{3}([  a,0.0,c])
        nodes_positions = [r̄p1,r̄p2,r̄p3,r̄p4,r̄p5,r̄p6,r̄p7]
        axes = [
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([0.0,1.0,0.0]),
            SVector{3}([0.0,1.0,0.0]),
            SVector{3}([0.0,1.0,0.0]),
            SVector{3}([0.0,1.0,0.0]),
            SVector{3}([0.0,1.0,0.0]),
            SVector{3}([0.0,1.0,0.0]),
        ]
        m = 1.0
        Ī = SMatrix{3,3}(Matrix(1.0I,3,3))
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

        state = RB.RigidBodyState(prop, ri, R, ṙo, ω)
        if coordsType isa Type{RB.NCF.NC}
            nmcs = RB.NCF.NC1P3V(ri, ro, R)
            ci = collect(1:12)
            cstr_idx = Int[]    
            coords = RB.NonminimalCoordinates(nmcs, ci, cstr_idx)
        else
            qcs = RB.QCF.QC(m,Ī)
            ci = collect(1:7)
            cstr_idx = Int[]    
            coords = RB.NonminimalCoordinates(qcs, ci, cstr_idx)
        end
        basemesh = load(RB.assetpath("装配体1.STL")) |> RB.make_patch(;
            # trans=[-1.0,0,0],
            rot = RotZ(π),
            scale=1/500,
            color = :silver,
        )

        RB.RigidBody(prop,state,coords,basemesh)
    end
    function make_slider(i;ri=SVo3,R=RotZ(0.0),)
        contactable = true
        visible = true
        mass_center_position = SVo3
        r̄p1 = SVector(0.0,0.0,c)
        r̄p2 = SVector(0.0,0.0,c/2)
        nodes_positions = [r̄p1,r̄p2]
        axes = [
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
        ]
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

        # @show q[1:3]
        # @show q[4:6]
        # @show q[7:9]
        # @show q[10:12]
        ro = ri
        ṙo = zero(ro)
        ω = zero(ro)
        state = RB.RigidBodyState(prop, ri, R, ṙo, ω)
        if coordsType isa Type{RB.NCF.NC}
            nmcs = RB.NCF.NC1P3V(ri, ro, R)
            coords = RB.NonminimalCoordinates(nmcs,)
        else
            qcs = RB.QCF.QC(m,Ī)
            coords = RB.NonminimalCoordinates(qcs,)
        end
        slidermesh = load(RB.assetpath("装配体2.2.STL")) |> RB.make_patch(;
            # trans=[-1.0,0,0],
            rot = begin
                if i == 2
                    RotYZ(π,π/2)
                else
                    RotY(π)
                end
            end,
            scale=1/500,
            color=:mediumpurple4
        )
        RB.RigidBody(prop,state,coords,slidermesh)
    end
    base = make_base(1)
    p1 = SVector(-b*cos(θ),0.0,0)
    p2 = SVector(0,b*sin(θ),0)
    slider1 = make_slider(2;ri = p1)
    slider2 = make_slider(3;ri = p2,)
    bar = make_3d_bar(
        4,
        p1 + SVector(0,0,c/2),
        p2 + SVector(0,0,c/2);
		loadmesh = true,
    )
    rbs = [base,slider1,slider2,bar]
    rigdibodies = TypeSortedCollection(rbs)

    ncables = 4
    original_restlens = zeros(ncables)
    original_restlens = zeros(ncables)
    ks = fill(100.0,ncables)
    ks[2] = 400.0
    cs = zeros(ncables)
    spring_dampers = [RB.DistanceSpringDamper3D(original_restlens[i],ks[i],cs[i];slack=false) for i = 1:ncables]

    connecting_matrix = [
        1 -1  0  0;
        2 -1  0  0;
        # 3  0 -1  0;
        # 4  0 -1  0;
        6 -1  0  0;
        7 -1  0  0;
    ]

    cables = RB.connect(rigdibodies, spring_dampers; connecting_matrix,)


    j1 = RB.PrismaticJoint(ncables+1,RB.Hen2Egg(RB.Signifier(base,5,1),RB.Signifier(slider1,1,1)))
    j2 = RB.PrismaticJoint(ncables+2,RB.Hen2Egg(RB.Signifier(base,5,2),RB.Signifier(slider2,1,1)))

    j3 = RB.PinJoint(ncables+3,RB.Hen2Egg(RB.Signifier(bar,1,1),RB.Signifier(slider1,2,1)))
    j4 = RB.PinJoint(ncables+4,RB.Hen2Egg(RB.Signifier(bar,2,1),RB.Signifier(slider2,2,1)))

    js = [j1,j2,j3,j4]

    apparatuses = TypeSortedCollection(
        vcat(cables,js)
    )
    
    # sm = [
    #     0 1 0 1;
    #     0 2 0 2;
    #     0 3 0 3;
    #     0 0 1 4;
    #     0 0 2 5;
    #     0 0 3 6;
    # ]
    # indexed = RB.index(rigdibodies,sm)
    st = RB.Structure(rigdibodies,apparatuses,)
    gauges = Int[]
    actuators = Int[]
    hub = RB.ControlHub(
        st,
        gauges,
        actuators,
        RB.Coalition(st,gauges,actuators)
    )

    RB.Robot(st,hub)
end