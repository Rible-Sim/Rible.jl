function slider_crank(;θ = 0, coordsType = RB.NCF.NC)
    SVo3 = SVector{3}([0.0,0.0,0.0])
    l = [0.153,0.306]
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
    ω0 = [
        150.0,
        -75.0,
    ]./10
    b = 0.05
    a = 0.025
    d = 0.05
    # θ0 = [
    #     deg2rad(  0),
    #     deg2rad( 45),
    #     deg2rad(-30),
    #     deg2rad( 15)
    # ]
    θ0 = zeros(4)
    function make_base(i)
        contactable = false
        visible = true
        r̄g = SVo3
        r̄p1 = SVector{3}([0.0, 0,  0])
        loci_positions = [r̄p1,]
        axes_normals = [
            SVector{3}([1.0,0.0,0.0]),
        ]
        m = mass[1]
        Ī = SMatrix{3,3}([
            J[1] 0 0;
            0 J[1] 0;
            0 0 J[1];
        ])
        R = RotX(θ0[i])
        ω = SVo3
        ri = SVo3
        ro = ri
        ṙo = zero(ro)

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
    function make_link(i;
            ro = SVo3
        )
        contactable = false
        visible = true
        r̄p1 = SVector{3}([0.0, -l[i-1]/2, 0])
        r̄p2 = SVector{3}([0.0,  l[i-1]/2, 0])
        r̄g  = SVector{3}([0.0, 0,     0])
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
        R = RotX(θ0[i])
        ω = SVector{3}(
            ω0[i-1],0,0,
        )
        ri = ro
        # @show ro
        if i == 2
            ṙo = [0,0,ω0[i-1]*l[i-1]/2]
        else
            ṙo = [0,0,ω0[i-2]*l[i-2] .+ ω0[i-1]*l[i-1]/2]
        end

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
                ifelse(
                    i == 2,
                    0.01*2,
                    0.01
                ),
                -l[i-1]/2,
                0
            ],
            # rot = RotZ(π),
            scale=1/1000,
            color = :slategrey,
        )

        RB.RigidBody(prop,state,coords,linkmesh)
    end
    function make_slider(i;
            ro=SVo3,
            R=RotZ(0.0),
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
        R = RotX(θ0[i])
        ri = ro
        ṙo = zero(ro)
        ω = zero(ro)
        
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
    base = make_base(1)
    d1 = SVector(0,l[1]*cos(θ0[2]),l[1]*sin(θ0[2]),)
    p1 = d1./2
    link1 = make_link(2;ro=p1)
    d2 = SVector(0,l[2]*cos(θ0[3]),l[2]*sin(θ0[3]),)
    p2 = d1 + d2./2
    link2 = make_link(3;ro=p2)
    p3 = d1 + d2
    slider1 = make_slider(4;ro = p3)
    rbs = [base,link1,link2,slider1,]
    rigdibodies = TypeSortedCollection(rbs)
    numbered = RB.number(rigdibodies)
    indexed = RB.index(rigdibodies,)

    ss = Int[]
    tensiles = (cables = ss,)
    connected = RB.connect(rbs,zeros(Int,0,0))
    tensioned = @eponymtuple(connected,)

    j1 = RB.FixedBodyConstraint(1,indexed,base)
    j2 = RB.RevoluteJoint(2,indexed,RB.Hen2Egg(2,RB.ID(base ,1,1),RB.ID(link1,1,1)))
    j3 = RB.RevoluteJoint(3,indexed,RB.Hen2Egg(3,RB.ID(link1,2,2),RB.ID(link2,1,1)))
    j4 = RB.RevoluteJoint(4,indexed,RB.Hen2Egg(4,RB.ID(link2,2,2),RB.ID(slider1,5,5)))

    js = [j1,j2,j3,j4]

    jointed = RB.join(js,indexed)
    cnt = RB.Connectivity(numbered,indexed,tensioned,jointed)
    st = RB.Structure(rigdibodies,tensiles,cnt)
    RB.Robot(st,)
end