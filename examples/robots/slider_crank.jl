function slider_crank(;θ = 0, coordsType = RB.NCF.NC)
    SVo3 = SVector{3}([0.0,0.0,0.0])
    l = [0.1,0.2]
    b = 0.04
    a = 0.02
    d = 0.05
    θs0 = [
        deg2rad(  0),
        deg2rad( 45),
        deg2rad(-30),
        deg2rad( 15)
    ]
    function make_base(i)
        movable = true
        constrained = false
        r̄g = SVo3
        r̄p1 = SVector{3}([0.0, 0,  0])
        loci_positions = [r̄p1,]
        axes_normals = [
            SVector{3}([1.0,0.0,0.0]),
        ]
        m = 1.0
        Ī = SMatrix{3,3}(Matrix(1.0I,3,3))
        R = RotX(θs0[i])
        ω = SVo3
        ri = SVo3
        ro = ri
        ṙo = zero(ro)

        prop = RB.RigidBodyProperty(
            i,movable,m,
            Ī,
            r̄g,
            loci_positions,
            axes_normals;
            constrained
        )

        state = RB.RigidBodyState(prop, ri, R, ṙo, ω)
        if coordsType isa Type{RB.NCF.NC}
            nmcs = RB.NCF.NC1P3V(ri, ro, R)
            pres_idx = Int[]
            cstr_idx = collect(1:12)
            coords = RB.NonminimalCoordinates(nmcs, pres_idx, cstr_idx)
        else
            qcs = RB.QCF.QC(m,Ī)
            pres_idx = Int[]
            cstr_idx = [1]
            coords = RB.NonminimalCoordinates(qcs, pres_idx, cstr_idx)
        end

        basemesh = load(RB.assetpath("装配体1.STL")) |> RB.make_patch(;
            # trans=[-1.0,0,0],
            rot = RotZ(π),
            scale=1/500,
            color = :silver,
        )

        RB.RigidBody(prop,state,coords,basemesh)
    end
    function make_link(i;
            ro = SVo3
        )
        movable = true
        constrained = false
        r̄p1 = SVector{3}([0.0, -l[i-1]/2, 0])
        r̄p2 = SVector{3}([0.0,  l[i-1]/2, 0])
        r̄g  = SVector{3}([0.0, 0,     0])
        loci_positions = [r̄p1,r̄p2,]
        axes_normals = [
            SVector{3}([1.0,0.0,0.0]),
            SVector{3}([1.0,0.0,0.0]),
        ]
        m = 1.0
        Ī = SMatrix{3,3}(Matrix(1.0I,3,3))
        R = RotX(θs0[i])
        ω = SVo3
        ri = ro
        @show ro
        ṙo = zero(ro)

        prop = RB.RigidBodyProperty(
                    i,movable,m,
                    Ī,
                    r̄g,
                    loci_positions,
                    axes_normals;
                    constrained
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

        basemesh = load(RB.assetpath("装配体1.STL")) |> RB.make_patch(;
            # trans=[-1.0,0,0],
            rot = RotZ(π),
            scale=1/500,
            color = :silver,
        )

        RB.RigidBody(prop,state,coords,basemesh)
    end
    function make_slider(i;
            ro=SVo3,
            R=RotZ(0.0),
        )
        movable = true
        constrained = false
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
        m = 1.0
        Ī = SMatrix{3,3}(Matrix(1.0I,3,3))
        R = RotX(θs0[i])
        ri = ro
        ṙo = zero(ro)
        ω = zero(ro)
        
        prop = RB.RigidBodyProperty(
            i,movable,m,
            Ī,
            r̄g,
            loci_positions,
            axes_normals;
            constrained
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
    d1 = SVector(0,l[1]*cos(θs0[2]),l[1]*sin(θs0[2]),)
    p1 = d1./2
    link1 = make_link(2;ro=p1)
    d2 = SVector(0,l[2]*cos(θs0[3]),l[2]*sin(θs0[3]),)
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

    j1 = RB.FixedBodyConstraint(1,base,indexed)
    j2 = RB.RevoluteJoint(2,RB.Hen2Egg(2,RB.ID(base ,1,1),RB.ID(link1,1,1)))
    j3 = RB.RevoluteJoint(3,RB.Hen2Egg(3,RB.ID(link1,2,2),RB.ID(link2,1,1)))
    j4 = RB.RevoluteJoint(4,RB.Hen2Egg(4,RB.ID(link2,2,2),RB.ID(slider1,5,5)))

    js = [j1,j2,j3,j4]

    jointed = RB.join(js,indexed)
    cnt = RB.Connectivity(numbered,indexed,tensioned,jointed)
    st = RB.Structure(rigdibodies,tensiles,cnt)
    RB.Robot(st,)
end