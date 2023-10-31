function dualtri(num_of_dof,;onedir=[1.0,0.0],θ=0.0,k=400.0,c=0.0,restlen=0.16)
    nbodies = num_of_dof + 1
    nbp = 2nbodies - num_of_dof
    n_lower = count(isodd,1:nbodies)
    n_upper = count(iseven,1:nbodies)
    lower_index = 1:2:nbodies
    upper_index = 2:2:nbodies
    a = zeros(nbodies)
    m = zeros(nbodies)
    Ia = zeros(nbodies)
    for (i,j) in enumerate(lower_index)
        a[j] = 12252e-5
        m[j] = 200.40071778e-3
        Ia[j] = 3738.8e-7
    end
    for (i,k) in enumerate(upper_index)
        a[k] = 12252e-5
        m[k] = 200.40071778e-3
        Ia[k] = 3738.8e-7
    end
    R = [
        cos(θ) -sin(θ);
        sin(θ) cos(θ)
    ]
    A = zeros(2,nbp)
    A[:,2] .= A[:,1] .+ a[1]*onedir
    for i in 3:nbp
        A[:,i] .= A[:,i-1] .+ a[i-1]*R^(i-2)*onedir
    end

    function rigidbody(i,m,a,Ia,ri,α)
        movable = true
        mass_center_position = SVector{2}([0.0,0.0])

        r̄p1 = SVector{2}([-a/2,0.0])
        r̄p2 = SVector{2}([ a/2,0.0])
        if isodd(i)
            r̄p3 = SVector{2}([0.0,-√3/2*a])
            r̄p4 = r̄p3 .+ SVector{2}([ a/2,tan(deg2rad(45))])
            r̄p5 = r̄p3 .+ SVector{2}([ a/2,tan(deg2rad(45))])
        else
            r̄p3 = SVector{2}([0.0,√3/2*a])
            r̄p4 = r̄p3 .- SVector{2}([ a/2,tan(deg2rad(45))])
            r̄p5 = r̄p3 .- SVector{2}([ a/2,tan(deg2rad(45))])
        end

        nodes_positions = [r̄p1,r̄p2,r̄p3,r̄p4,r̄p5]

        axes = [SVector{2}([1.0,0.0]),]
        Ī = SMatrix{2, 2}([
            0.99Ia 0
            0 0.01Ia
        ])
        
        if i == 1
            constrained = true
            ci = collect(1:6)
            cstr_idx = Int[]
        else
            constrained = false
            ci = Int[]
            cstr_idx = collect(1:3)
        end
        
        ω = 0.0
        ro = ri
        ṙo = zero(ro)

        prop = RB.RigidBodyProperty(
                    i,movable,m,
                    Ī,
                    mass_center_position,
                    nodes_positions,
                    axes;
                    constrained=constrained
                    )

        nmcs = RB.NCF.NC1P2V(SVector{2}(ri), ro, α)

        state = RB.RigidBodyState(prop, nmcs, ri, α, ṙo, ω, ci, cstr_idx)

        body = RB.RigidBody(prop,state)
    end
    rbs = [
        begin
            d = A[:,i+1]-A[:,i]
            rigidbody(
                i,m[i],a[i],
                Ia[i],
                (A[:,i]+A[:,i+1])./2,
                atan(d[2],d[1])
            )
        end
        for i = 1:nbodies
    ]

    rigdibodies = TypeSortedCollection(rbs)
    numbered = RB.number(rigdibodies)

    indexed = RB.index(rigdibodies)

    ncables = 2(nbodies-1)
    upstringlen = restlen
    lostringlen = restlen
    original_restlens = zeros(ncables)
    restlens = zeros(ncables)
    actuallengths = zeros(ncables)
    ks = zeros(ncables)
    cs = zeros(ncables)
    cs .= c
    for i = 1:ncables
        j = i % 4
        original_restlens[i] = ifelse(j∈[1,0],upstringlen,lostringlen)
        ks[i] = ifelse(j∈[1,0],k,k)
    end
    ss = [RB.Cable2D(i, original_restlens[i],ks[i],cs[i];slack=false) for i = 1:ncables]
    tensiles = (cables=ss,)

    matrix_cnt = zeros(Int,2(nbodies-1),nbodies)
    for i = 1:nbodies-1
        ul = [1,-3]
        dl = [3,-2]
        if iseven(i)
            ul,dl = dl,ul
        end
        matrix_cnt[2(i-1)+1,i:i+1] = ul
        matrix_cnt[2(i-1)+2,i:i+1] = dl
    end

    connected = RB.connect(rigdibodies, matrix_cnt)
    tensioned = @eponymtuple(connected,)

    function ganged_act(actid,id1,id2,original_restlens)
        ids = [id1,id2]
        original_values = original_restlens[[id1,id2]]
        RB.ManualActuator(actid,ids,original_values,RB.Ganged())
    end
    acs = [ifelse(isodd(i),ganged_act(i,2(i-1)+1,2i,original_restlens),
                           ganged_act(i,2i,2(i-1)+1,original_restlens)) for i = 1:nbodies-1]
    hub = (actuators=acs,)
    pjs = [
        RB.PinJoint(i,RB.Hen2Egg(i,RB.ID(rbs[i],2,1),RB.ID(rbs[i+1],1)))
        for i = 1:nbodies-1
    ]
    jointed = RB.join(pjs,indexed)

    cnt = RB.Connectivity(numbered,indexed,tensioned,jointed)
    st = RB.Structure(rigdibodies,tensiles,cnt)
    RB.Robot(st,hub)
end

function Tbars(;θ = 0)
    SVo3 = SVector{3}([0.0,0.0,0.0])
    a = 1.0
    b = 1.0
    c = 0.05
    function make_base(i)
        movable = true
        constrained = true
        mass_center_position = SVo3
        r̄p1 = SVector{3}([ -a, b,c])
        r̄p2 = SVector{3}([ -a,-b,c])
        r̄p3 = SVector{3}([0.0, b,c])
        r̄p4 = SVector{3}([0.0,-b,c])
        r̄p5 = SVo3
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
        ci = collect(1:12)
        cstr_idx = Int[]    
        R = RotZ(0.0)    
        ω = SVo3
        ri = SVo3
        ro = ri
        ṙo = zero(ro)

        prop = RB.RigidBodyProperty(
                    i,movable,m,
                    Ī,
                    mass_center_position,
                    nodes_positions,
                    axes;
                    constrained
                    )

        nmcs = RB.NCF.NC1P3V(ri, ro, R)

        state = RB.RigidBodyState(prop, nmcs, ri, R, ṙo, ω, ci, cstr_idx)
        basemesh = load(RB.assetpath("装配体1.STL")) |> make_patch(;
            # trans=[-1.0,0,0],
            rot = RotZ(π),
            scale=1/500,
            color = :silver,
        )
        RB.RigidBody(prop,state,basemesh)
    end
    function make_slider(i;ri=SVo3,R=RotZ(0.0),)
        movable = true
        constrained = false
        mass_center_position = SVo3
        r̄p1 = SVector(0.0,0.0,c)
        nodes_positions = [r̄p1,]
        axes = [SVector{3}([1.0,0.0,0.0]),]
        Ī = SMatrix{3,3}(Matrix(1.0I,3,3))
        ω = SVo3
        ro = ri
        ṙo = zero(ro)
        m = 1.0

        prop = RB.RigidBodyProperty(
            i,movable,m,
            Ī,
            mass_center_position,
            nodes_positions,
            axes;
            constrained
        )

        nmcs = RB.NCF.NC1P3V(ri, ro, R)
        # @show q[1:3]
        # @show q[4:6]
        # @show q[7:9]
        # @show q[10:12]
        ro = ri
        ṙo = zero(ro)
        ω = zero(ro)
        state = RB.RigidBodyState(prop, nmcs, ri, R, ṙo, ω)
        slidermesh = load(RB.assetpath("装配体2.2.STL")) |> make_patch(;
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
        RB.RigidBody(prop,state,slidermesh)
    end
    base = make_base(1)
    p1 = SVector(-b*cos(θ),0.0,0.0)
    p2 = SVector(0,b*sin(θ),0)
    slider1 = make_slider(2;ri = p1)
    slider2 = make_slider(3;ri = p2,)
    bar = make_3d_bar(
        4,
        p1,
        p2;
		loadmesh = true,
    )
    rbs = [base,slider1,slider2,bar]
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

    ncables = 4
    original_restlens = zeros(ncables)
    original_restlens = zeros(ncables)
    ks = fill(100.0,ncables)
    ks[2] = 400.0
    cs = zeros(ncables)
    ss = [RB.Cable3D(i, original_restlens[i],ks[i],cs[i];slack=false) for i = 1:ncables]
    tensiles = (cables=ss,)

    cm = [
        1 -1  0  0;
        2 -1  0  0;
        # 3  0 -1  0;
        # 4  0 -1  0;
        6 -1  0  0;
        7 -1  0  0;
    ]

    connected = RB.connect(rigdibodies, cm)
    tensioned = @eponymtuple(connected,)

    j1 = RB.PrismaticJoint(1,RB.Hen2Egg(1,RB.ID(base,5,1),RB.ID(slider1,1,1)))
    j2 = RB.PrismaticJoint(2,RB.Hen2Egg(2,RB.ID(base,5,2),RB.ID(slider2,1,1)))

    j3 = RB.PinJoint(3,RB.Hen2Egg(3,RB.ID(bar,1,1),RB.ID(slider1,1,1)))
    j4 = RB.PinJoint(4,RB.Hen2Egg(4,RB.ID(bar,2,1),RB.ID(slider2,1,1)))


    js = [
        j1,j2,
        j3,j4
    ]

    jointed = RB.join(js,indexed)
    cnt = RB.Connectivity(numbered,indexed,tensioned,jointed)
    st = RB.Structure(rigdibodies,tensiles,cnt)
    RB.Robot(st,)
end

function planar_parallel()
    SVo3 = SVector{3}([0.0,0.0,0.0])
    a = 1.0
    b = 0.5
    θ = 2π/3
    function make_base(i)
        movable = true
        constrained = true
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
                    i,movable,m,
                    Ī,
                    mass_center_position,
                    nodes_positions,
                    axes;
                    constrained
                    )

        nmcs = RB.NCF.NC1P3V(ri, ro, R)

        state = RB.RigidBodyState(prop, nmcs, ri, R, ṙo, ω, ci, cstr_idx)

        RB.RigidBody(prop,state)
    end
    function make_platform(i;ri=SVo3,R=RotZ(0.0),)
        movable = true
        constrained = false
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
                    i,movable,m,
                    Ī,
                    mass_center_position,
                    nodes_positions,
                    axes;
                    constrained
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
    ss = [RB.Cable3D(i, original_restlens[i],ks[i],cs[i];slack=false) for i = 1:ncables]
    tensiles = (cables=ss,)

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

    # j1 = RB.PrismaticJoint(1,RB.Hen2Egg(1,RB.ID(base,5,1),RB.ID(slider1,1,1)))
    # j2 = RB.PrismaticJoint(2,RB.Hen2Egg(2,RB.ID(base,5,2),RB.ID(slider2,1,1)))

    # js = [
    #     j1,j2
    # ]
    # jointed = RB.join(js,indexed)
    cnt = RB.Connectivity(numbered,indexed,tensioned,)
    st = RB.Structure(rigdibodies,tensiles,cnt)
    RB.Robot(st,)
end
