function quad(c=100.0;
            μ = 0.9,
            e = 0.0
        )
    z = 0.3
    b = 0.4
    θ = π/2
    c = 0.13135
    a1 = 0.09112999999999999
    a2 = 0.16244
    a3 = 0.10983
    d = 0.0046
    p0 = [0,0,z]
    p1_to_4    = [[i*c,k*d,z]        for i in [-1,1] for k = [-1,1]]
    p5_to_12   = [[i*c+j*a1,k*a2,z]  for i in [-1,1] for k = [-1,1] for j = [-1,1]]
    p13_to_p20 = [[i*c,k*a3,j*b/2+z] for i in [-1,1] for k = [-1,1] for j = [-1,1]]

    p = OffsetArray(SVector{3}.(vcat([p0],p1_to_4,p5_to_12,p13_to_p20)), 0:20)
    # fig,ax,plt = scatter(p)
    # # ax.aspect = DataAspect()
    # fig
    ṙo_ref = SVector(1.0,0,0)
    function rigidbase(i,p)

        movable = true
        constrained = false
        ci = Int[]
        Φi = collect(1:6)
        ri = p[0]
        ro = ri
        R = SMatrix{3,3}(Matrix(1.0I,3,3))
        ṙo = ṙo_ref
        ω = [0.1,0,0.1]
        r̄g = SVector{3}(0.0,0.0,0.0)
        # @show p
        r̄ps = [p[i]-p[0] for i = 1:12]
        m = inertia = 1.0
        Ī = SMatrix{3,3}(
            [
                inertia 0       0;
                0       inertia 0;
                0       0       inertia
            ]
        )
        prop = TR.RigidBodyProperty(
            i,
            movable,
            m,
            Ī,
            r̄g,
            r̄ps;
            constrained = constrained,
        )
        lncs, _ = TR.NCF.NC1P3V(ri, ro, R, ṙo, ω)
        state = TR.RigidBodyState(prop, lncs, ri, R, ṙo, ω, ci, Φi)
        trunk_mesh = load("身体.STL")
        rb = TR.RigidBody(prop, state, trunk_mesh)
    end

    function rigidbar(i,p)
        movable = true
        constrained = false
        ci = Int[]
        Φi = [1]

        ri = p[12+2(i-1)+1]
        rj = p[12+2(i-1)+2]
        ro = (ri + rj)/2
        u = rj - ri
        u /= norm(u)
        v,w = TR.NCF.HouseholderOrthogonalization(u)
        R = SMatrix{3,3}(hcat(u,v,w))
        ṙo = ṙo_ref
        ω = zero(ri)
        # if isodd(pos)
        b = norm(rj-ri)
        r̄g  = SVector{3}(   0,0.0,0.0)
        r̄p1 = SVector{3}(-b/2,0.0,0.0)
        r̄p2 = SVector{3}( b/2,0.0,0.0)
        m = inertia = 0.2
        Ī = SMatrix{3,3}([
            inertia 0 0;
            0       0 0;
            0       0 0
        ])
        r̄ps = [r̄p1,r̄p2]
        prop = TR.RigidBodyProperty(
            i,
            movable,
            m,
            Ī,
            r̄g,
            r̄ps;
            constrained = constrained,
        )
        lncs, _ = TR.NCF.NC3D2P(ri, rj, ro, R, ṙo, ω)
        state = TR.RigidBodyState(prop, lncs, ro, R, ṙo, ω, ci, Φi)
        leg_mesh = load("400杆.STL")
        rb = TR.RigidBody(prop, state, leg_mesh)
    end
    rbs = vcat(
        rigidbase(5, p),[
            rigidbar(i, p) for i = 1:4
        ]
    )
    rigdibodies = TypeSortedCollection(rbs)
    numberedpoints = TR.number(rigdibodies)
    indexedcoords = TR.index(rigdibodies)
    #
    ncables = 4*6
    original_restlens = zeros(ncables)
    ks = zeros(ncables)
    for j = 1:ncables
        lev,i = divrem(j-1,6)
        if i < 3
            original_restlens[j] = 0.1
            ks[j] = 200
        else
            original_restlens[j] = 0.1
            ks[j] = 180
        end
    end
    cables = [
        TR.Cable3D(i, original_restlens[i], ks[i], c) for i = 1:ncables
    ]
    #
    tensiles = (cables = cables,)
    acs = [
        TR.ManualActuator(
            i,
            collect(1:6), [original_restlens[6*(i-1)+j] for j = 1:6],
        ) for i = 1:4
    ]
    hub = (actuators = acs,)
    # #
    tripoints = [
        [1,5,6],
        [2,7,8],
        [3,9,10],
        [4,11,12]
    ]

    matrix_cnt = zeros(Int,ncables,5)
    for j = 1:ncables
        lev,i = divrem(j-1,6)
        if i < 3
            matrix_cnt[j, lev+1] =  1
            matrix_cnt[j,     5] = -tripoints[lev+1][i+1]
        else
            matrix_cnt[j, lev+1] =  2
            matrix_cnt[j,     5] = -tripoints[lev+1][i+1-3]
        end
    end
    # display(matrix_cnt)
    connected = TR.connect(rigdibodies, matrix_cnt)
    tensioned = @eponymtuple(connected,)
    #
    cnt = TR.Connectivity(numberedpoints, indexedcoords, tensioned)
    # #

    contacts = [TR.Contact(i,μ,e) for i = 1:4]

    tg = TR.TensegrityStructure(rigdibodies, tensiles, cnt, contacts)
    bot = TR.TensegrityRobot(tg, hub)
end

function rigidbar(i,ri,rj;
        ṙi,ṙj,
        m = 0.05,
        movable = true,
        constrained = false,
        ci = Int[],
        Φi = [1],
        isbody = false,
    )
    ro = (ri + rj)/2
    ṙo = (ṙi + ṙj)/2
    u = rj - ri
    u /= norm(u)
    v,w = Meshes.householderbasis(u)
    R = SMatrix{3,3}(hcat(u,v,w))
    # if isodd(pos)
    b = norm(rj-ri)
    r̄g  = SVector{3}(   0,0.0,0.0)
    r̄p1 = SVector{3}(-b/2,0.0,0.0)
    r̄p2 = SVector{3}( b/2,0.0,0.0)
    r̄p3 = SVector{3}(   0,0.0,0.0)
    inertia = m*b^2/12
    Ī = SMatrix{3,3}([
        inertia 0 0;
        0       0 0;
        0       0 0
    ])
    r̄ps = [r̄p1,r̄p2,r̄p3]
    ās = [SVector(1.0,0,0),SVector(1.0,0,0),SVector(1.0,0,0)]
    prop = TR.RigidBodyProperty(
        i,
        movable,
        m,
        Ī,
        r̄g,
        r̄ps,
        ās;
        constrained = constrained,
    )
    pretty_table(
        SortedDict(
            [
                ("id", i),
                ("bar length", b),
            ]
        )
    )
    if isbody
        @show isbody
        lncs, _ = TR.NCF.NC1P3V(ri, ro, R)
        ω = zero(ṙo)
    else
        lncs, _ = TR.NCF.NC3D1P1V(ri, u, ro, R)
        ω = TR.NCF.find_ω(lncs,vcat(ri,rj),vcat(ṙi,ṙj))
    end
    state = TR.RigidBodyState(prop, lncs, ro, R, ṙo, ω, ci, Φi)
    # leg_mesh = load("400杆.STL")
    barmesh = endpoints2mesh(r̄p1,r̄p2;radius = norm(r̄p2-r̄p1)/40)
    rb = TR.RigidBody(prop, state, barmesh)
end

function uni(c=100.0;
        μ = 0.9,
        e = 0.0,
        z0 = 0.2,
        ωz = 5.0,
        mbar = 0.05,
        isbody = false,
    )
    lₛ = 80e-3
    DE = 130.64e-3
    θ = 2π/3
    p = SVector{3}.(
        [   
            [0,0,0.0],
            lₛ/2 .*[cos(0θ), sin(0θ),0],
            lₛ/2 .*[cos(1θ), sin(1θ),0],
            lₛ/2 .*[cos(2θ), sin(2θ),0],
            [0,0,-DE/2],
            [0,0, DE/2]
        ]
    )
    # fig,ax,plt = scatter(p)
    # # ax.aspect = DataAspect()
    # fig
    ṙo_ref = SVector(0.0,0,0)
    function rigidbase(i,p)
        movable = true
        constrained = true
        ci = collect(1:12)
        Φi = Int[]
        ri = p[1]
        ro = ri
        R = SMatrix{3,3}(Matrix(1.0I,3,3))
        ṙo = ṙo_ref
        ω = [1.0,0,ωz]
        # ω = zeros(3)
        r̄g = SVector{3}(0.0,0.0,0.0)
        # @show p
        r̄ps = p[[2,3,4,1]]
        ās = [SVector(0,0,1.0) for i = 1:4]
        m = inertia = 0.1
        Ī = SMatrix{3,3}(
            [
                inertia 0       0;
                0       inertia 0;
                0       0       inertia
            ]
        )
        prop = TR.RigidBodyProperty(
            i,
            movable,
            m,
            Ī,
            r̄g,
            r̄ps,
            ās;
            constrained = constrained,
        )
        lncs, _ = TR.NCF.NC1P3V(ri, ro, R, ṙo, ω)
        state = TR.RigidBodyState(prop, lncs, ri, R, ṙo, ω, ci, Φi)
        # trunk_mesh = load("身体.STL")
        long = 25

        b1 = endpoints2mesh(r̄ps[1],r̄ps[2];radius = norm(r̄ps[2]-r̄ps[1])/long)
        b2 = endpoints2mesh(r̄ps[2],r̄ps[3];radius = norm(r̄ps[3]-r̄ps[2])/long)
        b3 = endpoints2mesh(r̄ps[3],r̄ps[1];radius = norm(r̄ps[3]-r̄ps[1])/long)
        trimesh = GB.merge([b1,b2,b3])
        TR.RigidBody(prop, state, trimesh)
    end

    rb1 = rigidbar(
        1, 
        RotX(0.0)*p[end-1],
        RotX(0.0)*p[end];
        ṙi = ṙo_ref,
        ṙj = ṙo_ref,
        m = mbar,
        movable = true,
        constrained = false,
        ci = Int[],
        Φi = collect(1:6),
        isbody,
    )
    rb2 = rigidbase(2, p)
    rbs = [
        rb1,
        rb2
    ]
    rigdibodies = TypeSortedCollection(rbs)
    numberedpoints = TR.number(rigdibodies)
    indexedcoords = TR.index(rigdibodies)
    #
    ncables = 6

    original_restlens = zeros(ncables)
    original_restlens[1:3] .= 80e-3
    original_restlens[4:6] .= 80e-3

    ks = zeros(ncables)
    ks[1:3] .= 1000.0
    ks[4:6] .= 1000.0

    cables = [
    TR.Cable3D(i, original_restlens[i], ks[i], c;slack=false) for i = 1:ncables
    ]
    #
    tensiles = (cables = cables,)
    acs = [
    TR.ManualActuator(
        i,
        collect(1:6), [original_restlens[6*(i-1)+j] for j = 1:6],
    ) for i = [1]
    ]
    hub = (actuators = acs,)
    # #
    matrix_cnt = zeros(Int,ncables,2)
    for j = 1:ncables
    if j < 4
        matrix_cnt[j, 1:2] = [1, -j]
    else
        matrix_cnt[j, 1:2] = [2, -(j-3)]
    end
    end
    # display(matrix_cnt)
    connected = TR.connect(rigdibodies, matrix_cnt)
    tensioned = @eponymtuple(connected,)
    
    # uj = TR.UniversalJoint(1,
    #     TR.End2End(
    #         1,
    #         TR.ID(rb2,4,4),
    #         TR.ID(rb1,3,3),
    #     )
    # )

    uj = TR.PrismaticUniversalJoint(1,
        TR.End2End(
            1,
            TR.ID(rb2,4,4),
            TR.ID(rb1,3,3),
        )
    )


    jointed = TR.join([uj],indexedcoords)

    cnt = TR.Connectivity(numberedpoints, indexedcoords, tensioned, jointed)
  
    contacts = [TR.Contact(i,μ,e) for i = 1:2]

    tg = TR.TensegrityStructure(rigdibodies, tensiles, cnt, contacts)
    bot = TR.TensegrityRobot(tg, hub)
end

function superball(c=0.0;
            μ = 0.9,
            e = 0.0,
            θ = atan(0.5,1),
            l = 1.7/2,
            d = l/2,
            z0 = l^2/(sqrt(5)*d) - 1e-7,
            k = 4000.0,
            constrained = false,
            addconst = Float64[],
        )
    p = Ref(RotY(θ)) .* SVector{3}.(
        [
            [ 0,  d,  l], [ 0,  d, -l],
            [ 0, -d, -l], [ 0, -d,  l],
            [ d,  l,  0], [ d, -l,  0],
            [-d, -l,  0], [-d,  l,  0],
            [ l,  0,  d], [-l,  0,  d],
            [-l,  0, -d], [ l,  0, -d],
        ]
    ).+Ref([0,0,z0])
    # p |> display
    ṗ = [
        SVector(5.0,1.0,0) + SVector(0.0,0.0,0.0)×r
        for r in p
    ]
    rbs = [
        rigidbar(i,
            p[2i-1],
            p[2i  ];
            ṙi=ṗ[2i-1],
            ṙj=ṗ[2i  ], 
            m = 5.0,
            constrained = ifelse(i==1,true,false)
            )
        for i = 1:6
    ]
    rigdibodies = TypeSortedCollection(rbs)
    numberedpoints = TR.number(rigdibodies)
    indexedcoords = TR.index(rigdibodies)
    #
    ncables = 24

    original_restlens = zeros(ncables)
    original_restlens .= 0.996
    ks = zeros(ncables)
    ks .= k

    pretty_table(
        SortedDict(
            [
                ("stiffness", k),
                ("bar length", 2l),
                ("bar distance", d)
            ]
        )
    )
    cables = [
        TR.Cable3D(i, original_restlens[i], ks[i], c;slack=false) for i = 1:ncables
    ]
    #
    tensiles = (cables = cables,)
    acs = [
        TR.ManualActuator(
            i,
            collect(1:6), [original_restlens[6*(i-1)+j] for j = 1:6],
        ) for i = [1]
    ]
    hub = (actuators = acs,)
    # #
    nb = CircularArray([
        [1,2],
        [3,4],
        [5,6]
    ])
    matrix_cnt = zeros(Int,ncables,6)
    for lev = 1:3
        matrix_cnt[8(lev-1)+1, [nb[lev][1],nb[lev+1][1]]] = [1, -1]
        matrix_cnt[8(lev-1)+2, [nb[lev][1],nb[lev+1][1]]] = [2, -1]
        matrix_cnt[8(lev-1)+3, [nb[lev][1],nb[lev+1][2]]] = [1, -2]
        matrix_cnt[8(lev-1)+4, [nb[lev][1],nb[lev+1][2]]] = [2, -2]
        matrix_cnt[8(lev-1)+5, [nb[lev][2],nb[lev+1][2]]] = [1, -1]
        matrix_cnt[8(lev-1)+6, [nb[lev][2],nb[lev+1][2]]] = [2, -1]
        matrix_cnt[8(lev-1)+7, [nb[lev][2],nb[lev+1][1]]] = [1, -2]
        matrix_cnt[8(lev-1)+8, [nb[lev][2],nb[lev+1][1]]] = [2, -2]
    end
    # display(matrix_cnt)
    connected = TR.connect(rigdibodies, matrix_cnt)
    tensioned = @eponymtuple(connected,)
    #
    
    if isempty(addconst)
        cnt = TR.Connectivity(numberedpoints, indexedcoords, tensioned)
    else
        cst1 = TR.LinearJoint(
            1,
            1,
            [0.0],
            addconst
        )
        jointed = TR.join([cst1],indexedcoords)
        cnt = TR.Connectivity(
            numberedpoints, 
            indexedcoords, 
            tensioned,
            jointed
            )
    end

    # #

    contacts = [TR.Contact(i,μ,e) for i = 1:12]

    tg = TR.TensegrityStructure(rigdibodies, tensiles, cnt, contacts)
    bot = TR.TensegrityRobot(tg, hub)
end
