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

    contactable = true
    visible = true
    ci = Int[]
    cstr_idx = collect(1:6)
    ri = p[0]
    ro = ri
    R = SMatrix{3,3}(Matrix(1.0I,3,3))
    ṙo = ṙo_ref
    ω = [0.1,0,0.1]
    mass_locus = SVector{3}(0.0,0.0,0.0)
    # @show p
    loci = [p[i]-p[0] for i = 1:12]
    m = inertia = 1.0
    Ī = SMatrix{3,3}(
        [
            inertia 0       0;
            0       inertia 0;
            0       0       inertia
        ]
    )
    prop = RB.RigidBodyProperty(
        i,
        contactable,
        m,
        Ī,
        mass_locus,
        loci;
        visible = visible,
    )
    nmcs = RB.NCF.NC1P3V(ri, ro, R)
    state = RB.RigidBodyState(prop, nmcs, ri, R, ṙo, ω, ci, cstr_idx)
    trunk_mesh = load("身体.STL")
    body = RB.RigidBody(prop, state, trunk_mesh)
    end

    function rigidbar(i,p)
    contactable = true
    visible = true
    ci = Int[]
    cstr_idx = [1]

    ri = p[12+2(i-1)+1]
    rj = p[12+2(i-1)+2]
    ro = (ri + rj)/2
    u = rj - ri
    u /= norm(u)
    v,w = RB.NCF.HouseholderOrthogonalization(u)
    R = SMatrix{3,3}(hcat(u,v,w))
    ṙo = ṙo_ref
    ω = zero(ri)
    # if isodd(pos)
    b = norm(rj-ri)
    mass_locus  = SVector{3}(   0,0.0,0.0)
    r̄p1 = SVector{3}(-b/2,0.0,0.0)
    r̄p2 = SVector{3}( b/2,0.0,0.0)
    m = inertia = 0.2
    Ī = SMatrix{3,3}([
        inertia 0 0;
        0       0 0;
        0       0 0
    ])
    loci = [r̄p1,r̄p2]
    prop = RB.RigidBodyProperty(
        i,
        contactable,
        m,
        Ī,
        mass_locus,
        loci;
        visible = visible,
    )
    nmcs = RB.NCF.NC3D2P(ri, rj, ro, R)
    state = RB.RigidBodyState(prop, nmcs, ro, R, ṙo, ω, ci, cstr_idx)
    leg_mesh = load("400杆.STL")
    body = RB.RigidBody(prop, state, leg_mesh)
    end
    rbs = vcat(
    rigidbase(5, p),[
        rigidbar(i, p) for i = 1:4
    ]
    )
    rigdibodies = TypeSortedCollection(rbs)
    numberedpoints = RB.number(rigdibodies)
    indexedcoords = RB.index(rigdibodies)
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
    RB.DistanceSpringDamper3D(i, original_restlens[i], ks[i], c) for i = 1:ncables
    ]
    #
    tensiles = (cables = cables,)
    acs = [
    RB.ManualActuator(
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
    connected = RB.connect(rigdibodies, matrix_cnt)
    tensioned = @eponymtuple(connected,)
    #
    cnt = RB.Connectivity(numberedpoints, indexedcoords, tensioned)
    # #
    st = RB.Structure(rigdibodies, tensiles, cnt,)
    bot = RB.Robot(st, hub)
end