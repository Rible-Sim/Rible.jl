function superball(c=0.0;
        z0 = l^2/(sqrt(5)*d) - 1e-7,
        origin_position = SVector(0,0,z0),
        θ = atan(0.5,1),
        R = RotY(θ),
        origin_velocity = SVector(0.0,0.0,0),
        ω = SVector(0.0,0.0,0.0),
        μ = 0.9,
        e = 0.0,
        l = 1.7/2,
        d = l/2,
        k = 4000.0,
        visible = true,
        addconst = Float64[],
        loadmesh = true,
    )
    p = Ref(R) .* SVector{3}.(
    [
        [ 0,  d,  l], [ 0,  d, -l],
        [ 0, -d, -l], [ 0, -d,  l],
        [ d,  l,  0], [ d, -l,  0],
        [-d, -l,  0], [-d,  l,  0],
        [ l,  0,  d], [-l,  0,  d],
        [-l,  0, -d], [ l,  0, -d],
    ]
    ).+Ref(origin_position)
    # p |> display
    ṗ = [
    origin_velocity + ω×(r-origin_position)
    for r in p
    ]
    rbs = [
    rigidbar(i,
        p[2i-1],
        p[2i  ];
        ṙi=ṗ[2i-1],
        ṙj=ṗ[2i  ], 
        m = 5.0,
        μ,
        e,
        visible = ifelse(i==1 && visible,true,false),
        ci = ifelse(i==1 && visible,collect(1:6),Int[]),
        cstr_idx = ifelse(i==1 && visible,Int[],[1]),
        loadmesh,
        isbody =false,
        )
    for i = 1:6
    ]
    rigdibodies = TypeSortedCollection(rbs)
    numberedpoints = RB.number(rigdibodies)
    indexedcoords = RB.index(rigdibodies)
    #
    ncables = 24

    original_restlens = zeros(ncables)
    original_restlens .= 1.199744871391589
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
    RB.Cable3D(i, original_restlens[i], ks[i], c;slack=false) for i = 1:ncables
    ]
    #
    tensiles = (cables = cables,)
    acs = [
    RB.ManualActuator(
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
    connected = RB.connect(rigdibodies, matrix_cnt)
    tensioned = @eponymtuple(connected,)
    #

    if isempty(addconst)
    cnt = RB.Connectivity(numberedpoints, indexedcoords, tensioned)
    else
    cst1 = RB.LinearJoint(
        1,
        1,
        [0.0],
        addconst
    )
    jointed = RB.join([cst1],indexedcoords)
    cnt = RB.Connectivity(
        numberedpoints, 
        indexedcoords, 
        tensioned,
        jointed
        )
    end

    st = RB.Structure(rigdibodies, tensiles, cnt, )
    bot = RB.Robot(st, hub)
end