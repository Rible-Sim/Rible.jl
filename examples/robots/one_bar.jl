function one_bar(k=0.0,c=0.0;ratio=0.8)
    n = 1

    bps = SVector{2}.(
        [
        [0.0, 0.0],
        [0.1, 0.0]
        ]
    )

    α1 = 0.0
    rb1 = build_2d_bar(1,bps[1],bps[2],α1)
    rbs = TypeSortedCollection((rb1,))
    numberedpoints = RB.number(rbs)
    matrix_sharing = zeros(Int,0,0)
    indexedcoords = RB.index(rbs,matrix_sharing)
    #
    ss = Int[]
    tensiles = (cables = ss,)
    hub = nothing
    #
    connections = RB.connect(rbs,zeros(Int,0,0))

    jointedmembers = RB.unjoin()

    cnt = RB.Connectivity(numberedpoints,indexedcoords,connections,jointedmembers)

    st = RB.Structure(rbs,tensiles,cnt)
    bot = RB.Robot(st,hub)
end