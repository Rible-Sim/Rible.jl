function one_bar_one_tri()
    n = 2
    b1 = 0.05*√2
    bps = SVector{2}.(
        [
           [  0.0,   0.0],
           [b1/√2,-b1/√2]
        ]
    )

    α1 = -π/4
    rb1 = build_2d_bar(1,bps[1],bps[2];α=α1,ci=collect(1:2))
    α2 =  0.0
    rb2 = build_2d_tri(2,bps[2];α=α2)
    rbs = TypeSortedCollection([rb1,rb2])
    numberedpoints = RB.number(rbs)
    matrix_sharing = [
        3 1;
        4 2;
    ]
    indexedcoords = RB.index(rbs,matrix_sharing)
    # cables
    cables = Int[]
    apparatuses = (cables = cables,)
    hub = nothing
    connected = RB.connect(rbs,zeros(Int,0,0))

    cnt = RB.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,))

    st = RB.Structure(rbs,apparatuses,cnt)
    bot = RB.Robot(st,hub)
end