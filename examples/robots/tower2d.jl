function tower2d(;k=100.0,c=0.0,ratio=0.8,ratio1=ratio,slack=true)
    n = 8

    bps_raw = [zeros(2) for i = 1:11]
    for i = 1:4
        bps_raw[2i-1] .= [0.0,0.1*(i-1)]
        bps_raw[2i  ] .= [0.1,0.1*(i-1)]
    end
    for i = 1:3
        bps_raw[8+i] .= [0.05,0.15+0.1*(i-1)]
    end
    bps = SVector{2}.(bps_raw)
    # display(bps)
    α_bar = π/2
    α_tri = 0.0
    rb1 = build_2d_bar(1,bps[1],bps[3];α=α_bar,ci=collect(1:2))
    rb2 = build_2d_bar(2,bps[2],bps[4];α=α_bar,ci=collect(1:2))
    # rb1 = build_2d_tri(1,bps[1],bps[3];α=α_bar,ci=collect(1:2))
    # rb2 = build_2d_tri(2,bps[2],bps[4];α=α_bar,ci=collect(1:2))
    # rb3 = build_2d_tri(3,bps[3],bps[4],bps[9];α=α_tri)
    rb3 = build_2d_tri(3,bps[3],;α=α_tri)
    rb4 = build_2d_bar(4,bps[9],bps[10];α=α_bar)
    rb5 = build_2d_tri(5,bps[5],bps[6];α=α_tri)
    rb6 = build_2d_bar(6,bps[10],bps[11];α=α_bar)
    rb7 = build_2d_tri(7,bps[7],bps[8],bps[11];α=α_tri)
    rb8 = build_2d_ground(8)
    rbs = TypeSortedCollection((rb1,rb2,rb3,rb4,rb5,rb6,rb7,rb8))
    # rbs = TypeSortedCollection((rb1,rb2,rb3,rb4,rb5,rb6,rb7))
    numberedpoints = RB.number(rbs)
    matrix_sharing = [
        3 0 1 0 0 0 0;
        4 0 2 0 0 0 0;
        # 0 3 3 0 0 0 0;
        # 0 4 4 0 0 0 0;
        0 0 0 0 0 3 5;
        0 0 0 0 0 4 6;
    ]
    indexedcoords = RB.index(rbs,matrix_sharing)
    #
    restlen4 = ratio1*0.1*√5
    restlen1 = ratio*0.1*√2
    restlen2 = ratio*0.1
    restlen3 = ratio*0.05*√2
    restlens = [
        restlen4,restlen4,
        restlen1,restlen1,
        restlen2,restlen2,
        restlen3,restlen3,
        restlen2,restlen2,
        restlen3,restlen3,
    ]
    ncables = length(restlens)
    naux = 2
    ness = 10
    ks = vcat(fill(k,naux),fill(k,ness))
    cs = fill(c,ncables)
    cables = [RB.DistanceSpringDamper2D(i,restlens[i],ks[i],cs[i];slack) for i = 1:ncables]
    acs = [RB.ManualActuator(1,collect(1:ncables),restlens[1:ncables])]
    tensiles = (cables = cables,)
    hub = (actuators = acs,)
    cnt_matrix_cables = [
        0  0 0 0 -1  0  0  2;
        0  0 0 0 -2  0  0  3;
        1 -2 0 0  0  0  0  0;
       -2  1 0 0  0  0  0  0;
        0  0 1 0 -1  0  0  0;
        0  0 2 0 -2  0  0  0;
        0  0 3 0 -1  0  0  0;
        0  0 3 0 -2  0  0  0;
         0  0 0 0  1  0 -1  0;
         0  0 0 0  2  0 -2  0;
         0  0 0 0  3  0 -1  0;
         0  0 0 0  3  0 -2  0;
    ]
    connected = RB.connect(rbs,cnt_matrix_cables)
    #

    cst1 = RB.PinJoint(RB.Hen2Egg(1,RB.ID(rb2,2),RB.ID(rb3,2)))
    cst2 = RB.PinJoint(RB.Hen2Egg(2,RB.ID(rb3,3),RB.ID(rb4,1)))
    cst3 = RB.PinJoint(RB.Hen2Egg(3,RB.ID(rb4,2),RB.ID(rb5,3)))
    cst4 = RB.PinJoint(RB.Hen2Egg(4,RB.ID(rb5,3),RB.ID(rb6,1)))
    jointedmembers = RB.join((cst1,cst2,cst3,cst4),indexedcoords)
    # jointedmembers = RB.join((cst1,cst2,cst3),indexedcoords)

    cnt = RB.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,),jointedmembers)
    # cnt = RB.Connectivity(numberedpoints,indexedcoords,connections)
    st = RB.Structure(rbs,tensiles,cnt)
    bot = RB.Robot(st,hub)
end