
function two_tri(;k=100.0,c=0.0,ratio=0.8)
    
    ro = SVector(0.0,0.0)
    α_tri1 = -π
    α_tri2 =  0
    rb1 = build_2d_tri(1,ro;
        α=α_tri1,
        ci=collect(1:6),
        cstr_idx=Int[]
    )
    rb2 = build_2d_tri(2,ro;
        α=α_tri2,
        ci=collect(1:2)
    )

    rbs = TypeSortedCollection((rb1,rb2))
    numbered = RB.number(rbs)
    matrix_sharing = [
        1 1;
        2 2;
    ]
    indexed = RB.index(rbs,matrix_sharing)
    #
    ncables = 4
    restlen1 = 0.05
    restlens = fill(restlen1,ncables)
    ks = fill(k,ncables)
    cs = fill(c,ncables)
    cables = [RB.DistanceSpringDamper2D(restlens[i],ks[i],cs[i];slack=true) for i = 1:ncables]
    acs = [RB.RegisterActuator(1,1:ncables,restlens,RB.Uncoupled())]
    apparatuses = (cables = cables,)
    hub = (actuators = acs,)
    cnt_matrix_cables = [
        3 -2 ;
       -2  3 ;
        5 -2 ;
       -4  3 ;
        # 5 -2 ;
        # 4 -3 ;
    ]
    connected = RB.connect(rbs,cnt_matrix_cables)
    #
    # cst1 = RB.PinJoint(
    #     1,
    #     RB.Hen2Egg(
    #         1,
    #         RB.Signifier(rb1,1),
    #         RB.Signifier(rb2,1)
    #     )
    # )

    # cst2 = RB.LinearJoint(
    #     2,
    #     2,
    #     zeros(2),
    #     begin
    #         A = zeros(2,12)
    #         A[1:2,1:2] = Matrix(1I,2,2)
    #         A
    #     end
    # )

    # jointedmembers = RB.join((cst1,),indexed)

    cnt = RB.Connectivity(
        numbered,
        indexed,
        @eponymtuple(connected,),
        # jointedmembers
    )

    st = RB.Structure(rbs,apparatuses,cnt)
    bot = RB.Robot(st,hub)
end