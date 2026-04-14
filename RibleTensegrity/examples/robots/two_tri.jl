
function two_tri(;k=100.0,c=0.0,ratio=0.8)
    
    ro = SVector(0.0,0.0)
    α_tri1 = -π
    α_tri2 =  0
    rb1 = build_2d_tri(1,ro;
        α=α_tri1,
        pres_idx=collect(1:6),
        cstr_idx=Int[]
    )
    rb2 = build_2d_tri(2,ro;
        α=α_tri2,
        pres_idx=collect(1:2)
    )

    rbs = TypeSortedCollection((rb1,rb2))
    #
    ncables = 4
    restlen1 = 0.05
    restlens = fill(restlen1,ncables)
    ks = fill(k,ncables)
    cs = fill(c,ncables)
    spring_dampers = [RB.DistanceSpringDamper2D(restlens[i],ks[i],cs[i];slack=true) for i = 1:ncables]

    connecting_matrix = [
        3 -2 ;
       -2  3 ;
        5 -2 ;
       -4  3 ;
        # 5 -2 ;
        # 4 -3 ;
    ]
    cm = RT.change_connecting_format(rbs, connecting_matrix)
    cables = RibleTensegrity.connect_spring(rbs,spring_dampers;cm,)
    #
    # cst1 = RB.PinJoint(
    #     1,
    #     RB.Hen2Egg(
    #         1,
    #         RB.Anchor(rb1,1),
    #         RB.Anchor(rb2,1)
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

    # apparatuses = TypeSortedCollection((cst1,),indexed)

    apparatuses = TypeSortedCollection(
        cables
    )

    
    sharing_matrix = [
        1 1;
        2 2;
    ]
    cnt = RB.PresFreeConnectivity(rbs,apparatuses;sharing_matrix)

    st = RT.TensegrityStructure(rbs,apparatuses,cnt)
    ## acs = [RB.RegisterActuator(1,1:ncables,restlens,RB.Uncoupled())]
    ## hub = (actuators = acs,)
    capta_gauges = Int[]
    error_gauges = Int[]
    actuators = Int[]
    hub = RB.ControlHub(
        st,
        capta_gauges,
        error_gauges,
        actuators,
        RB.Coalition(st,capta_gauges,error_gauges,actuators)
    )
    bot = RB.Robot(st,hub)
end