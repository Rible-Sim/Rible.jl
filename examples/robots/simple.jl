function simple(;
        c=0.0,
        m = 4,
        α = 2π/m,
        k = 100.0,
        z0 = 0.2,
        ωz = 5.0,
        mbar = 0.05,
        free = false,
    )
    lₛ = 15e-3
    DE = 300e-3
    d = (150*√2)*1e-3
    r = d/2
    bps = [
        [r*cos(i*α),r*sin(i*α),0.0]
        for i = 1:m
    ] .|> SVector{3}
    ro = SVector(0,0,z0)
    Rplate = RotX(π/4)
    Rbar = RotY(-π/12)
    p = Ref(Rbar) .* [
        [0,0, DE/2],
        [0,0,-DE/2]
    ] .+ Ref(ro) .|> SVector{3}

    # fig,ax,plt = scatter(p)
    # # ax.aspect = DataAspect()
    # fig
    ṙo_ref = SVector(0.0,0,0)

    rbs = [
        make_3d_bar(
            1, 
            p[end-1],
            p[end];
            # ṙi = ṙo_ref, 
            # ṙj = ṙo_ref,
            radius_ratio = 1/120,
            m = mbar,
            mat_name = "Teak",
            barcolor=:darkred,
        ),
        make_3d_plate(
            2, 
            bps,ro,
            Rplate,ro;
            m,
            radius = r,
            # visible = !free,	
            pres_idx = ifelse(free,Int[],collect(1:12)),
            cstr_idx = ifelse(free,collect(1:6),Int[]),
        )
    ]

    rigdibodies = TypeSortedCollection(rbs)
    #
    ncables = 2m

    original_restlens = zeros(ncables)
    original_restlens[  1: m] .= 0.05
    original_restlens[m+1:2m] .= 0.1

    ks = zeros(ncables)
    ks[  1: m] .= ks[m+1:2m] .= k

    spring_dampers = [
        RB.DistanceSpringDamper3D( original_restlens[i], ks[i], c;slack=true) for i = 1:ncables
    ]

    for i = 1:2m
        @info "No.$i, Rest length=$(original_restlens[i]), Stiffness = $(ks[i]), Damping = $c";
    end
 
    # #
    connecting_matrix = zeros(Int,ncables,2)
    for j = 1:ncables
    if j <= m
        connecting_matrix[j, 1:2] = [1, -j]
    else
        connecting_matrix[j, 1:2] = [2, -(j-m)]
    end
    end
    display(connecting_matrix)
    cm = RT.change_connecting_format(rigdibodies, connecting_matrix)
    cables = RT.connect_spring(rigdibodies, spring_dampers; cm )

    apparatuses = cables

    cnt = if free
        RB.Connectivity(rigdibodies,apparatuses)
    else
        RB.PresFreeConnectivity(rigdibodies,apparatuses)
    end
    
    # #
    st = RT.TensegrityStructure(rigdibodies, apparatuses, cnt, )
    
    capta_gauges = Int[]
    error_gauges = Int[]
    actuators = Int[]
    ## actuators = [
    ##     RB.RegisterActuator(
    ##         1,
    ##         collect(1:ncables_prism),
    ##         zeros(ncables_prism)
    ##     ),
    ## ]
    hub = RB.ControlHub(
        st,
        capta_gauges,
        error_gauges,
        actuators,
        RB.Coalition(st,capta_gauges,error_gauges,actuators)
    )

    bot = RB.Robot(st, hub)
end