# require rigidbar.jl
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
        orbital_position = origin_position,
        orbital_velocity = origin_velocity,
        use_orbit_apparatus=false,
        isbody = false,
        visible = true,
        constrained = false,
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
    @debug "p" p
    ṗ = [
        origin_velocity + ω×(r-origin_position)
        for r in p
    ]
    ṗ |> display
    rbs = [
        rigidbar(
            i,
            p[2i-1],
            p[2i  ];
            ṙi=ṗ[2i-1],
            ṙj=ṗ[2i  ], 
            m = 5.0,
            μ,
            e,
            visible,
            loadmesh,
            isbody =false,
        )
        for i = 1:6
    ]
    rigdibodies = TypeSortedCollection(rbs)
    #
    ncables = 24

    original_restlens = zeros(ncables)
    original_restlens .= 1.199744871391589
    original_restlens .= 0.996
    ks = zeros(ncables)
    ks .= k
    @debug "superball" stiffness=k bar_length = 2l bar_distance = d

    spring_dampers = [
        RB.DistanceSpringDamper3D( original_restlens[i], ks[i], c;slack=false) for i = 1:ncables
    ]
    #
    # #
    nb = CircularArray(
        [
            [1,2],
            [3,4],
            [5,6]
        ]
    )
    connecting_matrix = zeros(Int,ncables,6)
    for lev = 1:3
        connecting_matrix[8(lev-1)+1, [nb[lev][1],nb[lev+1][1]]] = [1, -1]
        connecting_matrix[8(lev-1)+2, [nb[lev][1],nb[lev+1][1]]] = [2, -1]
        connecting_matrix[8(lev-1)+3, [nb[lev][1],nb[lev+1][2]]] = [1, -2]
        connecting_matrix[8(lev-1)+4, [nb[lev][1],nb[lev+1][2]]] = [2, -2]
        connecting_matrix[8(lev-1)+5, [nb[lev][2],nb[lev+1][2]]] = [1, -1]
        connecting_matrix[8(lev-1)+6, [nb[lev][2],nb[lev+1][2]]] = [2, -1]
        connecting_matrix[8(lev-1)+7, [nb[lev][2],nb[lev+1][1]]] = [1, -2]
        connecting_matrix[8(lev-1)+8, [nb[lev][2],nb[lev+1][1]]] = [2, -2]
    end
    # display(connecting_matrix)
    cm = RT.change_connecting_format(rigdibodies, connecting_matrix)
    cables = RT.connect_spring(rigdibodies, spring_dampers; cm, istart=0)
    num_of_cables = length(cables)
    apparatuses = TypeSortedCollection(
        cables
    )
    #

    
    if constrained
        cnt = RB.PresFreeConnectivity(
        rigdibodies,
        apparatuses#; sharing_matrix
    )
    else
        cnt = RB.Connectivity(
            rigdibodies,
            apparatuses#; sharing_matrix
        )
    end
    st = RT.TensegrityStructure(rigdibodies, apparatuses, cnt,)
    
    pos_vel_errors = RB.ErrorGauge(
        1,
        RB.Signifier(rbs[1],2),
        RB.PositionCaptum(),
        [5.0,0,0],
        # RB.PosVelCaptum(),
        # [5.0,0,0,0,0,0],
        # RB.VelocityCaptum(),
        # zeros(3)
    )
    com_pos_vel_errors = RB.ErrorGauge(
        1,
        st,
        RB.CoMPosVelCaptum(),
        [5.0,0,0,0,0,0],
    )
    ## angle_err = RB.ErrorGauge(
    ##     RB.Gauge(
    ##         RB.Anchor(j3,3),
    ##         RB.AngularPositionCaptum()
    ##     ),
    ##     π/2
    ## )
    fullstate_captum = RB.CaptumGauge(
        1,
        st,
        RB.FullStateCaptum(),
    )
    capta_gauges = [fullstate_captum]
    error_gauges = [pos_vel_errors]
    actuators = [
        RB.RegisterActuator(
            1,
            [
                (cable=cable,val_id=cable.id) for cable in cables # signifier
            ],
             RB.FunctionOperator(
                (f, values, u) -> begin
                    @inbounds @. f = values + u
                end,
                (jac, values, u) -> nothing,
                copy(original_restlens), # func_vals
                Matrix{eltype(original_restlens)}(I(length(original_restlens))),  # Jac_vals
            ), # operator
            (
                values=original_restlens,
            ) # register
        ),
    ]
    hub = RB.ControlHub(
        st,
        capta_gauges,
        error_gauges,
        actuators,
        RB.Coalition(
            st,capta_gauges,error_gauges,actuators,
        ),
        # control,
    )
    bot = RB.Robot(st,hub)
end
