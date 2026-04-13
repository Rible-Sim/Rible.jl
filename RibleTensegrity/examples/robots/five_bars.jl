# includet("../bodies/make_3d_bar.jl")
function find_R(ri,rj,)
    u = rj - ri
    bar_length = norm(u)	
    û = u./bar_length
    v̂ = [0,1.0,0.0]
    ŵ = cross(û,v̂)
    R = SMatrix{3,3}(hcat(û,v̂,ŵ))
end

function inertia_tensor(mass,bar_length)
    bar_inertia = 1/12*mass*bar_length^2
    Īg = SMatrix{3,3}(
        [	
        2e-6bar_inertia 0.0 0.0;
        0.0 bar_inertia + 1e-6bar_inertia 0.0;
        0.0 0.0 bar_inertia + 1e-6bar_inertia

        ]
    )
end
function five_bars(;
        restitution_coefficient = 0.5,
        friction_coefficient = 0.5,
        coordsType = :NCF,
        err_ref = [0.5, 0.0, -2.0],
        L0s = [
            sqrt(2^2 + 1^2),  # m, rest length of spring 1 
            sqrt(2^2 + 0.5^2),  # m, rest length of spring 2
        ]
    )
    # Define parameters
    m1 = 1.0  # kg
    m2 = 1.5  # kg
    m3 = 1.0  # kg
    m4 = 1.5  # kg

    k1 = 100.0  # N/m
    k2 = 100.0  # N/m

    bar1p1, bar1p2 = [ 0.0, 0.0,  0.0], [-1.0, 0.0,  -1.0]
    bar2p1, bar2p2 = [-1.0, 0.0, -1.0], [ 0.5, 0.0,  -2.0]
    bar3p1, bar3p2 = [ 1.0, 0.0,  0.0], [ 2.0, 0.0, -1.0]
    bar4p1, bar4p2 = [ 2.0, 0.0, -1.0], [0.5, 0.0, -2.0]
    l1 = norm(bar1p1 .- bar1p2)
    l2 = norm(bar2p1 .- bar2p2)
    l3 = norm(bar3p1 .- bar3p2)
    l4 = norm(bar4p1 .- bar4p2)
    Īg1 = inertia_tensor(m1,l1)
    Īg2 = inertia_tensor(m2,l2)
    Īg3 = inertia_tensor(m3,l3)
    Īg4 = inertia_tensor(m4,l4)
    # Define the bars based on the diagram
    bar1 = make_3d_bar(
        1,
        bar1p1, bar1p2;
        m = m1, loadmesh = false,
        coordsType = RB.NCF.NC2P2V,
        R = find_R(bar1p1, bar1p2),
        Īg_override = Īg1
    )  # Bar between joint A and 1
    bar2 = make_3d_bar(
        2,
        bar2p1, bar2p2;
        m = m2, loadmesh = false, 
        coordsType = RB.NCF.NC2P2V,
        R = find_R(bar2p1, bar2p2),
        Īg_override = Īg2   
    )  # Bar between joint 1 and 2
    bar3 = make_3d_bar(
        3, 
        bar3p1, bar3p2;
        m = m3, loadmesh = false,
        coordsType = RB.NCF.NC2P2V,
        R = find_R(bar3p1, bar3p2),
        Īg_override = Īg3
    )  # Bar between joint 2 and 3
    bar4 = make_3d_bar(
        4, 
        bar4p1, bar4p2;
        m = m4, loadmesh = false,
        coordsType = RB.NCF.NC2P2V,
        R = find_R(bar4p1, bar4p2),
        loci_mode = :mass,
        Īg_override = Īg4
    )  # Bar between joint 3 and B

    rbs = [bar1,bar2,bar3,bar4]
    rigdibodies = TypeSortedCollection(rbs)

    ncables = 2
    original_restlens = L0s
    ks = [k1,k2]
    cs = zeros(ncables)
    spring_dampers = [RB.DistanceSpringDamper3D(original_restlens[i],ks[i],cs[i];slack=false) for i = 1:ncables]

    connecting_matrix = [
        2 0 -1 0;
        0 2 -1 0;
    ]
    cm = RT.change_connecting_format(rigdibodies, connecting_matrix)
    cables = RT.connect_spring(rigdibodies, spring_dampers; cm,)

    
    j1 = RB.FixedIndicesApparatus(
        3,
        bar1,
        [1,2,3,7,9],
        [0.0,0.0,0.0,0.0,0.0]
    )
    j2 = RB.FixedIndicesApparatus(
        4,
        bar3,
        [1,2,3,7,9],
        [ 1.0, 0.0,  0.0,0.0,0.0], 
    )

    apparatuses = TypeSortedCollection(
        vcat(cables,[j1, j2])
    )
    
    sharing_matrix = [
        4 1 0 0;
        5 2 0 0;
        6 3 0 0;
        7 7 0 0;
        9 9 0 0;
        0 0 4 1;
        0 0 5 2;
        0 0 6 3;
        0 0 7 7;
        0 0 9 9;
        0 4 0 4;
        # 0 5 0 5;
        0 6 0 6;
    ]
    cnt = RB.Connectivity(rigdibodies, apparatuses; sharing_matrix)
    st = RT.TensegrityStructure(rigdibodies,apparatuses, cnt)

    pos_vel_errors = RB.ErrorGauge(
        1,
        RB.Signifier(bar2,2),
        RB.PositionCaptum(),
        err_ref,
        # RB.PosVelCaptum(),
        # zeros(6),
        # RB.VelocityCaptum(),
        # zeros(3)
    )

    capta_gauges = Int[]
    error_gauges = TypeSortedCollection([pos_vel_errors,])
    actuators = [
        RB.RegisterActuator(
            i,
            [
                (cable=cables[j],val_id=1)
            ], # Device collection,
            # Operator, one control input
            RB.FunctionOperator(
                (f, values, u) -> begin
                    @inbounds @. f = values + u
                end,
                (jac, values, u) -> nothing,
                copy([L0s[j]]), # func_vals
                Matrix{eltype([L0s[j]])}(I(length([L0s[j]]))),  # Jac_vals
            ), # operator
            # Control input values, Reg
            (values = [L0s[j]],)
        )
        for (i,j) in enumerate([1,2])
    ]

    hub = RB.ControlHub(
        st,
        capta_gauges,
        error_gauges,
        actuators,
        RB.Coalition(st,capta_gauges,error_gauges,actuators)
    )

    bot = RB.Robot(st,hub)
    
    bot
end
