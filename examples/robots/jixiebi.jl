function build_jixiebi(n; ks1=80.0, restlens1=50.0,
    ks2=1.0, restlens2=100.0,
    kc=1.2, restlenc=160.0, pre=383.65,
    delta=-20)

    function RigidData(r̄g, m, i, r̄l, mesh=nothing)
        @eponymtuple(r̄g, m, i, r̄l, mesh)
    end

    nb = 4n + 1

    len1 = 234.0
    len11 = 45.33
    len12 = 195.33
    r1 = 100.0
    r2 = 120
    len2 = 160.0
    len21 = -112.8
    len22 = 47.2
    r3 = 35.0
    len3 = 416.76
    len31 = -233.56
    len32 = -73.56
    len33 = 33.21
    len34 = 183.21
    len4 = 199.5
    len41 = -163.99
    len42 = -3.99
    len43 = 18.77

    r̄g1 = [0.0, 0.0, 0.0]
    r̄g2 = [0.0, 0.0, 0.0]
    r̄g3 = [0.0, 0.0, 0.0]
    r̄g4 = [0.0, 0.0, 0.0]
    m1 = 0.85
    m2 = 0.05
    m3 = 1.0
    m4 = 0.15

    i1 = [4575.03 0 0; 0 5114.84 0; 0 0 4575.03]
    i2 = [372.01 0 0; 0 358.79 0; 0 0 372.01]
    i3 = [6225.41 0 0; 0 6223.91 0; 0 0 6225.41]
    i4 = [815.03 0 0; 0 108.664 0; 0 0 815.03]
    r̄l1 = [[len11, 0.0, -r1],
        [len11, r1 * sqrt(3) / 2, r1 / 2],
        [len11, -r1 * sqrt(3) / 2, r1 / 2],
        [len12, 0.0, -r1],
        [len12, r1 * sqrt(3) / 2, r1 / 2],
        [len12, -r1 * sqrt(3) / 2, r1 / 2],
        [len11, 0.0, r2],
        [len11, r2 * sqrt(3) / 2, -r2 / 2],
        [len11, -r2 * sqrt(3) / 2, -r2 / 2],
        [len12, 0.0, r2],
        [len12, r2 * sqrt(3) / 2, -r2 / 2],
        [len12, -r2 * sqrt(3) / 2, -r2 / 2]]
    r̄l2 = [[len21, 0.0, -r3],
        [len21, r3 * sqrt(3) / 2, r3 / 2],
        [len21, -r3 * sqrt(3) / 2, r3 / 2],
        [len22, 0.0, -r1],
        [len22, r1 * sqrt(3) / 2, r1 / 2],
        [len22, -r1 * sqrt(3) / 2, r1 / 2],
        [len22, 0.0, r2],
        [len22, r2 * sqrt(3) / 2, -r2 / 2],
        [len22, -r2 * sqrt(3) / 2, -r2 / 2]]
    r̄l3 = [[len31, 0.0, -r3],
        [len31, r3 * sqrt(3) / 2, r3 / 2],
        [len31, -r3 * sqrt(3) / 2, r3 / 2],
        [len32, 0.0, -r1],
        [len32, r1 * sqrt(3) / 2, r1 / 2],
        [len32, -r1 * sqrt(3) / 2, r1 / 2],
        [len33, 0.0, -r1],
        [len33, r1 * sqrt(3) / 2, r1 / 2],
        [len33, -r1 * sqrt(3) / 2, r1 / 2],
        [len34, 0.0, -r1],
        [len34, r1 * sqrt(3) / 2, r1 / 2],
        [len34, -r1 * sqrt(3) / 2, r1 / 2],
        [len32, 0.0, r2],
        [len32, r2 * sqrt(3) / 2, -r2 / 2],
        [len32, -r2 * sqrt(3) / 2, -r2 / 2],
        [len33, 0.0, r2],
        [len33, r2 * sqrt(3) / 2, -r2 / 2],
        [len33, -r2 * sqrt(3) / 2, -r2 / 2],
        [len34, 0.0, r2],
        [len34, r2 * sqrt(3) / 2, -r2 / 2],
        [len34, -r2 * sqrt(3) / 2, -r2 / 2]]
    r̄l4 = [[len41, 0.0, -r3],
        [len41, r3 * sqrt(3) / 2, r3 / 2],
        [len41, -r3 * sqrt(3) / 2, r3 / 2],
        [len42, 0.0, -r1],
        [len42, r1 * sqrt(3) / 2, r1 / 2],
        [len42, -r1 * sqrt(3) / 2, r1 / 2],
        [len42, 0.0, r2],
        [len42, r2 * sqrt(3) / 2, -r2 / 2],
        [len42, -r2 * sqrt(3) / 2, -r2 / 2]]

    mesh1 = load(RB.assetpath("jixiebi/m1.STL")) |> RB.make_patch(; rot=RotZ(-π / 2)) |> RB.make_patch(; rot=RotX(pi / 6))
    mesh2 = load(RB.assetpath("jixiebi/m2.STL")) |> RB.make_patch(; rot=RotZ(-π / 2)) |> RB.make_patch(; rot=RotX(pi / 6))
    mesh3 = load(RB.assetpath("jixiebi/m3.STL")) |> RB.make_patch(; rot=RotZ(-π / 2)) |> RB.make_patch(; rot=RotX(pi / 6))
    mesh4 = load(RB.assetpath("jixiebi/m4.STL")) |> RB.make_patch(; rot=RotZ(-π / 2)) |> RB.make_patch(; rot=RotX(pi / 6))

    rd1 = RigidData(r̄g1, m1, i1, r̄l1, mesh1)
    rd2 = RigidData(r̄g2, m2, i2, r̄l2, mesh2)
    rd3 = RigidData(r̄g3, m3, i3, r̄l3, mesh3)
    rd4 = RigidData(r̄g4, m4, i4, r̄l4, mesh4)
    rds = reduce(vcat, [[rd1], [rd2 for _ in 1:3], reduce(vcat, [vcat([rd3], [rd2 for _ in 1:3]) for _ in 1:n-1]), [rd4]])

    function get_r̄s(rds)
        r̄s = [[0.0, 0.0, 0.0] for _ in 1:nb]
        r̄s[1][1] = -len11
        for i in 1:nb-1
            if length(rds[i].r̄l) == 12
                r̄s[i+1][1] += r̄s[i][1] + len12 - len21 + delta
            elseif length(rds[i].r̄l) == 9
                if length(rds[i+1].r̄l) == 21
                    r̄s[i+1][1] += r̄s[i][1] + len22 - len31 + delta
                elseif i == nb - 1
                    r̄s[i+1][1] += r̄s[i][1] + len22 - len41 + delta
                else
                    r̄s[i+1][1] += r̄s[i][1] + len2 + delta
                end
            elseif length(rds[i].r̄l) == 21
                r̄s[i+1][1] += r̄s[i][1] + len34 - len21 + delta
            end
        end
        return r̄s
    end
    r̄s = get_r̄s(rds)

    function rigidbody(i, rdsi, r̄si)
        if i == 1
            contactable = false
            visible = true
            ci = Int[]
            Φi = collect(1:6)
        else
            contactable = true
            visible = true
            ci = Int[]
            Φi = collect(1:6)
        end

        ri = r̄si
        aps = rdsi.r̄l
        theta = 0
        α = [cos(theta) -sin(theta) 0;
            sin(theta) cos(theta) 0;
            0 0 1]
        r̄g = rdsi.r̄g
        m = rdsi.m
        I = rdsi.i
        nrp = length(aps)
        ṙo = zeros(3)
        ω = zeros(Float64, 3)
        prop = RB.RigidBodyProperty(
            i,
            contactable,
            m,
            SMatrix{3,3}(I),
            SVector{3}(r̄g),
            [SVector{3}(aps[i]) for i in 1:nrp];
            visible
        )
        ro = SVector{3}(ri)

        state = RB.RigidBodyState(prop, ri, α, ṙo, ω)
        nmcs = RB.NCF.NC1P3V(SVector{3}(ri))
        coords = RB.NonminimalCoordinates(nmcs, ci, Φi)
        rb = RB.RigidBody(prop, state, coords, rdsi.mesh)
    end
    rbs = [rigidbody(i, rds[i], r̄s[i]) for i in 1:nb]
    bodies = TypeSortedCollection(rbs)

    nstrings = 24n
    ks1 = ks1
    restlens1 = restlens1
    ks2 = ks2
    restlens2 = restlens2
    spring_dampers = [(i <= 12n) ? RB.DistanceSpringDamper3D(restlens1, ks1, 0.02) : RB.DistanceSpringDamper3D(restlens2, ks2, 0.02) for i in 1:nstrings]

    restlencs = zeros(Float64, 5n-1)
    restlenc1 = 266.77
    for i in 1:5n-1
        if (mod(i,5) == 0) 
            restlencs[i] = restlenc1
        else
            restlencs[i] = restlenc
        end
    end
    ncluster = 3
    cluster_segs = [[RB.DistanceSpringDamperSegment(restlencs[i], kc, prestress=pre) for i in 1:5n-1] for _ in 1:ncluster]
    cluster_sps = [[RB.SlidingPoint(0.01) for i in 1:5n-2] for _ in 1:ncluster]

    matrix_cnt_raw = Vector{Matrix{Int}}()

    for i in 1:n
        for j in 4i-3:4i
            s = zeros(3, 4)
            if (mod(j - 1, 4) == 0) & (j != 1)
                s[1, :] = [j, 10, j + 1, 1]
                s[2, :] = [j, 11, j + 1, 2]
                s[3, :] = [j, 12, j + 1, 3]
            else
                s[1, :] = [j, 4, j + 1, 1]
                s[2, :] = [j, 5, j + 1, 2]
                s[3, :] = [j, 6, j + 1, 3]
            end
            push!(matrix_cnt_raw, s)
        end
    end
    for i in 1:n
        for j in 4i-3:4i
            s = zeros(3, 4)
            if (mod(j - 1, 4) == 0) & (j != 1)
                s[1, :] = [j, 10, j + 1, 4]
                s[2, :] = [j, 11, j + 1, 5]
                s[3, :] = [j, 12, j + 1, 6]
            else
                s[1, :] = [j, 4, j + 1, 4]
                s[2, :] = [j, 5, j + 1, 5]
                s[3, :] = [j, 6, j + 1, 6]
            end
            push!(matrix_cnt_raw, s)
        end
    end
    connecting_matrix = reduce(vcat, matrix_cnt_raw)

    matrix_cnt_raw = [Vector{Matrix{Int}}() for _ in 1:ncluster]

    for i in 1:n
        j = 4i - 3
        s= i==1 ? zeros(4, 4, ncluster) : zeros(5, 4, ncluster)
        if (i == 1)
            for k in 0:ncluster-1
                s[1, :, k+1] = [j, 10+k, j + 1, 7+k]
                s[2, :, k+1] = [j + 1, 7+k, j + 2, 7+k]
                s[3, :, k+1] = [j + 2, 7+k, j + 3, 7+k]
                s[4, :, k+1] = [j + 3, 7+k, j + 4, 13+k]
            end
        elseif (i == n)
            for k in 0:ncluster-1
                s[1, :, k+1] = [j, 13+k, j, 19+k]
                s[2, :, k+1] = [j, 19+k, j + 1, 7+k]
                s[3, :, k+1] = [j + 1, 7+k, j + 2, 7+k]
                s[4, :, k+1] = [j + 2, 7+k, j + 3, 7+k]
                s[5, :, k+1] = [j + 3, 7+k, j + 4, 7+k]
            end
        else
            for k in 0:ncluster-1
                s[1, :, k+1] = [j, 13+k, j, 19+k]
                s[2, :, k+1] = [j, 19+k, j + 1, 7+k]
                s[3, :, k+1] = [j + 1, 7+k, j + 2, 7+k]
                s[4, :, k+1] = [j + 2, 7+k, j + 3, 7+k]
                s[5, :, k+1] = [j + 3, 7+k, j + 4, 13+k]
            end
        end
        for k in 1:ncluster
            push!(matrix_cnt_raw[k], s[:,:,k])
        end
    end
    connecting_cluster_matrix = [reduce(vcat, matrix_cnt_raw[k]) for k in 1:ncluster]
    cables, clusters = RB.connect_spring_and_clusters(bodies, spring_dampers, cluster_sps, cluster_segs, connecting_matrix, connecting_cluster_matrix, istart=ncluster+1)
    cst1 = RB.FixedBodyConstraint(ncluster+1,rbs[1])
    apparatuses = TypeSortedCollection(vcat([cst1],cables, clusters))

    indexed = RB.index(bodies, apparatuses;)
    numbered = RB.number(bodies, apparatuses)
    cnt = RB.Connectivity(numbered, indexed)
    structure = RB.Structure(bodies, apparatuses, cnt)
    # 测量器/传感器， 可测量误差， 之后用于定义cost函数， 实现优化或反馈控制， 目前暂时用不到。
    gauges = Int[]
    # 作动器， 依附于一个或多个bodies或apparatuses， 其类型分派给execute!来施加不同的驱动力
    actuators = [
        # 试下ExternalForceActuator堪不堪用， 不堪用的话再搞个新的Actuator
        RB.ExternalForceActuator(
            i,
            cluster,
            #操作员， 1个驱动量
            RB.NaiveOperator(1),
            #外力， 可能不需要用到
            0.0,
            #驱动量数值，
            [0.0],
        )
        for (i,cluster) in enumerate(clusters)
    ]
    # INFO 3 这个驱动量数值一直没变化，而且没用到。
    hub = RB.ControlHub(
        structure,
        gauges,
        actuators,
        RB.Coalition(structure, gauges, actuators,)
    )
    bot = RB.Robot(structure, hub)
end

