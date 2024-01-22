function RigidData(r̄g,m,i,r̄l,mesh=nothing)
    @eponymtuple(r̄g,m,i,r̄l,mesh)
end

function build_jixiebi(n; ks1 = 80.0, restlens1 = 50.0,
                          ks2 = .001, restlens2 = 100.0, 
                          kc = 1.2, restlenc = 160.0, pre = 10.0)
    nb = 4n+1

    len1 = 234.0; len11 = 45.33; len12 = 195.33; r1 = 100.0; r2 = 120;
    len2 = 160.0; len21 = -112.8; len22 = 47.2; r3 = 35.0
    len3 = 416.76; len31 = -233.56; len32 = -73.56; len33 = 33.21; len34 = 183.21
    len4 = 199.5; len41 = -163.99; len42 = -3.99; len43 = 18.77

    r̄g1 = [0.0, 0.0, 0.0]
    r̄g2 = [0.0, 0.0, 0.0]
    r̄g3 = [0.0, 0.0, 0.0]
    r̄g4 = [0.0, 0.0, 0.0]
    m1 = .85; m2=0.05; m3 = 1.0; m4 = 0.15

    i1 = [4575.03 0 0; 0 5114.84 0; 0 0 4575.03]
    i2 = [372.01 0 0;0 358.79 0; 0 0 372.01] 
    i3 = [6225.41 0 0; 0 6223.91 0; 0 0 6225.41] 
    i4 = [815.03 0 0; 0 108.664 0; 0 0 815.03] 
    r̄l1 = [ [len11, 0.0, -r1],
            [len11, r1*sqrt(3)/2, r1/2],
            [len11, -r1*sqrt(3)/2, r1/2],
            [len12, 0.0, -r1],
            [len12, r1*sqrt(3)/2, r1/2],
            [len12, -r1*sqrt(3)/2, r1/2],
            [len11, 0.0, r2],
            [len11, r2*sqrt(3)/2, -r2/2],
            [len11, -r2*sqrt(3)/2, -r2/2],
            [len12, 0.0, r2],
            [len12, r2*sqrt(3)/2, -r2/2],
            [len12, -r2*sqrt(3)/2, -r2/2] ]
    r̄l2 = [ [len21, 0.0, -r3],
            [len21, r3*sqrt(3)/2, r3/2],
            [len21, -r3*sqrt(3)/2, r3/2],
            [len22, 0.0, -r1],
            [len22, r1*sqrt(3)/2, r1/2],
            [len22, -r1*sqrt(3)/2, r1/2],
            [len22, 0.0, r2],
            [len22, r2*sqrt(3)/2, -r2/2],
            [len22, -r2*sqrt(3)/2, -r2/2] ]
    r̄l3 = [ [len31, 0.0, -r3],
            [len31, r3*sqrt(3)/2, r3/2],
            [len31, -r3*sqrt(3)/2, r3/2],
            [len32, 0.0, -r1],
            [len32, r1*sqrt(3)/2, r1/2],
            [len32, -r1*sqrt(3)/2, r1/2],
            [len33, 0.0, -r1],
            [len33, r1*sqrt(3)/2, r1/2],
            [len33, -r1*sqrt(3)/2, r1/2],
            [len34, 0.0, -r1],
            [len34, r1*sqrt(3)/2, r1/2],
            [len34, -r1*sqrt(3)/2, r1/2],
            [len32, 0.0, r2],
            [len32, r2*sqrt(3)/2, -r2/2],
            [len32, -r2*sqrt(3)/2, -r2/2],
            [len33, 0.0, r2],
            [len33, r2*sqrt(3)/2, -r2/2],
            [len33, -r2*sqrt(3)/2, -r2/2],
            [len34, 0.0, r2],
            [len34, r2*sqrt(3)/2, -r2/2],
            [len34, -r2*sqrt(3)/2, -r2/2] ]
    r̄l4 = [ [len41, 0.0, -r3],
            [len41, r3*sqrt(3)/2, r3/2],
            [len41, -r3*sqrt(3)/2, r3/2],
            [len42, 0.0, -r1],
            [len42, r1*sqrt(3)/2, r1/2],
            [len42, -r1*sqrt(3)/2, r1/2],
            [len42, 0.0, r2],
            [len42, r2*sqrt(3)/2, -r2/2],
            [len42, -r2*sqrt(3)/2, -r2/2] ]

    mesh1 = load("m1.STL") |> make_patch(;rot=RotZ(-π/2)) |> make_patch(;rot=RotX(pi/6))
    mesh2 = load("m2.STL") |> make_patch(;rot=RotZ(-π/2)) |> make_patch(;rot=RotX(pi/6))
    mesh3 = load("m3.STL") |> make_patch(;rot=RotZ(-π/2)) |> make_patch(;rot=RotX(pi/6))
    mesh4 = load("m4.STL") |> make_patch(;rot=RotZ(-π/2)) |> make_patch(;rot=RotX(pi/6))

    rd1 = RigidData(r̄g1, m1, i1, r̄l1, mesh1)
    rd2 = RigidData(r̄g2, m2, i2, r̄l2, mesh2)
    rd3 = RigidData(r̄g3, m3, i3, r̄l3, mesh3)
    rd4 = RigidData(r̄g4, m4, i4, r̄l4, mesh4)
    rds = reduce(vcat,[[rd1],[rd2 for _ in 1:3],reduce(vcat,[vcat([rd3],[rd2 for _ in 1:3]) for _ in 1:n-1]),[rd4]])

    function get_r̄s(rds)
        r̄s = [[0.0, 0.0, 0.0] for _ in 1:nb]
        r̄s[1][1] = -len11
        for i in 1:nb-1
            if length(rds[i].r̄l) == 12
                r̄s[i+1][1] += r̄s[i][1] + len12 - len21
            elseif length(rds[i].r̄l) == 9
                if length(rds[i+1].r̄l) == 21
                    r̄s[i+1][1] += r̄s[i][1] + len22 - len31
                elseif i == nb-1
                    r̄s[i+1][1] += r̄s[i][1] + len22 - len41
                else
                    r̄s[i+1][1] += r̄s[i][1] + len2 
                end
            elseif length(rds[i].r̄l) == 21
                r̄s[i+1][1] += r̄s[i][1] + len34 - len21
            end
        end
        return r̄s
    end
    r̄s = get_r̄s(rds)

    function rigidbody(i, rdsi, r̄si)
        if i == 1
            movable = false
            constrained = true
            ci = collect(1:12)
            Φi = Vector{Int}()
        else
            movable = true
            constrained = false
            ci = Vector{Int}()
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
        ṙo = zeros(3); ω = zeros(Float64, 3)
        prop = TR.RigidBodyProperty(
            i,
            movable,
            m,
            SMatrix{3, 3}(I),
            SVector{3}(r̄g),
            [SVector{3}(aps[i]) for i in 1:nrp],
            constrained=constrained
        )
        ro = SVector{3}(ri)

        lncs, _ = TR.NaturalCoordinates.NC1P3V(SVector{3}(ri))
        
        state = TR.RigidBodyState(prop, lncs, ri, α, ṙo, ω, ci, Φi)
        rb = TR.RigidBody(prop, state, rdsi.mesh)
    end
    rbs = [rigidbody(i, rds[i], r̄s[i]) for i in 1:nb]

    rigidbodies = TypeSortedCollection(rbs)
    numberedpoints = TR.number(rigidbodies)
    matrix_sharing_raw = Vector{Matrix{Int}}()
    matrix_sharing = matrix_sharing_raw
    indexedcoords = TR.index(rigidbodies, matrix_sharing)

    nstrings = 24n 
    ks1 = ks1; restlens1 = restlens1
    ks2 = ks2; restlens2 = restlens2
    cables = [(i<=12n) ? TR.Cable3D(i, restlens1, ks1, 0.02) : TR.Cable3D(i, restlens2, ks2, 0.02) for i in 1:nstrings]

    matrix_cnt_raw = Vector{Matrix{Int}}()

    for i in 1:n
        for j in 4i-3:4i
            s = zeros(3, nb)
            if (mod(j-1,4) == 0) & (j!=1)
                s[1, j] = 10; s[1, j+1] = -1
                s[2, j] = 11; s[2, j+1] = -2
                s[3, j] = 12; s[3, j+1] = -3
            else
                s[1, j] = 4; s[1, j+1] = -1
                s[2, j] = 5; s[2, j+1] = -2
                s[3, j] = 6; s[3, j+1] = -3
            end
            push!(matrix_cnt_raw, s)
        end
    end
    for i in 1:n
        for j in 4i-3:4i
            s = zeros(3, nb)
            if (mod(j-1,4) == 0) & (j!=1)
                s[1, j] = 10; s[1, j+1] = -4
                s[2, j] = 11; s[2, j+1] = -5
                s[3, j] = 12; s[3, j+1] = -6
            else
                s[1, j] = 4; s[1, j+1] = -4
                s[2, j] = 5; s[2, j+1] = -5
                s[3, j] = 6; s[3, j+1] = -6
            end
            push!(matrix_cnt_raw, s)
        end
    end
    matrix_cnt = reduce(vcat, matrix_cnt_raw)

    restlenc = restlenc; kc = kc; pre = pre
    cs = Vector{TR.ClusterCables}()
    restlencs = zeros(Float64, 4n)
    restlenc1 = 266.77; restlenc2 = 310.0
    for i in 1:4n
        if (mod(i,4) == 0) & (i!=4n)
            restlencs[i] = restlenc1
        elseif (mod(i-1,4) == 0) & (i!=1)
            restlencs[i] = restlenc2
        else
            restlencs[i] = restlenc
        end
    end
    c_section = StructArray([TR.CableSegment3D(i, restlencs[i], kc, prestress=pre) for i in 1:4n])
    cs1 = TR.ClusterCables(1, 4n-1, deepcopy(c_section); μ=0.01)
    cs2 = TR.ClusterCables(2, 4n-1, deepcopy(c_section); μ=0.01)
    cs3 = TR.ClusterCables(3, 4n-1, deepcopy(c_section); μ=0.01)
    cs = vcat(cs1, cs2, cs3)

    s = zeros(3, nb)
    for i in 1:nb
        if i == 1
            s[1, i] = 10; s[2, i] = 11; s[3, i] = 12
        elseif (mod(i-1,4) == 0) & (i!=nb)
            s[1, i] = 13; s[2, i] = 14; s[3, i] = 15
        else
            s[1, i] = 7; s[2, i] = 8; s[3, i] = 9
        end
    end

    matrix_cnt2 = s
    tensiles = (cables = cables, clustercables=cs)
    cc = TR.connect_and_cluster_for_jixiebi(rigidbodies, matrix_cnt, matrix_cnt2)
    cnt = TR.Connectivity(numberedpoints, indexedcoords, cc)
    tg = TR.TensegrityStructure(rigidbodies, tensiles, cnt)
    bot = TR.TensegrityRobot(tg,)

end