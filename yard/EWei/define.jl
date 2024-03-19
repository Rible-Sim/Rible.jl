function RigidData(mass_locus,m,i,r̄l,mesh=nothing)
    @eponymtuple(mass_locus,m,i,r̄l,mesh)
end

"""
不含滑动绳索
"""
function BuildTail()
    nb = 24
    r̄gC = [44.71/2, 0.0]
    r̄gA = [50/2, 0.0]
    r̄gB = [50/2, 0.0]
    r̄g1 = [0.0, 0.0]
    r̄g2 = [0.0, 0.0]
    r̄g3 = [0.0, 0.0]
    r̄g4⁺ = [50-7.788, 0.0]
    r̄g4 = [35-3.835, 0.0]
    r̄g5 = [35-1.841, 0.0]
    r̄g6 = [25-0.079, 0.0]
    r̄g7 = [20+0.688, 0.0]
    r̄g8 = [20+0.766, 0.0]
    r̄g9 = [20+0.865, 0.0]
    r̄g10 = [20+1.074, 0.0]

    mC = 1.; mA = 22.; mB = 26.; m1 = 60.; m2 = 59.
    m3 = 54.; m4⁺ = 78.; m4 = 48.; m5 = 48.; m6 = 45.
    m7 = 40.; m8 = 36.; m9 = 32.; m10 = 40.
    iC = 1.; iA = 6.071; iB = 8.646; i1 = 60.660;
    i2 = 54.204; i3 = 41.680; i4⁺ = 68.548; i4 = 39.220
    i5 = 47.819; i6 = 53.351; i7 = 45.136; i8 = 36.216
    i9 = 24.286; i10 = 13.179

    r̄lC = [0.0, 0.0, 0.0, 0.0]; r̄lA = deepcopy(r̄gA); r̄lB = deepcopy(r̄gB)
    r̄l1 = 40; r̄l2 = 40; r̄l3 = 40
    r̄l4⁺ = [50, 32, 15]; r̄l4 = [35, 32, 15]
    r̄l5 = [35, 28, 15]; r̄l6 = [25, 27, 15]
    r̄l7 = [20, 21, 15]; r̄l8 = [20, 20, 15]
    r̄l9 = [20, 15, 15]; r̄l10 = [20, 15, 15]

    function Get_apsAB(r̄, mass_locus)
        ri = r̄
        rj = r̄ + 2* [mass_locus[2], mass_locus[1]]
        aps = [[0.0, 0.0], 2r̄g, [15.1, 0.0]]
        return ri, rj, aps
    end

    function Get_apsC()
        ri = [0.0, 0.0]
        rj = [0.0, 44.71]
        aps = [[44.71, 0.0], [0.0, -40.0], [0.0, 40.0], [15.1, 0.0]]
        return ri, rj, aps
    end

    function Get_aps123(r̄, r̄l)
        ri = r̄
        rj = r̄ + [r̄l, 0.0]
        aps = [[r̄l, 0.0], -[r̄l, 0.0], [0.0, 0.0]]
        return ri, rj, aps
    end

    function Get_apsOther(r̄, r̄ls)
        down, right, up = r̄ls
        ri = r̄ - [0.0, down]
        rj = r̄ + [0.0, up]
        aps = [[0.0, 0.0], [down, right], [down, -right], [up+down, 0.0], [down, 0.0]]
        return ri, rj, aps
    end

    function Get_aps(r̄, r̄l)
        return @match length(r̄l) begin
            2 => Get_apsAB(r̄, r̄l)
            1 => Get_aps123(r̄, r̄l)
            3 => Get_apsOther(r̄, r̄l)
            4 => Get_apsC()
        end
    end

    rds = [
        RigidData(r̄gC, mC, iC, r̄lC),
        RigidData(r̄g1, m1, i1, r̄l1),
        RigidData(r̄gB, mB, iB, r̄lB),
        RigidData(r̄g1, m1, i1, r̄l1),
        RigidData(r̄gA, mA, iA, r̄lA),
        RigidData(r̄g1, m1, i1, r̄l1),
        RigidData(r̄gB, mB, iB, r̄lB),
        RigidData(r̄g1, m1, i1, r̄l1),
        RigidData(r̄gA, mA, iA, r̄lA),
        RigidData(r̄g2, m2, i2, r̄l2),
        RigidData(r̄gB, mB, iB, r̄lB),
        RigidData(r̄g2, m2, i2, r̄l2),
        RigidData(r̄gA, mA, iA, r̄lA),
        RigidData(r̄g2, m2, i2, r̄l2),
        RigidData(r̄gB, mB, iB, r̄lB),
        RigidData(r̄g3, m3, i3, r̄l3),
        RigidData(r̄g4⁺, m4⁺, i4⁺, r̄l4⁺),
        RigidData(r̄g4, m4, i4, r̄l4),
        RigidData(r̄g5, m5, i5, r̄l5),
        RigidData(r̄g6, m6, i6, r̄l6),
        RigidData(r̄g7, m7, i7, r̄l7),
        RigidData(r̄g8, m8, i8, r̄l8),
        RigidData(r̄g9, m9, i9, r̄l9),
        RigidData(r̄g10, m10, i10, r̄l10)]

    function Computer̄s(rds)
        n = length(rds)
        r̄s = [[0.0, 0.0] for _ in 1:n]
        for i in 1:n
            if length(rds[i].r̄l) == 1
                r̄s[i][2] += r̄s[i-1][2]
            elseif length(rds[i].r̄l) == 2
                len = rds[i].r̄l[1] * 2
                r̄s[i+1][2] += r̄s[i][2] + len
                if i != 1
                    r̄s[i][2] += r̄s[i-1][2]
                end
            elseif length(rds[i].r̄l) == 3
                len = rds[i].r̄l[1]
                r̄s[i][2] += r̄s[i-1][2] + len
                if i != n
                    r̄s[i+1][2] += rds[i].r̄l[3]
                end
            elseif length(rds[i].r̄l) == 4
                r̄s[i+1][2] += 44.71
            end
        end
        return r̄s
    end

    r̄s = Computer̄s(rds)
    rigidIndex = vcat([1], [2i for i in 1:8], [i for i in 17:24])
    barIndex = collect(3:2:13)
    function rigidbody(i, rdsi, r̄si; rigidIndex=rigidIndex, barIndex=barIndex)
        if i == 1
            contactable = false
            visible = true
            ci = collect(1:6)
            cstr_idx = Int[]
        else
            contactable = true
            if i in 2:3
                visible = true
                ci = collect(1:2)
                if i == 2
                    cstr_idx = collect(1:3)
                else
                    cstr_idx = [1]
                end
            else
                visible = true
                ci = Int[]
                if i in rigidIndex
                    cstr_idx = collect(1:3)
                else
                    cstr_idx = [1]
                end
            end
        end

        ri, rj, aps = Get_aps(r̄si, rdsi.r̄l)
        u = rj - ri
        α = atan(u[2], u[1])
        mass_locus = rdsi.mass_locus
        m = rdsi.m*1e-3
        ## I = rdsi.i
        @show m
        nrp = length(aps)
        ṙo = zeros(2); ω = 0.0
        prop = RB.RigidBodyProperty(
            i,
            contactable,
            m,
            SMatrix{2, 2}([
                m 0
                0 m
            ]),
            SVector{2}(mass_locus),
            [SVector{2}(aps[i]) for i in 1:nrp],
            visible=visible
        )
        ro = SVector{2}(ri)
        if i in rigidIndex
            @show i, "NC2P1V", ri, rj
            nmcs = RB.NCF.NC2P1V(SVector{2}(ri), SVector{2}(rj), ro, RB.rotation_matrix(α))
        else
            @show i, "NC2D2P"
            nmcs = RB.NCF.NC2D2P(SVector{2}(ri), SVector{2}(rj), ro, RB.rotation_matrix(α))
        end
        state = RB.RigidBodyState(prop, ri, α, ṙo, ω, )
        coords = RB.NonminimalCoordinates(nmcs, ci, cstr_idx)
        body = RB.RigidBody(prop, state,coords)
    end

    rbs = [rigidbody(i, rds[i], r̄s[i]) for i in 1:nb]
    rigdibodies = TypeSortedCollection(rbs)

    nstrings = 8 * 4 + 8 * 2
    ks = zeros(nstrings); restlens = zeros(nstrings)
    for i in 1:4*8
        j = i % 4
        restlens[i] = ifelse(j in [1, 0], 25, 20)
        ks[i] = ifelse(j in [1, 0], 0.3, 0.5)
    end
    for i in 4*8+1:nstrings
        restlens[i] = 25
        ks[i] = 0.3
    end
    spring_dampers = [RB.DistanceSpringDamper2D( restlens[i], ks[i], 0.0) for i in 1:nstrings]

    matrix_cnt_raw = Vector{Matrix{Int}}()
    s = zeros(4, nb)
    s[1, 1] = 3; s[1, 2] = -2
    s[2, 1] = 4; s[2, 2] = -2
    s[3, 1] = 4; s[3, 2] = -1
    s[4, 1] = 2; s[4, 2] = -1
    push!(matrix_cnt_raw, s)
    for i in 1:7
        s = zeros(4, nb)
        s[1, 2i] = 2; s[1, 2i+2] = -2
        s[2, 2i+1] = 3; s[2, 2i+2] = -2
        s[3, 2i+1] = 3; s[3, 2i+2] = -1
        s[4, 2i] = 1; s[4, 2i+2] = -1
        push!(matrix_cnt_raw, s)
    end
    s = zeros(2, nb)
    s[1, 16] = 1; s[1, 17] = -3
    s[2, 16] = 2; s[2, 17] = -2
    push!(matrix_cnt_raw, s)
    for i in 17:23
        s = zeros(2, nb)
        s[1, i] = 3; s[1, i+1] = -3
        s[2, i] = 2; s[2, i+1] = -2
        push!(matrix_cnt_raw, s)
    end
    connecting_matrix = reduce(vcat, matrix_cnt_raw)
    spring_dampers = RB.connect(rigdibodies, spring_dampers;connecting_matrix)
    apparatuses = TypeSortedCollection(
        spring_dampers
    )

    numbered = RB.number(rigdibodies, apparatuses)
    
    matrix_sharing_raw = Vector{Matrix{Int}}()
    for i in 1:nb-1
        s = zeros(2, nb)
        if i == 1
            s[1:2, 1] = 3:4
            s[1:2, 2] = 1:2
        elseif i in 2:2:16
            s[1:2, i] = 1:2
            s[1:2, i+1] = 1:2
        else
            s[1:2, i] = 3:4
            s[1:2, i+1] = 1:2
        end
        push!(matrix_sharing_raw, s)
    end
    sharing_matrix = reduce(vcat, matrix_sharing_raw)
    indexed = RB.index(rigdibodies, apparatuses; sharing_matrix)
    cnt = RB.Connectivity(
        numbered, 
        indexed,
    )
    st = RB.Structure(rigdibodies, apparatuses, cnt)
    gauges = Int[]
    actuators = Int[]
    hub = RB.ControlHub(
        st,
        gauges,
        actuators,
        RB.Coalition(st,gauges,actuators)
    )
    bot = RB.Robot(st,)
end

"""
包含滑动绳索的模型构建函数。
type=1,2,3,4分别对应四种不同的滑动绳索安装方式
"""
function BuildTail(type; β=1.0, μ=0.02)
    nb = 24
    r̄gC = [44.71/2, 0.0]
    r̄gA = [50/2, 0.0]
    r̄gB = [50/2, 0.0]
    r̄g1 = [0.0, 0.0]
    r̄g2 = [0.0, 0.0]
    r̄g3 = [0.0, 0.0]
    r̄g4⁺ = [50-7.788, 0.0]
    r̄g4 = [35-3.835, 0.0]
    r̄g5 = [35-1.841, 0.0]
    r̄g6 = [25-0.079, 0.0]
    r̄g7 = [20+0.688, 0.0]
    r̄g8 = [20+0.766, 0.0]
    r̄g9 = [20+0.865, 0.0]
    r̄g10 = [20+1.074, 0.0]

    mC = 1.; mA = 22.; mB = 26.; m1 = 60.; m2 = 59.
    m3 = 54.; m4⁺ = 78.; m4 = 48.; m5 = 48.; m6 = 45.
    m7 = 40.; m8 = 36.; m9 = 32.; m10 = 40.
    iC = 1.; iA = 6.071; iB = 8.646; i1 = 60.660;
    i2 = 54.204; i3 = 41.680; i4⁺ = 68.548; i4 = 39.220
    i5 = 47.819; i6 = 53.351; i7 = 45.136; i8 = 36.216
    i9 = 24.286; i10 = 13.179

    r̄lC = [0.0, 0.0, 0.0, 0.0]; r̄lA = deepcopy(r̄gA); r̄lB = deepcopy(r̄gB)
    r̄l1 = 40; r̄l2 = 40; r̄l3 = 40
    r̄l4⁺ = [50, 32, 15]; r̄l4 = [35, 32, 15]
    r̄l5 = [35, 28, 15]; r̄l6 = [25, 27, 15]
    r̄l7 = [20, 21, 15]; r̄l8 = [20, 20, 15]
    r̄l9 = [20, 15, 15]; r̄l10 = [20, 15, 15]

    function Get_apsAB(r̄, mass_locus)
        ri = r̄
        rj = r̄ + 2* [mass_locus[2], mass_locus[1]]
        aps = [[0.0, 0.0], 2r̄g, [15.1, 0.0]]
        return ri, rj, aps
    end

    function Get_apsC()
        ri = [0.0, 0.0]
        rj = [0.0, 44.71]
        aps = [[44.71, 0.0], [0.0, -40.0], [0.0, 40.0], [15.1, 0.0]]
        return ri, rj, aps
    end

    function Get_aps123(r̄, r̄l)
        ri = r̄
        rj = r̄ + [r̄l, 0.0]
        aps = [[r̄l, 0.0], -[r̄l, 0.0], [0.0, 0.0]]
        return ri, rj, aps
    end

    function Get_apsOther(r̄, r̄ls)
        down, right, up = r̄ls
        ri = r̄ - [0.0, down]
        rj = r̄ + [0.0, up]
        aps = [[0.0, 0.0], [down, right], [down, -right], [up+down, 0.0], [down, 0.0]]
        return ri, rj, aps
    end

    function Get_aps(r̄, r̄l)
        return @match length(r̄l) begin
            2 => Get_apsAB(r̄, r̄l)
            1 => Get_aps123(r̄, r̄l)
            3 => Get_apsOther(r̄, r̄l)
            4 => Get_apsC()
        end
    end
    asbMid = load("STL/装配体-中间部分.STL") |> make_patch(;trans=[0.0,0,0],rot=RotZ(deg2rad(-77)))
    crSec1 = load("STL/截面形状-1.STL") |> make_patch(;trans=[0.0,-5,0]) |> (x)->merge([x,asbMid])
    crSec2 = load("STL/截面形状-2.STL") |> make_patch(;trans=[0.0,-5,0]) |> (x)->merge([x,asbMid])
    crSec3 = load("STL/截面形状-3.STL") |> make_patch(;trans=[0.0,-5,0]) |> (x)->merge([x,asbMid])
    crSec4f = load("STL/截面形状-4-随动段.STL") |> make_patch(;trans=[40,0,0.],rot=RotZ(-π/2))
    crSec5f = load("STL/截面形状-5-随动段.STL") |> make_patch(;trans=[40,0,0.],rot=RotZ(-π/2))
    crSec6f = load("STL/截面形状-6-随动段.STL") |> make_patch(;trans=[30,0,0.],rot=RotZ(-π/2))
    crSec7f = load("STL/截面形状-7-随动段.STL") |> make_patch(;trans=[25,0,0.],rot=RotZ(-π/2))
    crSec8f = load("STL/截面形状-8-随动段.STL") |> make_patch(;trans=[25,0,0.],rot=RotZ(-π/2))
    crSec9f = load("STL/截面形状-9-随动段.STL") |> make_patch(;trans=[25,0,0.],rot=RotZ(-π/2))
    crSec10f = load("STL/截面形状-10-随动段.STL") |> make_patch(;trans=[30,0,0.],rot=RotZ(-π/2))

    asbMot = load("STL/装配体-电机.STL") |> make_patch(;trans=[45.0,0,0],rot=RotZ(π))

    asbB32 = load("STL/装配体-中间梁6.0-32.STL") |> make_patch(;trans=[25.0,0,0])
    asbB40 = load("STL/装配体-中间梁6.0-40.STL") |> make_patch(;trans=[25.0,0,0])
    asbLF = load("STL/装配体-主动-随动.STL") |> make_patch(;trans=[55,0,0.],rot=RotZ( π/2))

    rds = [
        RigidData(r̄gC, mC, iC, r̄lC, asbMot), #1
        RigidData(r̄g1, m1, i1, r̄l1, crSec1), #2
        RigidData(r̄gB, mB, iB, r̄lB, asbB32), #3
        RigidData(r̄g1, m1, i1, r̄l1, crSec1), #4
        RigidData(r̄gA, mA, iA, r̄lA, asbB40), #5
        RigidData(r̄g1, m1, i1, r̄l1, crSec1), #6
        RigidData(r̄gB, mB, iB, r̄lB, asbB32), #7
        RigidData(r̄g1, m1, i1, r̄l1, crSec1), #8
        RigidData(r̄gA, mA, iA, r̄lA, asbB40), #9
        RigidData(r̄g2, m2, i2, r̄l2, crSec2), #10
        RigidData(r̄gB, mB, iB, r̄lB, asbB32), #11
        RigidData(r̄g2, m2, i2, r̄l2, crSec2), #12
        RigidData(r̄gA, mA, iA, r̄lA, asbB40), #13
        RigidData(r̄g2, m2, i2, r̄l2, crSec2), #14
        RigidData(r̄gB, mB, iB, r̄lB, asbB32), #15
        RigidData(r̄g3, m3, i3, r̄l3, crSec3), #16
        RigidData(r̄g4⁺, m4⁺, i4⁺, r̄l4⁺, asbLF), #17
        RigidData(r̄g4, m4, i4, r̄l4, crSec4f), #18
        RigidData(r̄g5, m5, i5, r̄l5, crSec5f), #19
        RigidData(r̄g6, m6, i6, r̄l6, crSec6f), #20
        RigidData(r̄g7, m7, i7, r̄l7, crSec7f), #21
        RigidData(r̄g8, m8, i8, r̄l8, crSec8f), #22
        RigidData(r̄g9, m9, i9, r̄l9, crSec9f), #23
        RigidData(r̄g10, m10, i10, r̄l10, crSec10f) #24
    ]

    function Computer̄s(rds)
        n = length(rds)
        r̄s = [[0.0, 0.0] for _ in 1:n]
        for i in 1:n
            if length(rds[i].r̄l) == 1
                r̄s[i][2] += r̄s[i-1][2]
            elseif length(rds[i].r̄l) == 2
                len = rds[i].r̄l[1] * 2
                r̄s[i+1][2] += r̄s[i][2] + len
                if i != 1
                    r̄s[i][2] += r̄s[i-1][2]
                end
            elseif length(rds[i].r̄l) == 3
                len = rds[i].r̄l[1]
                r̄s[i][2] += r̄s[i-1][2] + len
                if i != n
                    r̄s[i+1][2] += rds[i].r̄l[3]
                end
            elseif length(rds[i].r̄l) == 4
                r̄s[i+1][2] += 44.71
            end
        end
        return r̄s
    end

    r̄s = Computer̄s(rds)
    rigidIndex = vcat([1], [2i for i in 1:8], [i for i in 17:24])
    barIndex = collect(3:2:13)
    function rigidbody(i, rdsi, r̄si; rigidIndex=rigidIndex, barIndex=barIndex)
        if i == 1
            contactable = false
            visible = true
            ci = collect(1:6)
            cstr_idx = Int[]
        else
            contactable = true
            if i in 2:3
                visible = true
                ci = collect(1:2)
                if i == 2
                    cstr_idx = collect(1:3)
                else
                    cstr_idx = [1]
                end
            else
                visible = true
                ci = Int[]
                if i in rigidIndex
                    cstr_idx = collect(1:3)
                else
                    cstr_idx = [1]
                end
            end
        end

        ri, rj, aps = Get_aps(r̄si, rdsi.r̄l)
        u = rj - ri
        α = atan(u[2], u[1])
        mass_locus = rdsi.mass_locus
        m = rdsi.m*1e-3
        I = rdsi.i
        nrp = length(aps)
        ṙo = zeros(2); ω = 0.0
        prop = RB.RigidBodyProperty(
            i,
            contactable,
            m,
            SMatrix{2, 2}([
                0.99I 0
                0 0.01I
            ]),
            SVector{2}(mass_locus),
            [SVector{2}(aps[i]) for i in 1:nrp],
            visible=visible
        )
        ro = SVector{2}(ri)
        if i in rigidIndex
            nmcs = RB.NCF.NC2P1V(SVector{2}(ri), SVector{2}(rj), ro, α)
        else
            nmcs = RB.NCF.NC2D2P(SVector{2}(ri), SVector{2}(rj), ro, α)
        end
        state = RB.RigidBodyState(prop, nmcs, ri, α, ṙo, ω, ci, cstr_idx)
        body = RB.RigidBody(prop, state, rdsi.mesh)
    end

    rbs = [rigidbody(i, rds[i], r̄s[i]) for i in 1:nb]
    rigdibodies = TypeSortedCollection(rbs)

    nstrings = 8 * 4 + 8 * 2
    ks = zeros(nstrings); restlens = zeros(nstrings)
    for i in 1:4*8
        j = i % 4
        restlens[i] = ifelse(j in [1, 0], 25, 20)
        ks[i] = ifelse(j in [1, 0], 0.3, 0.5)
    end
    for i in 4*8+1:nstrings
        restlens[i] = 25
        ks[i] = 0.3
    end
    spring_dampers = [RB.DistanceSpringDamper2D( restlens[i], ks[i], 0.0) for i in 1:nstrings]
    @match type begin
        1 => begin
            lens = [44.71, 50, 50]
            @show pre = [lens[i]*(1-β)*0.2 for i in 1:3]
            c_section = StructArray([RB.DistanceSpringDamperSegment(i, lens[i], 0.2, prestress=pre[i]) for i in 1:3])
            cs1 = RB.ClusterDistanceSpringDampers(1, 2, deepcopy(c_section); μ=μ)
            cs2 = RB.ClusterDistanceSpringDampers(2, 2, deepcopy(c_section); μ=μ)
            apparatuses = (spring_dampers=spring_dampers, clustercables=[cs1, cs2])
            s = zeros(Int, 2, nb)
            s[1, 1] = 3; s[1, 2] = 2; s[1, 4] = 2; s[1, 6] = 2
            s[2, 1] = 2; s[2, 2] = 1; s[2, 4] = 1; s[2, 6] = 1
            matrix_cnt2 = s
        end
        2 => begin
            lens = vcat([44.71], [50 for _ in 1:7])
            pre = [lens[i]*(1-β)*0.2 for i in 1:8]
            c_section = StructArray([RB.DistanceSpringDamperSegment(i, lens[i], 0.2, prestress=pre[i]) for i in 1:8])
            cs1 = RB.ClusterDistanceSpringDampers(1, 7, deepcopy(c_section); μ=μ)
            cs2 = RB.ClusterDistanceSpringDampers(2, 7, deepcopy(c_section); μ=μ)
            apparatuses = (spring_dampers=spring_dampers, clustercables=[cs1, cs2])
            s = zeros(Int, 2, nb)
            s[1, 1] = 3; s[2, 1] = 2
            for i in 1:8
                s[1, 2i] = 2
                s[2, 2i] = 1
            end
            matrix_cnt2 = s
        end
        3 => begin
            lens1 = [44.71, 50]; lens2 = [50.0, 50.0]
            pre1 = [lens1[i]*(1-β)*0.2 for i in 1:2]
            pre2 = [lens2[i]*(1-β)*0.2 for i in 1:2]
            c_section1 = StructArray([RB.DistanceSpringDamperSegment(i, lens1[i], 0.2, prestress=pre1[i]) for i in 1:2])
            c_section2 = StructArray([RB.DistanceSpringDamperSegment(i, lens2[i], 0.2, prestress=pre2[i]) for i in 1:2])
            cs1 = [RB.ClusterDistanceSpringDampers(i, 1, deepcopy(c_section1); μ=μ) for i in 1:2]
            cs2 = [RB.ClusterDistanceSpringDampers(i, 1, deepcopy(c_section2); μ=μ) for i in 3:8]
            apparatuses = (spring_dampers=spring_dampers, clustercables=vcat(cs1, cs2))
            s = zeros(Int, 8, nb)
            s[1, 1] = 3; s[1, 2] = 2; s[1, 4] = 2;
            s[2, 1] = 2; s[2, 2] = 1; s[2, 4] = 1;
            for i in 1:3
                for j in 1:3
                    s[2j+1, 4j+2i-2] = 2
                    s[2j+2, 4j+2i-2] = 1
                end
            end
            matrix_cnt2 = s
        end
        4 => begin
            lens1 = [44.71, 50, 50]; lens2 = [50.0, 50.0, 50.0]
            pre1 = [lens1[i]*(1-β)*0.2 for i in 1:3]
            pre2 = [lens2[i]*(1-β)*0.2 for i in 1:3]
            c_section1 = StructArray([RB.DistanceSpringDamperSegment(i, lens1[i], 0.2, prestress=pre1[i]) for i in 1:3])
            c_section2 = StructArray([RB.DistanceSpringDamperSegment(i, lens2[i], 0.2, prestress=pre2[i]) for i in 1:3])
            cs1 = [RB.ClusterDistanceSpringDampers(i, 2, deepcopy(c_section1); μ=μ) for i in 1:2]
            cs2 = [RB.ClusterDistanceSpringDampers(i, 2, deepcopy(c_section2); μ=μ) for i in 3:6]
            apparatuses = (spring_dampers=spring_dampers, clustercables=vcat(cs1, cs2))
            s = zeros(Int, 6, nb)
            s[1, 1] = 3; s[1, 2] = 2; s[1, 4] = 2; s[1, 6] = 2
            s[2, 1] = 2; s[2, 2] = 1; s[2, 4] = 1; s[2, 6] = 1
            for i in 1:4
                for j in 1:2
                    s[2j+1, 4j+2i-2] = 2
                    s[2j+2, 4j+2i-2] = 1
                end
            end
            matrix_cnt2 = s
        end
    end
    # acs = []
    # hub = (actuators = acs, )

    matrix_cnt_raw = Vector{Matrix{Int}}()
    s = zeros(4, nb)
    s[1, 1] = 3; s[1, 2] = -2
    s[2, 1] = 4; s[2, 2] = -2
    s[3, 1] = 4; s[3, 2] = -1
    s[4, 1] = 2; s[4, 2] = -1
    push!(matrix_cnt_raw, s)
    for i in 1:7
        s = zeros(4, nb)
        s[1, 2i] = 2; s[1, 2i+2] = -2
        s[2, 2i+1] = 3; s[2, 2i+2] = -2
        s[3, 2i+1] = 3; s[3, 2i+2] = -1
        s[4, 2i] = 1; s[4, 2i+2] = -1
        push!(matrix_cnt_raw, s)
    end
    s = zeros(2, nb)
    s[1, 16] = 1; s[1, 17] = -3
    s[2, 16] = 2; s[2, 17] = -2
    push!(matrix_cnt_raw, s)
    for i in 17:23
        s = zeros(2, nb)
        s[1, i] = 3; s[1, i+1] = -3
        s[2, i] = 2; s[2, i+1] = -2
        push!(matrix_cnt_raw, s)
    end
    connecting_matrix = reduce(vcat, matrix_cnt_raw)

    cc = RB.connect_and_cluster(rigdibodies, connecting_matrix, matrix_cnt2)

    matrix_sharing_raw = Vector{Matrix{Int}}()
    for i in 1:nb-1
        s = zeros(2, nb)
        if i == 1
            s[1:2, 1] = 3:4
            s[1:2, 2] = 1:2
        elseif i in 2:2:16
            s[1:2, i] = 1:2
            s[1:2, i+1] = 1:2
        else
            s[1:2, i] = 3:4
            s[1:2, i+1] = 1:2
        end
        push!(matrix_sharing_raw, s)
    end
    matrix_sharing = reduce(vcat, matrix_sharing_raw)

    numbered = RB.number(rigdibodies, apparatuses;)
    indexed = RB.index(rigdibodies, apparatuses; matrix_sharing)

    cnt = RB.Connectivity(numbered, indexed)

    st = RB.Structure(rigdibodies, apparatuses, cnt)

    gauges = Int[]
    actuators = Int[]
    hub = RB.ControlHub(
        st,
        gauges,
        actuators,
        RB.Coalition(st,gauges,actuators)
    )

    bot = RB.Robot(st,hub)
end
