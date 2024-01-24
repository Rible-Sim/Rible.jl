
function build_jixiezhua(n; ks1 = 80.0, restlens1 = 50.0,
    ks2 = 1.0, restlens2 = 100.0, 
    kc = 1.2, restlenc = 160.0, pre = 383.65,
    delta = -20)
nb = 4n+1

len1 = 234.0; len11 = 45.33; len12 = 195.33; r1 = 100.0; r2 = 120;
len2 = 160.0; len21 = -112.8; len22 = 47.2; r3 = 35.0
len3 = 416.76; len31 = -233.56; len32 = -73.56; len33 = 33.21; len34 = 183.21
len4 = 199.5; len41 = -229.32; len42 = -79.32; len43 = 18.77

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

mesh1 = load("m0.STL") |> make_patch(;rot=RotZ(-π/2)) 
mesh2 = load("m2.STL") |> make_patch(;rot=RotZ(-π/2)) |> make_patch(;rot=RotX(pi/6))
mesh3 = load("m3.STL") |> make_patch(;rot=RotZ(-π/2)) |> make_patch(;rot=RotX(pi/6))
mesh4 = load("m5.STL") |> make_patch(;rot=RotZ(-π/2)) |> make_patch(;rot=RotX(pi/6))

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
r̄s[i+1][1] += r̄s[i][1] + len12 - len21 + delta
elseif length(rds[i].r̄l) == 9
if length(rds[i+1].r̄l) == 21
r̄s[i+1][1] += r̄s[i][1] + len22 - len31+ delta
elseif i == nb-1
r̄s[i+1][1] += r̄s[i][1] + len22 - len41+ delta
else
r̄s[i+1][1] += r̄s[i][1] + len2 + delta
end
elseif length(rds[i].r̄l) == 21
r̄s[i+1][1] += r̄s[i][1] + len34 - len21+ delta
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
restlencs = zeros(Float64, 5n-1)
restlenc1 = 266.77
for i in 1:5n-1
if (mod(i,5) == 0) 
restlencs[i] = restlenc1
else
restlencs[i] = restlenc
end
end
c_section = StructArray([TR.CableSegment3D(i, restlencs[i], kc, prestress=pre) for i in 1:5n-1])
c_section2 = StructArray([TR.CableSegment3D(i, restlencs[i], kc, prestress=50.0) for i in 1:5])
cs1 = TR.ClusterCables(1, 5n-2, deepcopy(c_section); μ=0.01)
cs2 = TR.ClusterCables(2, 5n-2, deepcopy(c_section); μ=0.01)
cs3 = TR.ClusterCables(3, 5n-2, deepcopy(c_section); μ=0.01)
cs4 = TR.ClusterCables(4, 4, deepcopy(c_section2); μ=0.01)
cs5 = TR.ClusterCables(5, 4, deepcopy(c_section2); μ=0.01)
cs6 = TR.ClusterCables(6, 4, deepcopy(c_section2); μ=0.01)
cs7 = TR.ClusterCables(7, 4, deepcopy(c_section2); μ=0.01)
cs8 = TR.ClusterCables(8, 4, deepcopy(c_section2); μ=0.01)
cs9 = TR.ClusterCables(9, 4, deepcopy(c_section2); μ=0.01)
cs = vcat(cs1, cs2, cs3, cs4, cs5, cs6, cs7, cs8, cs9)

s = zeros(9, nb)
for i in 1:nb
    if i == 1
        s[1, i] = 10; s[2, i] = 11; s[3, i] = 12
        elseif (mod(i-1,4) == 0) & (i!=nb)
        s[1, i] = 1319; s[2, i] = 1420; s[3, i] = 1521
        else
        s[1, i] = 7; s[2, i] = 8; s[3, i] = 9
    end
end
s[4,end-4:end] = s[1,end-4:end]
s[5,end-4:end] = s[2,end-4:end]
s[6,end-4:end] = s[3,end-4:end]
s[7,5:9] = s[1,5:9]
s[8,5:9] = s[2,5:9]
s[9,5:9] = s[3,5:9]

matrix_cnt2 = s
tensiles = (cables = cables, clustercables=cs)
cc = TR.connect_and_cluster_for_jixiebi(rigidbodies, matrix_cnt, matrix_cnt2)
cnt = TR.Connectivity(numberedpoints, indexedcoords, cc)
tg = TR.TensegrityStructure(rigidbodies, tensiles, cnt)
bot = TR.TensegrityRobot(tg,)

end

function build_jixiebi2(n; ks1 = 80.0, restlens1 = 50.0,
                          ks2 = 1.0, restlens2 = 100.0, 
                          kc = 1.2, restlenc = 160.0, pre = 383.65,
                          delta = -20)
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
                r̄s[i+1][1] += r̄s[i][1] + len12 - len21 + delta
            elseif length(rds[i].r̄l) == 9
                if length(rds[i+1].r̄l) == 21
                    r̄s[i+1][1] += r̄s[i][1] + len22 - len31+ delta
                elseif i == nb-1
                    r̄s[i+1][1] += r̄s[i][1] + len22 - len41+ delta
                else
                    r̄s[i+1][1] += r̄s[i][1] + len2 + delta
                end
            elseif length(rds[i].r̄l) == 21
                r̄s[i+1][1] += r̄s[i][1] + len34 - len21+ delta
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
    restlencs = zeros(Float64, 5n-1)
    restlenc1 = 266.77
    for i in 1:5n-1
        if (mod(i,5) == 0) 
            restlencs[i] = restlenc1
        else
            restlencs[i] = restlenc
        end
    end
    

    c_section = StructArray([TR.CableSegment3D(i, restlencs[i], kc, prestress=pre) for i in 1:5n-1])
    c_section2 = StructArray([TR.CableSegment3D(i, restlencs[i], kc, prestress=50.0) for i in 1:5])
    cs1 = TR.ClusterCables(1, 5n-2, deepcopy(c_section); μ=0.01)
    cs2 = TR.ClusterCables(2, 5n-2, deepcopy(c_section); μ=0.01)
    cs3 = TR.ClusterCables(3, 5n-2, deepcopy(c_section); μ=0.01)
    cs4 = TR.ClusterCables(4, 4, deepcopy(c_section2); μ=0.01)
    cs5 = TR.ClusterCables(5, 4, deepcopy(c_section2); μ=0.01)
    cs6 = TR.ClusterCables(6, 4, deepcopy(c_section2); μ=0.01)
    cs7 = TR.ClusterCables(7, 4, deepcopy(c_section2); μ=0.01)
    cs8 = TR.ClusterCables(8, 4, deepcopy(c_section2); μ=0.01)
    cs9 = TR.ClusterCables(9, 4, deepcopy(c_section2); μ=0.01)
    cs = vcat(cs1, cs2, cs3, cs4, cs5, cs6, cs7, cs8, cs9)
    s = zeros(9, nb)
    for i in 1:nb
        if i == 1
            s[1, i] = 10; s[2, i] = 11; s[3, i] = 12
            elseif (mod(i-1,4) == 0) & (i!=nb)
            s[1, i] = 1319; s[2, i] = 1420; s[3, i] = 1521
            else
            s[1, i] = 7; s[2, i] = 8; s[3, i] = 9
        end
    end
    s[4,end-4:end] = s[1,end-4:end]
    s[5,end-4:end] = s[2,end-4:end]
    s[6,end-4:end] = s[3,end-4:end]
    s[7,11:15] = s[1,11:15]
    s[8,11:15] = s[2,11:15]
    s[9,11:15] = s[3,11:15]

    matrix_cnt2 = s
    tensiles = (cables = cables, clustercables=cs)
    cc = TR.connect_and_cluster_for_jixiebi(rigidbodies, matrix_cnt, matrix_cnt2)
    cnt = TR.Connectivity(numberedpoints, indexedcoords, cc)
    tg = TR.TensegrityStructure(rigidbodies, tensiles, cnt)
    bot = TR.TensegrityRobot(tg,)

end

function build_jixiebi3(n; ks1 = 80.0, restlens1 = 50.0,
                          ks2 = 1.0, restlens2 = 100.0, 
                          kc = 1.2, restlenc = 160.0, pre = 383.65,
                          delta = -20)
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
                r̄s[i+1][1] += r̄s[i][1] + len12 - len21 + delta
            elseif length(rds[i].r̄l) == 9
                if length(rds[i+1].r̄l) == 21
                    r̄s[i+1][1] += r̄s[i][1] + len22 - len31+ delta
                elseif i == nb-1
                    r̄s[i+1][1] += r̄s[i][1] + len22 - len41+ delta
                else
                    r̄s[i+1][1] += r̄s[i][1] + len2 + delta
                end
            elseif length(rds[i].r̄l) == 21
                r̄s[i+1][1] += r̄s[i][1] + len34 - len21+ delta
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
    restlencs = zeros(Float64, 5n-1)
    restlenc1 = 266.77
    for i in 1:5n-1
        if (mod(i,5) == 0) 
            restlencs[i] = restlenc1
        else
            restlencs[i] = restlenc
        end
    end
    

    c_section = StructArray([TR.CableSegment3D(i, restlencs[i], kc, prestress=pre) for i in 1:5n-1])
    c_section2 = StructArray([TR.CableSegment3D(i, restlencs[i], kc, prestress=50.0) for i in 1:5])
    cs1 = TR.ClusterCables(1, 5n-2, deepcopy(c_section); μ=0.01)
    cs2 = TR.ClusterCables(2, 5n-2, deepcopy(c_section); μ=0.01)
    cs3 = TR.ClusterCables(3, 5n-2, deepcopy(c_section); μ=0.01)
    cs4 = TR.ClusterCables(4, 4, deepcopy(c_section2); μ=0.01)
    cs5 = TR.ClusterCables(5, 4, deepcopy(c_section2); μ=0.01)
    cs6 = TR.ClusterCables(6, 4, deepcopy(c_section2); μ=0.01)
    cs7 = TR.ClusterCables(7, 4, deepcopy(c_section2); μ=0.01)
    cs8 = TR.ClusterCables(8, 4, deepcopy(c_section2); μ=0.01)
    cs9 = TR.ClusterCables(9, 4, deepcopy(c_section2); μ=0.01)
    cs = vcat(cs1, cs2, cs3, cs4, cs5, cs6, cs7, cs8, cs9)
    s = zeros(9, nb)
    for i in 1:nb
        if i == 1
            s[1, i] = 10; s[2, i] = 11; s[3, i] = 12
            elseif (mod(i-1,4) == 0) & (i!=nb)
            s[1, i] = 1319; s[2, i] = 1420; s[3, i] = 1521
            else
            s[1, i] = 7; s[2, i] = 8; s[3, i] = 9
        end
    end
    s[4,2:6] = s[1,2:6]
    s[5,2:6] = s[2,2:6]
    s[6,2:6] = s[3,2:6]
    s[7,6:10] = s[1,6:10]
    s[8,6:10] = s[2,6:10]
    s[9,6:10] = s[3,6:10]

    matrix_cnt2 = s
    tensiles = (cables = cables, clustercables=cs)
    cc = TR.connect_and_cluster_for_jixiebi(rigidbodies, matrix_cnt, matrix_cnt2)
    cnt = TR.Connectivity(numberedpoints, indexedcoords, cc)
    tg = TR.TensegrityStructure(rigidbodies, tensiles, cnt)
    bot = TR.TensegrityRobot(tg,)

end