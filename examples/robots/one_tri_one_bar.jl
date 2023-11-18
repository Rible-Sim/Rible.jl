
function one_tri_one_bar(;k=nothing)
    a = 0.255 #tetrahedra base triangle length
    h = 0.39
    l = 2.0 #legs lengths
    θ = 2π/4 #triangle angle
    d = 0.72
    α = π/3
    R = RotY(-π/2)
    P = Ref(R).*vcat(
        [
            [
                2a/sqrt(2)*cos(i*θ+θ/2),
                2a/sqrt(2)*sin(i*θ+θ/2),
                -h + 0.01
            ] 
            for i = 0:3
        ],
        [
            [
                2a/sqrt(3)*cos((i+0.5)*θ),
                2a/sqrt(3)*sin((i+0.5)*θ),
                h
            ] 
            for i = -1:1
        ],
        [
            RotXY(-π/18,π/24)*[0.0,0.0,d],
            [0.0,0.0,0.0],
        ]
    ) .|> SVector{3}
    
    # scatter(P)
    
    base = make_3d_tri(
        1,
        P[[8,1,2,3,4]],
        zero(P[8]),
        R,
        P[8];
        contactable = true,
        visible = true,
        ci = collect(1:12),
        cstr_idx = Int[],
        radius = 0.04,
        loadmesh=true
    )

    bar = make_3d_bar(
        2,
        P[end],
        P[end-1],;		
        radius_ratio = 1/40,
        ci = Int[],
        loadmesh2 = true,
    )
    

    # # #
    nb = 2
    rbs = TypeSortedCollection([base,bar,])
    numberedpoints = RB.number(rbs)
    indexedcoords = RB.index(rbs,)
    
    # 
    cnt_matrix_elas = ElasticArray{Int}(undef, nb, 0)
    cnt_circular = CircularArray([2,3,4])
    # left
    for i = 1:4	
        row = zeros(Int,nb)
        row[1] = -(i+1)
        row[2] = 2
        append!(cnt_matrix_elas,row)
    end

    display(cnt_matrix_elas)
    cnt_matrix = Matrix(cnt_matrix_elas')
    ncables = size(cnt_matrix,1)
    if k isa Nothing
        cables = [RB.Cable3D(i,0.0,100.0,0.0;slack=false) for i = 1:ncables]
    else
        cables = [RB.Cable3D(i,0.0,k[i],0.0;slack=false) for i = 1:ncables]
    end
    tensiles = (cables = cables,)
    connected = RB.connect(rbs,cnt_matrix)

    cst1 = RB.PinJoint(
        1,
        RB.Hen2Egg(
            1,
            RB.ID(base,1),
            RB.ID(bar,1),
        )
    )
    
    jointedmembers = RB.join((cst1,),indexedcoords)

    cnt = RB.Connectivity(
        numberedpoints,
        indexedcoords,
        @eponymtuple(connected,),
        jointedmembers
    )

    st = RB.Structure(rbs,tensiles,cnt)
    bot = RB.Robot(st,)
end