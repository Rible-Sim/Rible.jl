
function tower(;k=nothing)
    a = 1.0 #tetrahedra base triangle length
    h = 0.6 #tegrahedra height
    l = 2.0 #legs lengths
    θ = 2π/3 #triangle angle
    d = 1.0
    α = π/3
    P = vcat(
        [
            [
                2a/sqrt(3)*cos(i*θ),
                2a/sqrt(3)*sin(i*θ),
                -0.5h
            ] 
            for i = -1:1
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
            [0.0,0.0,d],
            [0.0,0.0,0.0],
        ]
    ) .|> SVector{3}
    
    # scatter(P)
    
    base = make_3d_tri(
        1,
        P[[8,1,2,3]],
        zero(P[8]),
        RotX(0.0),
        P[8];
        contactable = true,
        visible = true,
        ci = collect(1:12),
        cstr_idx = Int[],
        radius = 0.04,
    )

    bar = make_3d_bar(
        2,
        P[end],
        P[end-1],;		
        radius_ratio = 1/40,
        ci = Int[]
    )

    
    top = make_3d_tri(
        3,
        P[[7,4,5,6]],
        zero(P[7]),
        RotX(0.0),
        P[7];
        radius = 0.04,
        contactable = true,
        visible = true,
    )

    # # #
    nb = 3
    rbs = TypeSortedCollection([base,bar,top])
    numbered = RB.number(rbs)
    indexed = RB.index(rbs,)
    
    # 
    cnt_matrix_elas = ElasticArray{Int}(undef, nb, 0)
    cnt_circular = CircularArray([2,3,4])
    # left
    for i = 1:3	
        row = zeros(Int,nb)
        row[1] = -(i+1)
        row[2] = 2
        append!(cnt_matrix_elas,row)
    end
    for i = 1:3	
        row = zeros(Int,nb)
        row[2] = 1
        row[3] = -(i+1)
        append!(cnt_matrix_elas,row)
    end
    for i = 1:3	
        row = zeros(Int,nb)
        row[1] = cnt_circular[i]
        row[3] = -cnt_circular[i]
        append!(cnt_matrix_elas,row)
        row = zeros(Int,nb)
        row[1] = cnt_circular[i]
        row[3] = -cnt_circular[i-1]
        append!(cnt_matrix_elas,row)
    end	

    display(cnt_matrix_elas)
    cnt_matrix = Matrix(cnt_matrix_elas')
    ncables = size(cnt_matrix,1)
    if k isa Nothing
        cables = [RB.DistanceSpringDamper3D(0.0,100.0,0.0;slack=false) for i = 1:ncables]
    else
        cables = [RB.DistanceSpringDamper3D(0.0,k[i],0.0;slack=false) for i = 1:ncables]
    end
    apparatuses = (cables = cables,)
    connected = RB.connect(rbs,cnt_matrix)

    cst1 = RB.PinJoint(
        1,
        RB.Hen2Egg(
            1,
            RB.ID(base,1),
            RB.ID(bar,1),
        )
    )
    
    cst2 = RB.PinJoint(
        2,
        RB.Hen2Egg(
            2,
            RB.ID(bar,2),
            RB.ID(top,1)
        )
    )

    jointedmembers = RB.join((cst1,cst2,),indexed)

    cnt = RB.Connectivity(
            numbered,
            indexed,
            @eponymtuple(connected,),
            jointedmembers
        )

    st = RB.Structure(rbs,apparatuses,cnt)
    bot = RB.Robot(st,)
end