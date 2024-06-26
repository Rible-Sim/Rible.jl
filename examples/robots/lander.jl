function lander(;k=nothing)
    a = 1.0 #tetrahedra base triangle length
    h = 1.5 #tegrahedra height
    l = 2.0 #legs lengths
    θ = 2π/3 #triangle angle
    d = 1.0
    α = π/3
    P = vcat(
        [[2a/sqrt(3)*cos(i*θ),2a/sqrt(3)*sin(i*θ),0.0] for i = -1:1],
        [[sin(α).*l*cos(i*θ+π/3),sin(α).*l*sin(i*θ+π/3),-d-cos(α)*l] for i = -1:1],
        [[0.0,0.0,-h],[0.0,0.0,-d]]
    ) .|> SVector{3}
    
    # scatter(P)
    legs = [
        make_3d_bar(
            i,
            P[end],
            P[i+3],; ci = Int[]
        )
        for i = 1:3
    ]

    base = make_3d_tri(
        4,
        P[[7,1,2,3]],
        zero(P[7]),
        RotX(0.0),P[7];
        contactable = true,
        visible = true,
        ci = collect(1:12),
        cstr_idx = Int[],
    )
    
    # # #
    nb = 4
    sharing_elas = ElasticArray{Int}(undef, nb, 0)
    # left
    for i = 2:3
        for j = 1:3
            row = zeros(Int,nb)
            row[i] = j
            row[1] = j
            append!(sharing_elas,row)
        end
    end
    sharing = Matrix(sharing_elas')
    display(sharing)
    rbs = TypeSortedCollection(vcat(legs,[base]))
    numbered = RB.number(rbs)
    indexed = RB.index(rbs,sharing)
    
    # 
    cnt_matrix_elas = ElasticArray{Int}(undef, nb, 0)
    cnt_circular = CircularArray([2,3,4])
    # left
    for i = 1:3		
        row = zeros(Int,nb)
        row[i] = 2
        row[nb] = -1
        append!(cnt_matrix_elas,row)
        for j = 0:1		
            row = zeros(Int,nb)
            row[i] = 2
            row[nb] = -cnt_circular[i+j]
            append!(cnt_matrix_elas,row)
        end		
        row = zeros(Int,nb)
        row[i] = 2
        row[cnt_circular[i+1]-1] = -2
        append!(cnt_matrix_elas,row)
    end
    for i = 2:4		
        row = zeros(Int,nb)
        row[1] = 1
        row[nb] = -i
        append!(cnt_matrix_elas,row)
    end	
    row = zeros(Int,nb)
    row[1] = 1
    row[nb] = -1
    append!(cnt_matrix_elas,row)

    display(cnt_matrix_elas)
    cnt_matrix = Matrix(cnt_matrix_elas')
    ncables = size(cnt_matrix,1)
    hncables = ncables
    if k isa Nothing
        cables = [RB.DistanceSpringDamper3D(0.0,100.0,0.0;slack=false) for i = 1:ncables]
    else
        cables = [RB.DistanceSpringDamper3D(0.0,k[i],0.0;slack=false) for i = 1:ncables]
    end
    acs = [
        RB.RegisterActuator(i,[i],zeros(1))
        for i = 1:hncables
    ]
    hub = (actuators = acs,)
    apparatuses = (cables = cables,)
    connected = RB.connect(rbs,cnt_matrix)
    # #
    #
    # cst1 = RB.PinJoint(RB.Hen2Egg(RB.Signifier(rb1_to_3[1],2),RB.Signifier(rb4,1)))
    # cst2 = RB.PinJoint(RB.Hen2Egg(RB.Signifier(rb1_to_3[2],2),RB.Signifier(rb4,2)))
    # cst3 = RB.PinJoint(RB.Hen2Egg(RB.Signifier(rb1_to_3[3],2),RB.Signifier(rb4,3)))
    # jointedmembers = RB.join((cst1,cst2,cst3),indexed)
    #

    cnt = RB.Connectivity(numbered,indexed,@eponymtuple(connected,))

    st = RB.Structure(rbs,apparatuses,cnt)
    bot = RB.Robot(st,hub)
end