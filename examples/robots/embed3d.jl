
function embed3d(;
        κ = nothing,
        r = 1.0,
        r1 = 1.0,
        b = 2.0,
        m = 4,
        α = 2π/m,
        θ = 1.25α, 
        n = 1,
        outer = false,
        isprism = false,
    )
    # hh = 0.5
    @assert θ>α
    c = sqrt(
        r^2 + r1^2 -2*r*r1*cos(θ)
    )
    h² = b^2-c^2
    @assert h² > 0
    h = sqrt(b^2-c^2)
    @show r,r1
    bps = [
        [
            [r*cos(i*α),r*sin(i*α),j*2h]
            for i = 1:m
        ] .|> SVector{3} |> CircularArray
        for j = 0:n
    ]
    midps = [
        [
            [r1*cos(i*α+θ),r1*sin(i*α+θ),j*2h-h]
            for i = 1:m
        ] .|> SVector{3} |> CircularArray
        for j = 1:n
    ]
    # fig = Figure()
    # ax = Axis3(fig[1,1])
    # scatter!(ax,reduce(vcat,bps))
    # scatter!(ax,reduce(vcat,midps))
    # fig
    bars = [
        vcat(
            [
                make_3d_bar(
                    2m*(j-1)+i,
                    bps[j][i],
                    midps[j][i]; ci = ifelse(j==1,[1,2,3],Int[])
                ) for i = 1:m
            ],
            [
                make_3d_bar(
                    2m*(j-1)+m+i,
                    midps[j][i-1],
                    bps[j+1][i-1];
                ) for i = 1:m
            ]
        )
        for j = 1:n
    ]
    plates = [
        begin
            if j == 1				
                ci = collect(1:12)
                cstr_idx = Int[]
            else
                ci = Int[]
                cstr_idx = collect(1:6)
            end
            id = (2n)*m+j
            ro = SVector(0.0,0.0,(j-1)*2h)
            ri = ro
            R = RotZ(0.0)
            nodes_raw = bps[1] |> Array
            nodes_ext = [
                SVector(3r̄p[1],3r̄p[2],0.0)
                for r̄p in nodes_raw
            ]
            loci = vcat(nodes_raw,nodes_ext)
            make_3d_plate(
                id,
                loci,
                ro,
                R,
                ri;
                radius=3r,
                contactable = true,
                m,
                height=1e-2,
                ci,
                cstr_idx,
                loadmesh = false,
                meshvisible = !isprism,
            )
        end
        for j = 1:n+1
    ]
    nb = sum(length.(bars)) + length(plates)
    bodies = vcat(reduce(vcat,bars),plates) |> TypeSortedCollection

    sharing_elas = ElasticArray{Int}(undef, nb, 0)
    cm = CircularArray(collect(1:m))
    for j = 1:n
        is = 2m*(j-1)
        for i = 1:m
            for k = 1:1
                row = zeros(Int,nb)
                row[is+cm[i]  ]   = k+3
                row[is+m+cm[i+1]] = k
                append!(sharing_elas,row)		
            end
        end
    end	
    for j = 2:n
        is = 2m*(j-2)+m
        for i = 1:m
            for k = 1:1
                row = zeros(Int,nb)
                row[is+cm[i]  ]   = k+3
                row[is+m+cm[i+2]] = k
                append!(sharing_elas,row)		
            end
        end
    end
    sharing_matrix = Matrix(sharing_elas')

    connecting_elas = ElasticArray{Int}(undef, nb, 0)
    for j = 1:n
        is = 2m*(j-1)
        # lower cross
        for i = 1:m
            row = zeros(Int,nb)
            row[is+cm[i  ]] =  1
            row[is+cm[i-1]] = -2
            append!(connecting_elas,row)
        end
        # additional lower cross
        # for i = 1:m
        # 	row = zeros(Int,nb)
        # 	row[is+cm[i  ]] =  1
        # 	row[is+cm[i-2]] = -2
        # 	append!(connecting_elas,row)
        # end
        # upper cross
        for i = 1:m
            row = zeros(Int,nb)
            row[is+cm[i+1]] = -2
            row[is+m+cm[i]] =  2
            append!(connecting_elas,row)
        end
        # addtional upper cross
        # for i = 1:m
        # 	row = zeros(Int,nb)
        # 	row[is+cm[i+1]] = -2
        # 	row[is+m+cm[i+1]] =  2
        # 	append!(connecting_elas,row)
        # end
        # mid
        for i = 1:m
            row = zeros(Int,nb)
            row[is+cm[i  ]] =  2
            row[is+cm[i+1]] = -2
            append!(connecting_elas,row)
        end
        # # upper
        # for i = 1:m
        # 	row = zeros(Int,nb)
        # 	row[is+m+cm[i  ]] =  2
        # 	row[is+m+cm[i+1]] = -2
        # 	append!(connecting_elas,row)
        # end
    end
    ncables_prism = size(connecting_elas,2)
    # outer
    ncables_outer = 0
    for j = 1:n
        if outer
            is = (2n)*m
            for i = 1:m
                row = zeros(Int,nb)
                row[is+j] =  m+i
                row[is+j+1] = -i-m
                append!(connecting_elas,row)
            end
            ncables_outer = n*m
        end
    end
    if isprism
        for j = 1:n
            is = 0
            for i = 1:m
                row = zeros(Int,nb)
                row[is+cm[i]] =    1
                row[is+cm[i+1]] = -1
                append!(connecting_elas,row)
                row = zeros(Int,nb)
                row[is+m+cm[i]] =    2
                row[is+m+cm[i+1]] = -2
                append!(connecting_elas,row)
            end
        end
        ncables_prism += (n+1)*m
    end
    connecting_matrix = Matrix(connecting_elas')
    # display(connecting_matrix)
    ncables = size(connecting_matrix,1)
    # @assert ncables == ncables_prism + ncables_outer

    mat_cable = filter(
        row->row.name == "Nylon 66",
        RB.material_properties
    )[1]
    diameter = 1e-3Unitful.m
    cable_length = 0.1Unitful.m
    κ = (mat_cable.modulus_elas)*π*(diameter/2)^2/cable_length
    # @show κ
    @show uconvert(Unitful.N/Unitful.m,κ),ustrip(Unitful.N/Unitful.m,κ)
    spring_dampers_prism = [RB.DistanceSpringDamper3D(0.0,   ustrip(Unitful.N/Unitful.m,κ),0.0;slack=true) for i = 1:ncables_prism]
    spring_dampers_outer = [
        RB.DistanceSpringDamper3D([0.95,0.85,0.75][((i-ncables_prism) % 3)+1]*0.2,ustrip(Unitful.N/Unitful.m,κ),0.0;slack=true) 
        for i = ncables_prism+1:ncables
    ]
    @show 2h
    spring_dampers = vcat(
        spring_dampers_prism,
        spring_dampers_outer
    )

    pin_joints = [
        RB.PinJoint(i+m*(j-1),RB.Hen2Egg(RB.Signifier(bars[j][m+cm[i]],2),RB.Signifier(plates[j+1],cm[i-1])))
        for j = 1:n for i = 1:m
    ]

    cable_joints = RB.connect(bodies,spring_dampers;connecting_matrix,istart=n*m)

    apparatuses = TypeSortedCollection(
        vcat(pin_joints,cable_joints)
    )
    # display(sharing_matrix)
    indexed = RB.index(bodies,apparatuses;sharing_matrix)
    
    numbered = RB.number(bodies,apparatuses)

    cnt = RB.Connectivity(numbered,indexed)

    st = RB.Structure(bodies,apparatuses,cnt)

    actuators = [
        RB.RegisterActuator(
            1,
            cable_joints[1:ncables_prism],
            RB.ManualOperator(),
            (values=zeros(ncables_prism),matrix=I(ncables_prism))
        ),
        RB.RegisterActuator(
            2,
            cable_joints[ncables_prism+1:ncables],
            RB.ManualOperator(),
            (values=zeros(ncables_outer),matrix=I(ncables_outer))
        ),
    ]

    gauges = Int[]
    hub = RB.ControlHub(
        st,
        gauges,
        actuators,
        RB.Coalition(st,gauges,actuators)
    )
    bot = RB.Robot(st,hub)
end