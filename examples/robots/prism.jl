function prisms(;
    κ = nothing,
    r = 1.0,
    r1 = 1.0,
    b = 2.0,
    m = 4,
    α = 2π/m,
    h = nothing,
    n = 1,
    θ = 1.25α,
    hasplate = true,
)
if h isa Nothing
    @assert θ>α
    c = sqrt(
        r^2 + r1^2 -2*r*r1*cos(θ)
    )
    h² = b^2-c^2
    @assert h² > 0
    h = sqrt(h²)
end
@show r,r1
d, l = divrem(n,2)
@show θ |> rad2deg
bps = [
    begin
        if j |> iseven
            ps = [
                [r*cos(i*α),r*sin(i*α),j*h]
                for i = 1:m
            ]
        else
            ps = [
                [r1*cos(i*α+θ),r1*sin(i*α+θ),j*h]
                for i = 1:m
            ]
        end
        ps .|> SVector{3} |> CircularArray
    end
    for j = 0:n
]
# fig = Figure()
# ax = Axis3(fig[1,1])
# scatter!(ax,reduce(vcat,bps))
# scatter!(ax,reduce(vcat,midps))
# fig	
# @myshow length(bps), length(midps)
bars = [
    begin
        if j |> iseven
            [
                make_3d_bar(
                    m*(j-1)+i,
                    bps[j][i],
                    bps[j+1][i]; 
                    #ci = ifelse(j==1,[1,2,3],Int[]),
                    radius_ratio = 1/60,
                ) for i = 1:m
            ]
        else
            [
                make_3d_bar(
                    m*(j-1)+i,
                    bps[j][i-1],
                    bps[j+1][i-1];
                    radius_ratio = 1/60,
                ) for i = 1:m
            ]
        end
    end
    for j = 1:n
]

nbars = length.(bars) |> sum
nb = nbars + 1

ro = SVector(0.0,0.0,0.0)
nodes_raw = bps[1] |> Array
loci = vcat(nodes_raw,)
# nodes_ext = [
# 	SVector(3r̄p[1],3r̄p[2],0.0)
# 	for r̄p in nodes_raw
# ]
# loci = vcat(nodes_raw,nodes_ext)
if hasplate
    plate = make_3d_plate(
        nb,
        loci,
        ro,
        RotZ(0.0),
        ro;
        radius=r,
        contactable = true,
        m=3,
        height=1e-2,
        ci = collect(1:12),
        cstr_idx = Int[],
        loadmesh = false,
        meshvisible = true,
    )
    rbs = vcat(
        reduce(vcat,bars),
        plate
    ) |> TypeSortedCollection
else
    rbs = vcat(
        reduce(vcat,bars),
    ) |> TypeSortedCollection
end
numberedpoints = RB.number(rbs)
cm = CircularArray(collect(1:m))
# sharing_elas = ElasticArray{Int}(undef, nb, 0)
# for j = 1:n
#     is = 2m*(j-1)
#     for i = 1:m
#         for k = 1:p
#             row = zeros(Int,nb)
#             row[is+cm[i]  ]   = k+3
#             row[is+m+cm[i+1]] = k
#             append!(sharing_elas,row)		
#         end
#     end
# end	
# for j = 2:n
#     is = 2m*(j-2)+m
#     for i = 1:m
#         for k = 1:p
#             row = zeros(Int,nb)
#             row[is+cm[i]  ]   = k+3
#             row[is+m+cm[i+2]] = k
#             append!(sharing_elas,row)		
#         end
#     end
# end
# sharing = Matrix(sharing_elas')
# display(sharing)
# indexedcoords = RB.index(rbs,sharing)
indexedcoords = RB.index(rbs,)
    
@myshow nb
connecting_elas = ElasticArray{Int}(undef, nb, 0)
for j = 0:n-1
    is = m*j
    
    # # lower
    # for i = 1:m
    #     row = zeros(Int,nb)
    #     row[is+cm[i  ]] =  1
    #     row[is+cm[i+1]] = -1
    #     append!(connecting_elas,row)
    # end

    # cross
    if j |> iseven
        for i = 1:m
            row = zeros(Int,nb)
            row[is+cm[i  ]] =  1
            row[is+cm[i-1]] = -2
            append!(connecting_elas,row)
        end
        # additional cross
        # for i = 1:m
        # 	row = zeros(Int,nb)
        # 	row[is+cm[i  ]] =  1
        # 	row[is+cm[i-2]] = -2
        # 	append!(connecting_elas,row)
        # end
    else
        for i = 1:m
            row = zeros(Int,nb)
            row[is+cm[i+1]] = -2
            row[is+cm[i  ]] =  1
            append!(connecting_elas,row)
        end
        # addtional cross
        # for i = 1:m
        # 	row = zeros(Int,nb)
        # 	row[is+cm[i+1]] = -2
        # 	row[is+m+cm[i+1]] =  2
        # 	append!(connecting_elas,row)
        # end
    end

    # upper 
    for i = 1:m
        row = zeros(Int,nb)
        row[is+cm[i  ]] =  2
        row[is+cm[i+1]] = -2
        append!(connecting_elas,row)
    end
end
ncables_prism = size(connecting_elas,2)
connecting = Matrix(connecting_elas')
# display(connecting)
connected = RB.connect(rbs,connecting)

ncables = size(connecting,1)
# @assert ncables == ncables_prism + ncables_outer
@show ncables
mat_cable = filter(
    row->row.name == "Nylon 66",
    material_properties
)[1]
diameter = 1e-3Unitful.m
cable_length = 0.1Unitful.m
κ = (mat_cable.modulus_elas)*π*(diameter/2)^2/cable_length
# @show κ
@show uconvert(Unitful.N/Unitful.m,κ),ustrip(Unitful.N/Unitful.m,κ)
cables_prism = [RB.DistanceSpringDamper3D(0.0,   ustrip(Unitful.N/Unitful.m,κ),0.0;slack=true) for i = 1:ncables_prism]
cables = cables_prism
acs = [
    RB.ManualActuator(
        1,
        collect(1:ncables_prism),
        zeros(ncables_prism)
    ),
]
apparatuses = (cables = cables,)
hub = (actuators = acs,)


csts_bar2bar = [
    begin 
        hen = bars[j+1][cm[i+2]]
        @show hen.prop.id
        egg = bars[j][cm[i]]  
        @show egg.prop.id
        RB.PinJoint(
            i+m*(j-1),
            RB.Hen2Egg(
                i+m*(j-1),
                RB.ID(hen,1,1),
                RB.ID(egg,2,1)
            )
        )
    end
    for j = 1:d for i = 1:m
]

if hasplate
    csts_bar2plate = [
        RB.PinJoint(
            m*d+i,
            RB.Hen2Egg(
                m*d+i,
                RB.ID(plate,i,1),
                RB.ID(bars[1][cm[i]],1,1),
            )
        )
        for i = 1:m
    ]
    csts = vcat(
        csts_bar2plate,
        csts_bar2bar,
    )
else
    csts = csts_bar2bar
end

jointedmembers = RB.join(csts,indexedcoords)

cnt = RB.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,),jointedmembers)

st = RB.Structure(rbs,apparatuses,cnt)
bot = RB.Robot(st,hub)
end