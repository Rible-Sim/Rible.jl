
function twotre3d(;
    κ = nothing,
    r = 1.0,
    r1 = 1.0,
    b = 2.0,
    m = 4,
    α = 2π/m,
    θ = 1.25α, 
    n = 1,
    outer = false,
    loadmesh = true,
    isspine = false,
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
    vcat(
        [
            [0,0,h]
        ],
        [
            [3r*cos(i*α),3r*sin(i*α),0.0]
            for i = 1:m
        ]
    ) .|> SVector{3}
]
# fig = Figure()
# ax = Axis3(fig[1,1])
# scatter!(ax,reduce(vcat,bps))
# scatter!(ax,reduce(vcat,midps))
# fig
plates = [
    begin
        if j == 1				
            ci = collect(1:12)
            cstr_idx = Int[]
        else
            ci = Int[]
            cstr_idx = collect(1:6)
        end
        id = j
        if isspine
            ro = SVector(0.0,0.0,0.9(j-1)*h)
            R = RotX(0.0)
        else
            ro = SVector(0.0,0.0,(j-1)*2h)
            R = RotX(((id-1)*π))
        end
        ri = ro
        nodes_raw = bps[1] |> Array
        # nodes_ext = [
        # 	SVector(3r̄p[1],3r̄p[2],0.0)
        # 	for r̄p in nodes_raw
        # ]
        loci = vcat(nodes_raw,)
        make_3d_tri(
            id,
            loci,ro,
            R,ri,;
            # movable = true,
            # constrained = false,
            # ci = Int[],
            # cstr_idx = collect(1:6),
            color = :slategrey,
            loadmesh,
        )
    end
    for j = 1:n+1
]
nb =  length(plates)
rbs = plates |> TypeSortedCollection
numberedpoints = RB.number(rbs)
# indexedcoords = RB.index(rbs,sharing)
indexedcoords = RB.index(rbs,)

connecting_elas = ElasticArray{Int}(undef, nb, 0)
# outer
if isspine
    ncables = 2m
    for i = 2:m+1
        row = zeros(Int,nb)
        row[1] =  i
        row[2] = -i
        append!(connecting_elas,row)
    end
    for i = 2:m+1
        row = zeros(Int,nb)
        row[1] =  1
        row[2] = -i
        append!(connecting_elas,row)
    end
else
    ncables = m
    row = zeros(Int,nb)
    row[1] =   2
    row[2] = -(3)
    append!(connecting_elas,row)
    row = zeros(Int,nb)
    row[1] =   3
    row[2] = -(2)
    append!(connecting_elas,row)
    row = zeros(Int,nb)
    row[1] =   4
    row[2] = -(4)
    append!(connecting_elas,row)
end
connecting = Matrix(connecting_elas')
display(connecting)
connected = RB.connect(rbs,connecting)

ncables = size(connecting,1)
# @assert ncables == ncables_prism + ncables
κ0 = 72e9*π*(3e-3)^2/1.0
@show κ0

cables = [RB.Cable3D(i,0.8*2h,1e3,0.0;slack=true) for i = 1:ncables]

acs = [
    RB.ManualActuator(
        2,
        collect(1:ncables),
        zeros(ncables)
    ),
]
tensiles = (cables = cables,)
hub = (actuators = acs,)

# jointedmembers = RB.join(csts,indexedcoords)
cnt = RB.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,),)

st = RB.Structure(rbs,tensiles,cnt)
RB.Robot(st,hub)
end