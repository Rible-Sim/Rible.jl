function simple(;
    c=0.0,
    m = 4,
    α = 2π/m,
    k = 100.0,
    z0 = 0.2,
    ωz = 5.0,
    mbar = 0.05,
    free = false,
)
lₛ = 15e-3
DE = 300e-3
d = (150*√2)*1e-3
r = d/2
bps = [
    [r*cos(i*α),r*sin(i*α),0.0]
    for i = 1:m
] .|> SVector{3}
ro = SVector(0,0,z0)
Rplate = RotX(π/4)
Rbar = RotY(-π/12)
p = Ref(Rbar) .* [
    [0,0, DE/2],
    [0,0,-DE/2]
] .+ Ref(ro) .|> SVector{3}

# fig,ax,plt = scatter(p)
# # ax.aspect = DataAspect()
# fig
ṙo_ref = SVector(0.0,0,0)

rbs = [
    make_3d_bar(
        1, 
        p[end-1],
        p[end];
        # ṙi = ṙo_ref, 
        # ṙj = ṙo_ref,
        m = mbar,
        mat_name = "Teak"
    ),
    make_3d_plate(
        2, 
        bps,ro,
        Rplate,ro;
        m,
        radius = r,
        # constrained = !free,	
        ci = ifelse(free,Int[],collect(1:12)),
        cstr_idx = ifelse(free,collect(1:6),Int[]),
    )
]

rigdibodies = TypeSortedCollection(rbs)
numberedpoints = RB.number(rigdibodies)
indexedcoords = RB.index(rigdibodies)
#
ncables = 2m

original_restlens = zeros(ncables)
original_restlens[  1: m] .= 0.05
original_restlens[m+1:2m] .= 0.1

ks = zeros(ncables)
ks[  1: m] .= ks[m+1:2m] .= k

cables = [
    RB.Cable3D(i, original_restlens[i], ks[i], c;slack=true) for i = 1:ncables
]

pretty_table(
    reduce(vcat,[
            [
                "No.$i" "Rest length" "$(original_restlens[i])" "Stiffness" "$(ks[i])" "Damping" "$c";
            ] 
            for i = 1:2m
        ]
    )
)
#
tensiles = (cables = cables,)
acs = [
RB.ManualActuator(
    i,
    collect(1:2m), [original_restlens[2m*(i-1)+j] for j = 1:6],
) for i = [1]
]
hub = (actuators = acs,)
# #
matrix_cnt = zeros(Int,ncables,2)
for j = 1:ncables
if j <= m
    matrix_cnt[j, 1:2] = [1, -j]
else
    matrix_cnt[j, 1:2] = [2, -(j-m)]
end
end
# display(matrix_cnt)
connected = RB.connect(rigdibodies, matrix_cnt)
tensioned = @eponymtuple(connected,)
#
cnt = RB.Connectivity(numberedpoints, indexedcoords, tensioned)
# #
st = RB.Structure(rigdibodies, tensiles, cnt, )
bot = RB.Robot(st, hub)
end