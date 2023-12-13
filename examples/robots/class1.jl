
function class1(;
    k=500.0,c=0.0,
    d = 0.1*√2/2, r2 = 0.11,
    ijkl=1,
    R1 = RotX(0),
    )
r1 = 0.1
b = 0.22
h = 0.1*√2/2
γ = acos((d^2+r1^2+r2^2-b^2)/(2r1*r2))
θ =  γ - 2π/3
@show rad2deg.([γ,θ])
deg120 = deg2rad(120)
nodess = [
    SVector{3}.([
           [           0,           0, 1.5h],
        r.*[           1,           0,  0],
        r.*[cos( deg120), sin( deg120), 0],
        r.*[cos(-deg120), sin(-deg120), 0]
    ])
    for r in [r1,r1]
]
loci = nodess[1]

ro_by_rbid = [
    SVector(0.0,0.0,0.0),
    SVector(0.0,0.0,1.0h)
]
R_by_rbid = [
    SMatrix(RotZ(0.0)),
    SMatrix(R1),
]
rirjrkrl_by_rbid = [
    Ref(ro_by_rbid[i]) .+ Ref(R_by_rbid[i]).*nodess[i] for i = 1:2
]

cycle3 = [2,3,4,2]
rb1 = make_3d_tri(1,loci,ro_by_rbid[1],R_by_rbid[1],rirjrkrl_by_rbid[1][1];			
    visible = true,
    ci = collect(1:12),
    cstr_idx = Int[]
)

rb2 = make_3d_tri(2,loci,ro_by_rbid[2],R_by_rbid[2],rirjrkrl_by_rbid[2][1:ijkl]...)

rbs = TypeSortedCollection([rb1,rb2])
numbered = RB.number(rbs)
# matrix_sharing = [
# ]
# indexed = RB.index(rbs,matrix_sharing)
indexed = RB.index(rbs)
# #
ncables = 6
restlen = 0.01
restlens = fill(restlen,ncables)
ks = fill(k,ncables)
cs = fill(c,ncables)
cables = [RB.DistanceSpringDamper3D(restlens[i],ks[i],cs[i];slack=true) for i = 1:ncables]
acs = [RB.ManualActuator(1,collect(1:ncables),restlens[1:ncables])]
apparatuses = (cables = cables,)
hub = (actuators = acs,)
cnt_matrix_cables = [
    # Outer
    2 -2;
    3 -3;
    4 -4;
    # Inner
    -1 2;
    -1 3;
    -1 4;
]
connected = RB.connect(rbs,cnt_matrix_cables)
#
cnt = RB.Connectivity(numbered,indexed,@eponymtuple(connected,))

st = RB.Structure(rbs,apparatuses,cnt)
RB.Robot(st,hub)
end
