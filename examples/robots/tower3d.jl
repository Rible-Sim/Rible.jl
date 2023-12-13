function tower3d(;
    k=500.0,k1=1000.0,c=0.0,
    d = 0.1*√2/2, r2 = 0.11,		
    α = π/6,
    ijkl=1
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
           [           0,           0,  h],
        r.*[           1,           0,  0],
        r.*[cos( deg120), sin( deg120), 0],
        r.*[cos(-deg120), sin(-deg120), 0]
    ])
    for r in [r1,r2,r1,r1,r1,r1]
]
for i = 4:6
    nodess[i] .-= Ref(nodess[i][1])
end
ro_by_rbid = [
    SVector(0.0,0.0, 0),
    SVector(0.0,0.0, d),
    SVector(0.0,0.0,2d),
    SVector(0.0,0,2d+1.5h),
    SVector(0.0,0,2d+1.5h),
    SVector(0.0,-0.5h*sin(α),2d+1.5h+0.5h*cos(α))
]
R_by_rbid = [
    SMatrix(RotX(0.0)),
    SMatrix(RotZ(θ)),
    SMatrix(RotZ(2θ)),
    SMatrix(RotZ(2θ)),
    SMatrix(RotZ(2θ)*RotX(π)*RotY(α)),
    SMatrix(RotZ(2θ)*RotX(π)*RotY(α)),
]
rirjrkrl_by_rbid = [
    Ref(ro_by_rbid[i]) .+ Ref(R_by_rbid[i]).*nodess[i] for i = 1:6
]
# @show rirjrkrl_by_rbid[1]
cycle3 = [2,3,4,2]
# @show rirjrkrl_by_rbid[1]
rb1_to_3 = [
      make_3d_bar(i,rirjrkrl_by_rbid[1][cycle3[i  ]],
                    rirjrkrl_by_rbid[2][cycle3[i+1]]; ci = [1,2,3]) for i = 1:3
]

# rb4 = make_3d_tri(4,loci,ro_by_rbid[1],R_by_rbid[1],rirjrkrl_by_rbid[1][1])
# rb4 = make_3d_tri(4,loci,ro_by_rbid[1],R_by_rbid[1],rirjrkrl_by_rbid[1][1:2]...)
# rb4 = make_3d_tri(4,loci,ro_by_rbid[1],R_by_rbid[1],rirjrkrl_by_rbid[1][1:3]...)
# rb4 = make_3d_tri(4,loci[2:4],ro_by_rbid[2],R_by_rbid[2],rirjrkrl_by_rbid[2][2:4]...)
rb4_to_6 = [
    make_3d_bar(i+3,rirjrkrl_by_rbid[2][cycle3[i  ]],
                    rirjrkrl_by_rbid[3][cycle3[i+1]]; ci = Int[]) for i = 1:3
]

rb7  = make_3d_tri( 7,nodess[3],ro_by_rbid[3],R_by_rbid[3],rirjrkrl_by_rbid[3][1:4]...)

rb8  = make_3d_tri( 8,nodess[4],ro_by_rbid[4],R_by_rbid[4],rirjrkrl_by_rbid[4][1])

rb9  = make_3d_tri( 9,nodess[5],ro_by_rbid[5],R_by_rbid[5],rirjrkrl_by_rbid[5][1])

rb10 = make_3d_tri(10,nodess[6],ro_by_rbid[6],R_by_rbid[6],rirjrkrl_by_rbid[6][1:ijkl]...)

rbs = TypeSortedCollection(vcat(rb1_to_3,rb4_to_6,[rb7,rb8,rb9,]))
numbered = RB.number(rbs)
matrix_sharing = [
    4 0 0 0 1 0 0 0 0 0;
    5 0 0 0 2 0 0 0 0 0;
    6 0 0 0 3 0 0 0 0 0;
    0 4 0 0 0 1 0 0 0 0;
    0 5 0 0 0 2 0 0 0 0;
    0 6 0 0 0 3 0 0 0 0;
    0 0 4 1 0 0 0 0 0 0;
    0 0 5 2 0 0 0 0 0 0;
    0 0 6 3 0 0 0 0 0 0;

    # # #
    0 0 0 0 0 4 4 0 0 0;
    0 0 0 0 0 5 5 0 0 0;
    0 0 0 0 0 6 6 0 0 0;
    0 0 0 4 0 0 7 0 0 0;
    0 0 0 5 0 0 8 0 0 0;
    0 0 0 6 0 0 9 0 0 0;
    0 0 0 0 4 0 10 0 0 0;
    0 0 0 0 5 0 11 0 0 0;
    0 0 0 0 6 0 12 0 0 0;
    # #
    0 0 0 0 0 0 0 1 1 0;
    0 0 0 0 0 0 0 2 2 0;
    0 0 0 0 0 0 0 3 3 0;
]
indexed = RB.index(rbs,matrix_sharing[begin:end,begin:end-1])
# indexed = RB.index(rbs)
# #
ndcables = 9
nocables = 6
nvcables = 3
nocables = 6
ncables = ndcables + nocables + nvcables #+ nocables
restlend = 0.05
restleno = 0.01
restlenv = 0.05
restleno = 0.01
restlens = vcat(
    fill(restlend,ndcables),
    fill(restleno,nocables),
    fill(restlenv,nvcables),
    # fill(restleno,nocables),
)
ks = vcat(
    fill(k1,ndcables),
    fill(k,nocables),
    fill(k,nvcables),
    # fill(k,nocables),
)
cs = fill(c,ncables)
pretty_table(
    SortedDict(
        [
            ("k1 for first $ndcables", k1),
            ("k", k),
            ("c", c),
        ]
    )
)
cables = [RB.DistanceSpringDamper3D(restlens[i],ks[i],cs[i];slack=true) for i = 1:ncables]
acs = [RB.ManualActuator(1,collect(1:ncables),restlens[1:ncables])]
apparatuses = (cables = cables,)
hub = (actuators = acs,)
cnt_matrix_cables = [
    # triplex 1
    1 0 0 -1  0  0  0 0 0 0;
    0 1 0  0 -1  0  0 0 0 0;
    0 0 1  0  0 -1  0 0 0 0;
    # triplex intersect
    0 0 0  1 -1  0 0 0 0 0;
    0 0 0  0  1 -1 0 0 0 0;
    0 0 0 -1  0  1 0 0 0 0;
    # triplex 2
    0 0 0  1  0  0 -2 0 0 0;
    0 0 0  0  1  0 -3 0 0 0;
    0 0 0  0  0  1 -4 0 0 0;
    # Inner
    0 0 0 0 0 0 1 -2 0 0;
    0 0 0 0 0 0 1 -3 0 0;
    0 0 0 0 0 0 1 -4 0 0;
    # Outer
    0 0 0 0 0 0 2 -2 0 0;
    0 0 0 0 0 0 3 -3 0 0;
    0 0 0 0 0 0 4 -4 0 0;
    #
    0 0 0 0 0 0 0 2 -2 0;
    0 0 0 0 0 0 0 3 -4 0;
    0 0 0 0 0 0 0 4 -3 0;
    # # Outer
    # 0 0 0 0 0 0 0 0 2 -2;
    # 0 0 0 0 0 0 0 0 3 -3;
    # 0 0 0 0 0 0 0 0 4 -4;
    # # Inner
    # 0 0 0 0 0 0 0 0 2 -1;
    # 0 0 0 0 0 0 0 0 3 -1;
    # 0 0 0 0 0 0 0 0 4 -1;
]
connected = RB.connect(rbs,cnt_matrix_cables[begin:end,begin:end-1])
#
#
# cst1 = RB.PinJoint(RB.Hen2Egg(RB.ID(rb1_to_3[1],2),RB.ID(rb4,1)))
# cst2 = RB.PinJoint(RB.Hen2Egg(RB.ID(rb1_to_3[2],2),RB.ID(rb4,2)))
# cst3 = RB.PinJoint(RB.Hen2Egg(RB.ID(rb1_to_3[3],2),RB.ID(rb4,3)))
# jointedmembers = RB.join((cst1,cst2,cst3),indexed)
#
cnt = RB.Connectivity(numbered,indexed,@eponymtuple(connected,))

st = RB.Structure(rbs,apparatuses,cnt)
bot = RB.Robot(st,hub)
end