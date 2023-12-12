function scissor_tower(;
    θ = π/4,
    β = 0.0,#-π/4,
    k = nothing,
)
b = 1.0
a = sqrt(2)/2*b
c = cos(θ)*b
d = sin(θ)*b
n = 3
bps_raw = [zeros(2) for i = 1:2n*3]
for j = 1:n
    for i = 1:3
        is = 3(j-1)
        if j == n
            bps_raw[is+1] .= [-a/2,d*(j-1)]
            bps_raw[is+2] .= [ a/2,d*(j-1)]
            bps_raw[is+3] .= [ 0.0,d*(j-1)+a/2]
        else
            bps_raw[is+1] .= [ a/2-c,d*(j-1)]
            bps_raw[is+2] .= [-a/2+c,d*(j-1)]
            bps_raw[is+3] .= [ 0.0,d*(j-1)+a/2]
        end
    end
end
iss = n*3
R = RB.rotation_matrix(β)
x = SVector(0.0,a/2)
for j = 1:n
    for i = 1:3
        is = 3(j-1)
        if j |> isodd
            bps_raw[iss+is+1] .= bps_raw[iss].+R*(x+[-a/2,d*(j-1)])
            bps_raw[iss+is+2] .= bps_raw[iss].+R*(x+[ a/2,d*(j-1)])
            bps_raw[iss+is+3] .= bps_raw[iss].+R*(x+[ 0.0,d*(j-1)+a/2])
        else
            bps_raw[iss+is+1] .= bps_raw[iss].+R*(x+[ a/2-c,d*(j-1)])
            bps_raw[iss+is+2] .= bps_raw[iss].+R*(x+[-a/2+c,d*(j-1)])
            bps_raw[iss+is+3] .= bps_raw[iss].+R*(x+[ 0.0,d*(j-1)+a/2])
        end
    end
end
bps = SVector{2}.(bps_raw)
display(bps)
# fig,_,_=scatter(bps)
# return fig

# α_tri = 0.0
# bars = [
#     build_2d_bar(1,bps[1],bps[3];α=α_bar,ci=collect(1:2))
#     build_2d_bar(1,bps[1],bps[3];α=α_bar,ci=collect(1:2))
#     for i = 1:2
# ]

bar1 = build_2d_bar(1,bps[1],bps[5];α=θ,ci=collect(1:2))
bar2 = build_2d_bar(2,bps[2],bps[4];α=π-θ,ci=collect(1:2))
bar3 = build_2d_bar(3,bps[4],bps[8];α=θ,)
bar4 = build_2d_bar(4,bps[5],bps[7];α=π-θ,)
tri1 = build_2d_tri(5,bps[7],bps[8],bps[9];
    α = 0.0,
    b = b/sqrt(2),
    h = b/sqrt(2)/2,
    b1 = b/sqrt(2)/2
)
tri2 = build_2d_tri(6,bps[11],bps[10],bps[9];
    α = π+β,
    b = b/sqrt(2),
    h = b/sqrt(2)/2,
    b1 = b/sqrt(2)/2
)

# rb2 = build_2d_bar(2,bps[2],bps[4];α=α_bar,ci=collect(1:2))
# rb1 = build_2d_tri(1,bps[1],bps[3];α=α_bar,ci=collect(1:2))
# rb2 = build_2d_tri(2,bps[2],bps[4];α=α_bar,ci=collect(1:2))
# rb3 = build_2d_tri(3,bps[3],bps[4],bps[9];α=α_tri)
# rb3 = build_2d_tri(3,bps[3],;α=α_tri)
# rb4 = build_2d_bar(4,bps[9],bps[10];α=α_bar)
# rb5 = build_2d_tri(5,bps[5],bps[6];α=α_tri)
# rb6 = build_2d_bar(6,bps[10],bps[11];α=α_bar)
# rb7 = build_2d_tri(7,bps[7],bps[8],bps[11];α=α_tri)
# rb8 = build_2d_ground(8)
# rbs = TypeSortedCollection((rb1,rb2,rb3,rb4,rb5,rb6,rb7,rb8))
rbs = TypeSortedCollection((bar1,bar2,bar3,bar4,tri1,tri2))
numberedpoints = RB.number(rbs)
sharing = [
    3 0 0 1 0 0;
    4 0 0 2 0 0;
    0 3 1 0 0 0;
    0 4 2 0 0 0;
    0 0 0 3 1 0;
    0 0 0 4 2 0;
    0 0 3 0 3 0;
    0 0 4 0 4 0;
]
indexedcoords = RB.index(rbs,sharing)
# #
# restlen4 = ratio1*0.1*√5
# restlen1 = ratio*0.1*√2
# restlen2 = ratio*0.1
# restlen3 = ratio*0.05*√2
# restlens = [
#     restlen4,restlen4,
#     restlen1,restlen1,
#     restlen2,restlen2,
#     restlen3,restlen3,
#     restlen2,restlen2,
#     restlen3,restlen3,
# ]
# ncables = length(restlens)
# naux = 2
# ness = 10
# ks = vcat(fill(k,naux),fill(k,ness))
# cs = fill(c,ncables)
# cables = [RB.DistanceSpringDamper2D(restlens[i],ks[i],cs[i];slack) for i = 1:ncables]
# acs = [RB.ManualActuator(1,collect(1:ncables),restlens[1:ncables])]
# apparatuses = (cables = cables,)
# hub = (actuators = acs,)
cnt_matrix = [
     1 -2  0  0 0  0;
       -2  1  0  0 0  0;
     0  0  1 -2 0  0;
     0  0 -2  1 0  0;
     2 -2  0  0 0  0;
     0  0  0  0 1 -2;
     0  0  0  0 2 -1;
]
connected = RB.connect(rbs,cnt_matrix)
# #

# cst1 = RB.PinJoint(RB.Hen2Egg(RB.ID(rb2,2),RB.ID(rb3,2)))
# cst2 = RB.PinJoint(RB.Hen2Egg(RB.ID(rb3,3),RB.ID(rb4,1)))
# cst3 = RB.PinJoint(RB.Hen2Egg(RB.ID(rb4,2),RB.ID(rb5,3)))
# cst4 = RB.PinJoint(RB.Hen2Egg(RB.ID(rb5,3),RB.ID(rb6,1)))
# jointedmembers = RB.join((cst1,cst2,cst3,cst4),indexedcoords)
# # jointedmembers = RB.join((cst1,cst2,cst3),indexedcoords)

ncables = size(cnt_matrix,1)
hncables = ncables
if k isa Nothing
    cables = [RB.DistanceSpringDamper2D(0.0,1000.0,0.0;slack=false) for i = 1:ncables]
else
    cables = [RB.DistanceSpringDamper2D(0.0,k[i],0.0;slack=false) for i = 1:ncables]
end
acs = [
    RB.ManualActuator(i,[i],zeros(1))
    for i = 1:hncables
]
hub = (actuators = acs,)
apparatuses = (cables = cables,)
connected = RB.connect(rbs,cnt_matrix)


cnt = RB.Connectivity(numberedpoints,indexedcoords,@eponymtuple(connected,),)
# # cnt = RB.Connectivity(numberedpoints,indexedcoords,connections)
st = RB.Structure(rbs,apparatuses,cnt)
RB.Robot(st,hub)
end