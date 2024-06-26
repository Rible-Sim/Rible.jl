
function spine3d(n;
    c=0.0,
    RR=RotX(0.0)
)
nbp = 4*n

a = 0.04 #m
h = 0.04 #m
θ = π/4

mass = 0.1 #kg
#inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-8 # N/m^2
inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-2
mass_locus = [0.0, 0.0, 0.0]

r̄p_x = cos(θ)*a
r̄p_y = sin(θ)*a
r̄p5 = SVector{3}([  0.0,   0.0,   0.0])
r̄p4 = SVector{3}([  0.0,  r̄p_x,  r̄p_y])
r̄p3 = SVector{3}([  0.0, -r̄p_y,  r̄p_y])
r̄p2 = SVector{3}([ r̄p_x,   0.0, -r̄p_y])
r̄p1 = SVector{3}([-r̄p_x,   0.0, -r̄p_y])
loci = [r̄p1,r̄p2,r̄p3,r̄p4,r̄p5]

contactable = ones(Bool,n)
contactable[1] = false
visible = zeros(Bool,n)
visible[1] = true

props = [RB.RigidBodyProperty(i,contactable[i],mass,
            SMatrix{3,3}(inertia),
            SVector{3}(mass_locus),
            loci;visible=visible[i]) for i = 1:n]

rs = [[0.0,0.0,i*h] for i = 0:n-1]
Rs = [Matrix(RR)^i for i = 0:n-1]
ṙs = [zeros(3) for i = 0:n-1]
ωs = [zeros(3) for i = 0:n-1]

function rigidbody(i,prop,loci,r,R,ṙ,ω)
    if i == 1
        ci = collect(1:12)
        cstr_idx = Int[]
    else
        ci = Int[]
        cstr_idx = collect(1:6)
    end

    # ri,rj,rk,rl = [r+R*r̄p for r̄p in loci]
    # nmcs = RB.NCF.NC4P(ri,rj,rk,rl,r,R,ṙ,ω)
    # nmcs = RB.NCF.NC3P1V(ri,rk,rl,r,R,ṙ,ω)
    # nmcs = RB.NCF.NC2P2V(rk,rl,r,R,ṙ,ω)
    nmcs = RB.NCF.NC1P3V(r,r,R,ṙ,ω)
    state = RB.RigidBodyState(prop,r,R,ṙ,ω,ci,cstr_idx)
    radius = norm(loci[1]-loci[5])/25
    vertmesh = merge([
        endpoints2mesh(loci[i],loci[5];radius)
        for i in 1:4
    ])
    body = RB.RigidBody(prop,state,vertmesh)
end
rbs = [rigidbody(i,props[i],loci,rs[i],Rs[i],ṙs[i],ωs[i]) for i = 1:n]
rigdibodies = TypeSortedCollection(rbs)
numbered = RB.number(rigdibodies)
indexed = RB.index(rigdibodies)

ncables = 8*(n-1)
cablelenH = 0.6h
cablelenR = 0.04
cablelens = repeat(vcat(fill(cablelenH,4),fill(cablelenR,4)),n-1)
kH = 400e1
kR = 400e1
ks = repeat(vcat(fill(kH,4),fill(kR,4)),n-1)
# c = 0.0
cs = repeat(fill(c,8),n-1)
cables = [RB.DistanceSpringDamper3D(cablelens[i],ks[i],cs[i]) for i = 1:ncables]
apparatuses = (cables=cables,)	
acs = [RB.RegisterActuator(1,collect(1:ncables),cablelens[1:ncables])]
hub = (actuators=acs,)

matrix_cnt_raw = Vector{Matrix{Int}}()
for i = 2:n
    s = zeros(Int,8, n)
    for j = 1:4
        s[j,   i-1] = j;  s[j, i] = -j
    end
    s[5, i-1] = 3;  s[5, i] = -1
    s[6, i-1] = 4;  s[6, i] = -1
    s[7, i-1] = 3;  s[7, i] = -2
    s[8, i-1] = 4;  s[8, i] = -2

    # s = zeros(Int,6, n)
    # for j = 1:2
    # 	s[j,   i-1] = j+2;  s[j, i] = -(j+2)
    # end
    # s[3, i-1] = 3;  s[3, i] = -1
    # s[4, i-1] = 4;  s[4, i] = -1
    # s[5, i-1] = 3;  s[5, i] = -2
    # s[6, i-1] = 4;  s[6, i] = -2

    push!(matrix_cnt_raw, s)
end
matrix_cnt = reduce(vcat, matrix_cnt_raw)
# display(matrix_cnt)
connected = RB.connect(rigdibodies, matrix_cnt)
tensioned = @eponymtuple(connected,)
cnt = RB.Connectivity(numbered, indexed, tensioned)
st = RB.Structure(rigdibodies,apparatuses,cnt)
bot = RB.Robot(st,hub)
end
