function build_2d_tri(id,ri,rj=nothing,rk=nothing;
    α = 0.0, 
    ci = Int[], 
    cstr_idx = collect(1:3),        
    h = 0.05,
    b = 0.1,
    b1 = 0.05,
    m = 0.1271425597,
)
contactable = true
if isempty(ci)
    visible = true
else
    visible = true
end
# free_idx = collect(1:6)
b2 = b - b1
a1 = 1.0b; b1 = √(2b^2-a1^2)
a2 = 1.2b; b2 = √(2b^2-a2^2)
mass_center  = SVector{2}([ (b1+b)/3, h/3])
r̄p1 = SVector{2}([      0.0, 0.0])
r̄p2 = SVector{2}([       -a1, b2])
r̄p3 = SVector{2}([        a1, b2])
r̄p4 = SVector{2}([      -3a2, b2])
r̄p5 = SVector{2}([       3a2, b2])
loci_positions = [r̄p1,r̄p2,r̄p3,r̄p4,r̄p5]
@myshow loci_positions
axes = [ SVector{2}([      1.0, 0.0])]
Īgx = (b*h^3)/36
Īgy = h*b*(b^2 - b*b1 + b1^2)/36
A = b*h/2
ρ = m/A
Īgxy = -(b1^2*h^2)/24 + (b2^2*h^2)/24
dx = -2/3*(b/2-b1)
dy = -h/3
Īgxy = Īgxy + A*dx*dy
Īg = SMatrix{2,2}(
    ρ.*[
        Īgy Īgxy;
        Īgxy Īgx
    ]
)
@show norm(mass_center),m,Īg,tr(Īg)
@show mass_center,atan(mass_center[2],mass_center[1])
prop = RB.RigidBodyProperty(
    id,
    contactable,
    m,Īg,
    mass_center,
    loci_positions,
    ;visible=visible
)
ω = 0.0
ro = ri
# ṙo = [0.0,0.0001]
ṙo = zero(ro)
if id == 7
    ṙo = SVector(0.01,0.0)
# elseif id == 2
    ṙo = SVector(0.0,0.0001)
end
if rj isa Nothing
    nmcs = RB.NCF.NC1P2V(ri,ro,α)
elseif rk isa Nothing
    nmcs = RB.NCF.NC2P1V(ri,rj,ro,α)
else
    nmcs = RB.NCF.NC3P(ri,rj,rk,ro,α)
end
state = RB.RigidBodyState(prop,ro,α,ṙo,ω)
coords = RB.NonminimalCoordinates(nmcs,ci,cstr_idx)
trimesh = begin
    if id == 1
        load(RB.assetpath("零件1 - 副本.STL")) |> RB.make_patch(;scale=1/1000,color=:mediumpurple4)
    else
        load(RB.assetpath("零件1.STL")) |> RB.make_patch(;scale=1/1000,color=:slategray4)
    end
end
body = RB.RigidBody(prop,state,coords,trimesh)
end