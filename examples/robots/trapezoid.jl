
function trapezoid(id,ri,ro=ri,
    α = 0.0, 
    ci = Int[],
    cstr_idx = collect(1:3)

)
movable = true
if isempty(ci)
    constrained = false
else
    constrained = true
end
# free_idx = collect(1:6)
b2 = b - b1
r̄p1 = SVector{2}([ -r, 0.0])
r̄p2 = SVector{2}([  r, 0.0])
r̄p3 = SVector{2}([ -b/2,  h])
r̄p4 = SVector{2}([  b/2,  h])
loci = [r̄p1,r̄p2,r̄p3,r̄p4]
ρ = 1.0
m1 = ρ*2r
m2 = ρ*h
m3 = ρ*b
m = m1+m2+m3
r̄gy = (m2*h/2+m3*h)/m
mass_locus  = SVector{2}([  0, r̄gy])
I1_c = m1*(2r)^2/12
I2_c = m2*h^2/12
I3_c = m3*b^2/12
I1_g = I1_c + m1*r̄gy^2
I2_g = I2_c + m2*(h/2-r̄gy)^2
I3_g = I3_c + m3*(h-r̄gy)^2
I_g = I1_g + I2_g + I3_g
Īg = SMatrix{2,2}(
    [
        0.99I_g 0.0;
        0.0 0.01I_g
    ]
)

prop = RB.RigidBodyProperty(id,movable,m,Īg,
            mass_locus,loci;constrained=constrained
            )
ω = 0.0
# ṙo = [0.0,0.0001]
ṙo = zero(ro)
nmcs = RB.NCF.NC1P2V(ri,ro,α)
state = RB.RigidBodyState(prop,nmcs,ro,α,ṙo,ω,ci,cstr_idx)
body = RB.RigidBody(prop,state)

end