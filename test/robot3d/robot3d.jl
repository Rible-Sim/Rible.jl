using LinearAlgebra
using StaticArrays
using Revise
using TensegrityRobot
const TR = TensegrityRobot

mass = 1.0
inertia = SMatrix{3,3}(1.0I)
r̄g = zeros(3)
r = [0.0,0.0,-1.0]
R = Matrix(1.0I,3,3)
ṙ = [0.0,0.0,0.0]
ω = [0.0,0.0,0.0]
movable = true

a = 0.5 #m
h = 1.0 #m
θ = 2π/3
l = √(a^2+h^2)
b = √3a

# Configuration
offset = [0.0, 0.0, h]
ap1 = SVector{3}([a, 0.0, 0.0] + offset)
ap2 = SVector{3}([a*cos(θ), a*sin(θ), 0.0] + offset)
ap3 = SVector{3}([a*cos(θ), -a*sin(θ), 0.0] + offset)
ap4 = SVector{3}([0.0, 0.0, -h] + offset)

prop = TR.RigidBody3DProperty(1,movable,mass,
            SMatrix{3,3}(inertia),
            SVector(r̄g...),
            [ap1,ap2,ap3,ap4])
 
state = TR.RigidBody3DState(prop,r,R,ṙ,ω,Val(:NC))

TR.RigidBody(prop,state)

mass = 1.0 #kg
#inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-8 # N/m^2
inertia = Matrix(Diagonal([45.174,45.174,25.787]))*1e-1
r̄g = [0.0, 0.0, 17.56] .* 1e-4 # m
rb1 = RigidBody(:rb1,mass = mass, inertia = inertia, r̄g = r̄g,
                     r = [0.0,0.0, -1.0], movable = false)
