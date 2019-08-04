using StaticArrays
using Revise
using TensegrityRobotSimulator
CD = TensegrityRobotSimulator.CollisionDetection
s = CD.Sphere(1.0)
p = CD.Plane([1.0,1.0,0.0],0.4)
hs = CD.Halfspace([0.4,0.5,0.0],0.4)

s = CD.sphere_object(1.0,[2.0,1.0,2.0])
r = 2.0
normal = [0.0,0.0,1.0]
ground = CD.halfspace_object(normal,0.0)
plane = CD.plane_object(normal,0.0)
CD.isintersect(s,plane)
CD.isintersect(s,ground)
CD.isinside(s,ground)
s1 = CD.sphere_object(1.0,[0.0,1.0,2.0])
s2 = CD.sphere_object(2.0,[2.5,1.0,2.0])
c = CD.contact(s1,s2)
