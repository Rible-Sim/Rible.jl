using EnhancedGJK
import GeometryTypes: HyperRectangle, HyperSphere, Vec, Point

c1 = HyperRectangle(Vec(0., 0, 0), Vec(1., 1, 1))
c2 = HyperRectangle(Vec(3., 0, 0), Vec(1., 1, 1))
result = gjk(c1, c2)
