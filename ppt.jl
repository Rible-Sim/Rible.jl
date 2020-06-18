

rb1 = TRS.RigidBody(prop,state)

s1 = TRS.LinearString(k_inner,s0_inner,1,8)

j1 = TRS.BallJoint(1, 2)

connectivity = [
    [1,1], [1,2], [1,3], [1,4],
    [2,1], [2,2], [2,3], [2,4],
    [3,1], [3,2], [3,3], [3,4],
    ]
tgsys = TRS.TensegritySystem(rbs,sts,jts,connectivity)

s = 1
tab = SPARKTableau(s)
tspan = (0.0,10.0)
cache = SPARKCache(36,18,0.01,tspan,s,(A,Φ,∂T∂q̇!,F!,M!,jacs))
state = SPARKsolve!(q0,q̇0,λ0,cache,tab)
