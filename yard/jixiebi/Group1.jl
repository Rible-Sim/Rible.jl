using TestEnv
TestEnv.activate()

using Revise #jl
import Rible as RB
include(joinpath(pathof(RB),"../../yard/tensegrity.jl"))
using AbbreviatedStackTraces #jl
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true #jl
ENV["JULIA_STACKTRACE_MINIMAL"] = true #jl

figdir::String = joinpath(pathof(RB),"../../tmp")
include(joinpath(pathof(RB),"../../test/vis.jl"))
includet(joinpath(pathof(RB),"../../test/vis.jl")) #jl

include(joinpath(pathof(RB),"../../examples/robots/jixiebi.jl"))
includet(joinpath(pathof(RB),"../../examples/robots/jixiebi.jl")) #jl

tspan = (0.0,6.0)
dt = 2e-3
dt = 1e-2

policy = RB.TimePolicy(
    #目前还不清楚要用到什么， 所以先用namedtuple打包
    @eponymtuple(
        #驱动量是且仅是时间的函数
        f = (t) -> [0.0, t<15 ? t : 15, t<15 ? t : 15]
        # f = (t) -> [0.0, 0.0, 0.0]
    )
)

# 一块由法向量定义的地面（半空间）
# θ=0.0
# ground_plane = RB.StaticContactSurfaces(
#     [
#         # RB.HalfSpace([-tan(θ),0,1],[0.0,0.0,-150.0]),
#         RB.Sphere(0.400, [2.2000,0.0,-0.50+1e-3 ]) 
#     ]
# );

rigid_contacts = RB.ContactRigidBodies(
    [
        RB.RigidSphere(18)
    ]
)

# parameters
restitution_coefficients = 0.5
friction_coefficients = 0.1
bot_small = build_small_jixiebi(4)

# bot = build_jixiebi(4)

# RB.solve!(
#     RB.DynamicsProblem(
#         bot_small,policy,
#         RB.EulerEytelwein(),
#     ),
#     RB.DynamicsSolver(
#         RB.Zhong06(),
#         RB.MonolithicApparatusSolver(
#             RB.SmoothedFischerBurmeister()
#         ),
#     );
#     tspan=(0.0,0.3),
#     dt,
#     ftol=1e-7,verbose=true
# )

# bot = build_jixiebi(4)
RB.solve!(
    RB.DynamicsProblem(
        bot_small,policy,
        # ground_plane,
        rigid_contacts,
        RB.RestitutionFrictionCombined(
            RB.NewtonRestitution(),
            RB.CoulombFriction(),
        ),
        RB.EulerEytelwein(),
    ),
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        ),
        RB.MonolithicApparatusSolver(
            RB.SmoothedFischerBurmeister()
        ),
    );
    tspan=(0.0, 20.0),
    dt,
    ftol=1e-7,verbose=true,
    maxiters=50
)

# RB.solve!(
#     RB.DynamicsProblem(
#         bot, policy,
#         ground_plane,
#         RB.RestitutionFrictionCombined(
#             RB.NewtonRestitution(),
#             RB.CoulombFriction(),
#         ),
#     ),
#     RB.DynamicsSolver(
#         RB.Zhong06(),
#         RB.InnerLayerContactSolver(
#             RB.InteriorPointMethod()
#         ),
#     );
#     tspan = (0.0, 3.19),
#     dt,
#     ftol=1e-7, verbose=true
# )

plot_traj!(bot_small;
    # ground=ground_plane,
    showground  = false,
    xlims=(-1.0, 4.0),
    ylims = (-1.,1.),
    zlims = (-2.,1.),
    showlabels = false,
    showpoints = false,
)
