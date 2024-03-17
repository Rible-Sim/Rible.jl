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

dt = 1e-2

policy = RB.TimePolicy(
    #目前还不清楚要用到什么， 所以先用namedtuple打包
    @eponymtuple(
        #驱动量是且仅是时间的函数
        f = (t) -> reduce(vcat, [[0.0, 2t,  2t] for _ in 1:4])
        ## f = (t) -> reduce(vcat, [[0.0, 2t, 2t, 0.0, 0.5t, 0.5t] for _ in 1:2])
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
        RB.RigidSphere(53)
    ]
)

# parameters
restitution_coefficients = 0.5
friction_coefficients = 0.1
bot_4 = build_4_jixiebis(3)
# bot_small = build_small_jixiebi(3)

# bot = build_jixiebi(4)

# RB.solve!(
#     RB.DynamicsProblem(
#         bot_4,policy,
#         RB.EulerEytelwein(),
#     ),
#     RB.DynamicsSolver(
#         RB.Zhong06(),
#         RB.MonolithicApparatusSolver(
#             RB.SmoothedFischerBurmeister()
#         ),
#     );
#     tspan=(0.0,1.0),
#     dt,
#     ftol=1e-7,verbose=true
# )

# bot = build_jixiebi(4)
RB.solve!(
    RB.DynamicsProblem(
        bot_4,policy,
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
    tspan=(0.0, 2.2),
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



scale = 1e3
m0_yellow = load("yard/jixiebi/m0_yellow.STL") |> RB.make_patch(;rot=RotY(pi/2),trans=[-0.50, 0.0, 0.8],scale=1/scale,color=:yellow)
m0_blue1 = load("yard/jixiebi/m0_blue.STL") |> RB.make_patch(;rot=RotY(pi/2), trans=[-0.5,0.23,0.37], scale=1/scale,color=:blue)
m0_blue2 = load("yard/jixiebi/m0_blue.STL") |> RB.make_patch(;rot=RotY(pi/2), trans=[0.5,0.23,0.37], scale=1/scale,color=:blue)
m0_blue3 = load("yard/jixiebi/m0_blue.STL") |> RB.make_patch(;rot=RotY(pi/2)*RotX(-pi/2), trans=[-0.23,0.5,0.37], scale=1/scale,color=:blue)
m0_blue4 = load("yard/jixiebi/m0_blue.STL") |> RB.make_patch(;rot=RotY(pi/2)*RotX(-pi/2), trans=[-0.23,-0.5,0.37], scale=1/scale,color=:blue)
sat_mesh = GB.merge([m0_yellow,m0_blue1,m0_blue2,m0_blue3,m0_blue4])

plot_traj!(bot_4;
    ## ground=ground_plane,
    showground  = false,
    ## AxisType = Axis3,
    show_axis = false,
    xlims=(-1.0, 4.0),
    ylims = (-1.,1.),
    zlims = (-2.,1.),
    showlabels = false,
    showpoints = false,
    showinfo = false,
    sup! = (ax,tgob,subgrid_idx) -> begin
        mesh!(ax,sat_mesh)
        ## hidedecorations!(ax)
        ## hidespines!(ax)
    end
)

jldsave("bot_4_two.jld2"; bot_4)
jldsave("bot_4_four.jld2"; bot_4)