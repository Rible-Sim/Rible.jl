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

tspan = (0.0,30.0)
dt = 1e-2

policy = RB.TimePolicy(
    #目前还不清楚要用到什么， 所以先用namedtuple打包
    @eponymtuple(
        #驱动量是且仅是时间的函数
        f = (t) -> [t < 10.0 ? 20.0t : 200.0, 20.0t, 0.0]
    )
)

# 一块由法向量定义的地面（半空间）
θ=0.0
ground_plane = RB.StaticContactSurfaces(
    [
        RB.HalfSpace([-tan(θ),0,1],zeros(3)),
    ]
);

# parameters
restitution_coefficients = 0.5
friction_coefficients = 0.1

bot = build_jixiebi(4)

RB.solve!(
    RB.DynamicsProblem(
        bot,policy,
        RB.EulerEytelwein(),
    ),
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.MonolithicApparatusSolver(
            RB.SmoothedFischerBurmeister()
        ),
    );
    tspan,
    dt,
    ftol=1e-7,verbose=true
)

#TODO 增加接触碰撞， 参考src/mechanics/dynamics_solvers/Zhong06_family/Zhong06_CCP_constant_mass.jl
#TODO 第一步： 增加接触检测
# （1） activate_frictional_contacts!， get_frictional_directions_and_positions!
# （2） 检查返回的接触法向等是否正确
#TODO 第二步： 求解碰撞量
# （1) 增加内层迭代 IPM!(Λₖ,na,nΛ,Λₖini,yₖini,(1/α0).*𝐍 .+ L,𝐡;ftol,Nmax)
# （2) 检查接触力是否满足摩擦碰撞定律

RB.solve!(
    RB.DynamicsProblem(
        bot,policy,
        ground_plane,
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
    tspan,
    dt,
    ftol=1e-7,verbose=true
)
plot_traj!(bot)
