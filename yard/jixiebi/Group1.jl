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

bot = build_jixiebi(4)
plot_traj!(bot)

tspan = (0.0,1.0)
dt = 1e-2

RB.solve!(
    RB.DynamicsProblem(bot,RB.EulerEytelwein()),
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.MonolithicApparatusSolver(
            RB.SmoothedFischerBurmeister()
            # 没有使用，只用来多重派发
        )
    );
    tspan,
    dt,
    ftol=1e-7,verbose=true
)
plot_traj!(bot)
