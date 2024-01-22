# TestEnv.activate()
using Revise #jl
import Rible as RB
include(joinpath(pathof(RB),"../../test/yard/tensegrity.jl"))
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

# 以下未调试

prob1 = TR.SimProblem(bot,dynfuncs1)
dt = 1e-2; T = 9.0

TR.solve!(prob1,TR.FBZhong06(),
        (
            actuate! = actuate11!,
            prescribe! = nothing
            );
        dt=dt,tspan=(0.0,T),ftol=1e-7,verbose=true)


with_theme(theme_pub;
    Axis3 = (
        azimuth = -π/2,
        elevation = 0.0,
    ),
    Mesh = (
        # color = :black,
        transparency = false,
    ),
    ) do
plot_traj!(bot;
    # AxisType=LScene,
    # doslide=false,
    doslide=true,
    showinfo=true,
    axtitle=false,
    # hidezdecorations=true,
    hideydecorations=true,
    fontsize=10,
    # gridsize=(2,2),
    # attimes=[0.0,2.25,3.75,6.0],
    # atsteps=nothing,
    showground=false,
    showlabels=false,
    rigidcolor=:grey,
    xlims=(-100,3400),
    ylims=(-800,800),
    zlims=(-200,1800),
)
end