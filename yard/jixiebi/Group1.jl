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

nq̌ = length(bot.traj.q̌[1])
temp = zeros(nq̌, nq̌)
temp1 = zeros(nq̌, nq̌)
RB.build_tangent_stiffness_matrix!(temp1, bot.structure)



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


# ----------------------------------------------------------------

s0 = 1:38
RB.update_s!(bot.structure, s0)

bot.structure.state.members[1].s[1] = 20

foreach(bot.structure.apparatuses) do appar
    if (isa(appar, RB.Apparatus{<:RB.ClusterJoint}))
        @show  appar.joint.sps[1].s
        appar.joint.sps[1].s = 20
    end
end
foreach(bot.structure.apparatuses) do appar
    if (isa(appar, RB.Apparatus{<:RB.ClusterJoint}))
        @show  appar.joint.sps[1].s
        appar.joint.sps[1].s = 20
    end
end


# 视滑动绳索为普通绳索
RB.solve!(
    RB.DynamicsProblem(bot,),
    RB.DynamicsSolver(
        RB.Zhong06()
    );
    tspan,
    dt,
    ftol=1e-7,verbose=true
)

plot_traj!(bot)


#todo 将原来的算法src/mechanics/dynamics_solvers/Zhong06_family/Zhong06_sliding_cable_FB.jl适配到Zhong06_constant_mass_cluster_cables.jl的代码中



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