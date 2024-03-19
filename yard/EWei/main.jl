
using Revise #jl
import Rible as RB

include(joinpath(pathof(RB),"../../yard/stability_stiffness.jl"))
using AbbreviatedStackTraces #jl
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true #jl
ENV["JULIA_STACKTRACE_MINIMAL"] = true #jl

figdir::String = ""

include(joinpath(pathof(RB),"../../test/vis.jl"))
includet(joinpath(pathof(RB),"../../test/vis.jl")) #jl

include(joinpath(pathof(RB),"../../examples/robots/EWei.jl"))
includet(joinpath(pathof(RB),"../../examples/robots/EWei.jl")) #jl


bot1 = BuildTail()

bot1 = BuildTail(1; β=.9)

# rigid body mesh
RB.update!(bot1.structure)
rbs = RB.get_bodies(bot1.structure)
rbs[5].state.origin_frame.axes
rbs[5].coords.nmcs
RB.viz(rbs[2];showmesh=true,showupdatemesh=false)

# linearized dynamics
ω², δq̌ = RB.undamped_eigen(bot1.structure)
@show ω = [sqrt(ωi) for ωi in ω²[1:end]]

plot_traj!(bot1)

# nonlinear dynamics
function dynfuncs(bot)
    (;st) = bot
    function F!(F,q,q̇,s,t)
        RB.clear_forces!(st)
        RB.update_bodies!(st,q,q̇)
        RB.distribute_s̄!(st,s)
        RB.update_apparatuses!(st)
        ## RB.apply_gravity!(st)
        RB.assemble_forces!(st)
        RB.get_force!(F,st)
        ## F .= 0
    end
    # function apply_acu!(st, t; dt=1e-2)
    #     function inner(t;dt=dt)
    #         a = 10.0
    #         if 0<t<1
    #             return -a*t*dt
    #         elseif 1<=t<5
    #             return -a*dt
    #         elseif 5<=t<6
    #             return a*t*dt - 6a*dt
    #         else
    #             return 0
    #         end
    #     end
    #     st.apparatuses.clustercables[1].segs[1].state.restlen += inner(t)
    # end
    Jac_F! = true
    @eponymtuple(F!, Jac_F!)
end
prob = RB.DynamicsProblem(bot1,dynfuncs)

function actuate!(st, t; dt=1e-2)
    function inner(t;dt=dt)
        a = 10.0
        if 0<t<1
            return -a*t*dt
        elseif 1<=t<5
            return -a*dt
        elseif 5<=t<6
            return a*t*dt - 6a*dt
        else
            return 0
        end
    end
    st.apparatuses.clustercables[1].segs[1].state.restlen += inner(t)
end

RB.solve!(prob,RB.FBZhong06(),
        (actuate! = actuate!, prescribe! = nothing);
        dt=0.01,tspan=(0.0,5.0),ftol=1e-7,verbose=true)

# try LScene
plot_traj!(bot1;
    AxisType=LScene,
    showinfo=false,
    showground=false,
    showlabels=false,
)

# try Axis3
with_theme(theme_try;
        Axis3 = (
            azimuth = 5.10553063332698,
            elevation = 0.9397963267948968
        ),
        Mesh = (
            # color = :black,
            transparency = false,
        ),
    ) do
    plot_traj!(bot1;
        showinfo=true,
        showground=false,
        showlabels=false,
        rigidcolor=:blue,
        xlims=(-100,100),
        ylims=(-100,1000),
        zlims=(-100,200)
    )
end

# publications
with_theme(theme_pub;
        Axis3 = (
            azimuth = 5.4355306333269855,
            elevation = 0.3126990816987244,
        ),
        Mesh = (
            # color = :black,
            transparency = false,
        ),
        Label = (
            font = "CMU Serif",
        ),
    ) do
    plot_traj!(bot1;
        figsize=(0.9tw,1.2tw),
        doslide=false,
        showinfo=false,
        gridsize=(2,2),
        attimes=[0,1,2,3],
        atsteps=nothing,
        showground=false,
        showlabels=false,
        rigidcolor=:white,
        xlims=(-100,100),
        ylims=(-100,1000),
        zlims=(-100,200),
        figname="EWei"
    )
end
