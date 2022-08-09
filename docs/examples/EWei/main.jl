using LinearAlgebra
using SparseArrays
using StaticArrays
using BenchmarkTools
using FileIO, MeshIO
using Printf
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using Unitful
using NLsolve
using Revise
using StructArrays
using EponymTuples
using LaTeXStrings
using StaticArrays
using GeometryBasics
using Rotations
using CoordinateTransformations
using Meshing
using TypeSortedCollections
using Match
import TensegrityRobots as TR
cd("examples/EWei")
includet("define.jl")
includet("../vis.jl")

bot1 = BuildTail()

bot1 = BuildTail(1; β=.9)

# rigid body mesh
TR.update!(bot1.tg)
TR.update_orientations!(bot1.tg)
rbs = TR.get_rigidbodies(bot1.tg)
rbs[5].state.R
rbs[5].state.cache.funcs.lncs
plot_rigid(rbs[2];showmesh=true,showupdatemesh=false)

# linearized dynamics
ω², δq̌ = TR.undamped_eigen(bot1.tg)
@show ω = [sqrt(ωi) for ωi in ω²[1:end]]

# nonlinear dynamics
function dynfuncs(bot)
    (;tg) = bot
    function F!(F,q,q̇,s,t)
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.distribute_s̄!(tg,s)
        TR.update_tensiles!(tg)
        ## TR.apply_gravity!(tg)
        TR.generate_forces!(tg)
        TR.get_force!(F,tg)
        ## F .= 0
    end
    # function apply_acu!(tg, t; dt=1e-2)
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
    #     tg.tensiles.clustercables[1].segs[1].state.restlen += inner(t)
    # end
    Jac_F! = true
    @eponymtuple(F!, Jac_F!)
end
prob = TR.SimProblem(bot1,dynfuncs)

function actuate!(tg, t; dt=1e-2)
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
    tg.tensiles.clustercables[1].segs[1].state.restlen += inner(t)
end

TR.solve!(prob,TR.FBZhong06(),
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
