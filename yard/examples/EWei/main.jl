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
import Rible as RB
cd("examples/EWei")
includet("define.jl")
includet("../vis.jl")

bot1 = BuildTail()

bot1 = BuildTail(1; β=.9)

# rigid body mesh
RB.update!(bot1.st)
RB.update_orientations!(bot1.st)
rbs = RB.get_bodies(bot1.st)
rbs[5].state.R
rbs[5].state.cache.funcs.nmcs
plot_rigid(rbs[2];showmesh=true,showupdatemesh=false)

# linearized dynamics
ω², δq̌ = RB.undamped_eigen(bot1.st)
@show ω = [sqrt(ωi) for ωi in ω²[1:end]]

# nonlinear dynamics
function dynfuncs(bot)
    (;st) = bot
    function F!(F,q,q̇,s,t)
        RB.clear_forces!(st)
        RB.update_rigids!(st,q,q̇)
        RB.distribute_s̄!(st,s)
        RB.update_tensiles!(st)
        ## RB.apply_gravity!(st)
        RB.generate_forces!(st)
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
    #     st.tensiles.clustercables[1].segs[1].state.restlen += inner(t)
    # end
    Jac_F! = true
    @eponymtuple(F!, Jac_F!)
end
prob = RB.SimProblem(bot1,dynfuncs)

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
    st.tensiles.clustercables[1].segs[1].state.restlen += inner(t)
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
