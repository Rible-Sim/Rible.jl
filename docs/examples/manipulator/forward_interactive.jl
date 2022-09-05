using LinearAlgebra
using SparseArrays
using StaticArrays
using StructArrays
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using LaTeXStrings
using RecursiveArrayTools
using Cthulhu
using XLSX, DataFrames
using Unitful
using GeometryBasics
using TypeSortedCollections
using EponymTuples
using Match
using HomotopyContinuation
using Random
using Evolutionary
using Printf
using Revise
import TensegrityRobots as TR
cd("examples/manipulator")
includet("man_define.jl")
includet("../analysis.jl")
includet("../vis.jl")
includet("man_plotting.jl")


k = 453.09; c = 100.0;
restlen = 0.174-3.9061/k
man_inv = man_ndof(8,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=true); man_inv_plot = deepcopy(man_inv)

plot_traj!(man_inv)

function get_interactive(bot_input;g0=[0.0])
    bot = deepcopy(bot_input)
    (;tg) = bot
    (;nrigids) = tg
    (;mem2num,num2ID,num2sys) = tg.connectivity.numbered
    tg_end = deepcopy(tg)
    nsegs = nrigids-1
    startsols,parameters0 = TR.get_start_system(bot,TR.DeformMode())

    u0 = parameters0.u
    d0 = parameters0.d
    g1 = g0

    indx_d = reduce(
        vcat,
        [
            [3(i-1)+1,3(i-1)+2] for i = 1:nrigids
        ]
    )
    indx_u = collect(1:length(u0))
    d1 = copy(d0)
    d1[indx_d].+=0.0001
    u1 = copy(u0)
    u1[indx_u].+=0.0001

    dummy_parameters1 = (
        d=d1,
        u=u1,
        g=g0
    )
    dummy_parameter_points = [parameters0,dummy_parameters1]

    P,var_lens,parameters = TR.forward_system(tg,TR.DeformMode();)
    Psys, ide_pindx = TR.find_diff_system(P,parameters,dummy_parameter_points)

    # precompile
    TR.forward_once(Psys,
                        var_lens,
                        [[startsols.q̌; startsols.s; startsols.λ]],
                        reduce(vcat,(d=d0[indx_d], u=u0, )),
                        reduce(vcat,(d=d1[indx_d], u=u1, ))
        )

    set_theme!(theme_try;fontsize=6 |> pt2px)
    fig = Figure(resolution=(1920,1080))
    ax = Axis(fig[1,1])
    ax.aspect = DataAspect()

    bars,cables = plotstructure!(ax,tg)

    sg_d = SliderGrid(fig[1,2],
        [
            (
                label = latexstring("d_{$i}"),
                range = LinRange(d0[i]-0.01,d0[i]+0.01,11),
                startvalue = d0[i]
            )
            for i in indx_d
        ]...
    )

    sg_a = SliderGrid(fig[1,3],
        [
            (
                label = latexstring("\\Delta u_{$i}"),
                range = LinRange(-0.01,0.01,11),
                startvalue = 0.0
            )
            for i in collect(1:2:2nsegs)
        ]...
    )
    colsize!(fig.layout,1,Relative(0.6))

    sg_values = vcat(
        [s.value for s in sg_d.sliders],
        [s.value for s in sg_a.sliders]
    )
    sg_lens = length.([sg_d.sliders,sg_a.sliders])
    onany(sg_values...) do vals_tuple...
        vals = [v for v in vals_tuple]
        dv, av = TR.split_by_lengths(vals,sg_lens)
        d1 .= d0
        d1[indx_d] = dv
        u1 .= u0
        # a1 = 0.01
        u1[collect(1:2:2nsegs)] .+= av
        u1[collect(2:2:2nsegs)] .-= av

        rst = TR.forward_once(Psys,
                            var_lens,
                            [[startsols.q̌; startsols.s; startsols.λ]],
                            reduce(vcat,(d=d0[indx_d], u=u0, )),
                            reduce(vcat,(d=d1[indx_d], u=u1, ))
            )
        rst_rc = TR.recover(rst,tg)
        update_scene!(tg,bars,cables,rst_rc.q)
    end
    xlims!(ax,-0.1,1.2)
    ylims!(ax,-1.1,0.2)
    fig
end
get_interactive(man_inv)
