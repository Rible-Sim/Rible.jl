using LinearAlgebra
using Statistics
using SparseArrays
using StaticArrays
using TypeSortedCollections
using Rotations
# import GeometricalPredicates as GP
using CoordinateTransformations
using OffsetArrays
import DifferentialEquations as DE
using BenchmarkTools
using RecursiveArrayTools
using Interpolations
# using CubicSplines
# import FLOWMath
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
import Meshes
using LaTeXStrings
using TexTables
using Unitful
using EzXML
using CSV, Tables
using Printf
using StructArrays
using EponymTuples
using GeometryBasics, Meshing
# using FloatingTableView
using Match
using Cthulhu
using Revise
import Rible as RB
cd("examples/test_prescribed")

include("../vis.jl")
include("../analysis.jl")
include("../dyn.jl")
include("dyn.jl")
include("def.jl")
include("def3d.jl")

includet("../vis.jl")
includet("../analysis.jl")
includet("../dyn.jl")
includet("dyn.jl")
includet("def.jl")
includet("def3d.jl")

# figdir::String = raw"C:\Users\luo22\OneDrive\Papers\DynamicTensegrity\CS"
figdir::String = raw"C:\Users\luo22\OneDrive\Papers\Ph.D.Thesis\ff"


function inv!(bots;)
    for (i, bot) in enumerate(bots)
        μ = RB.inverse_for_restlength(bot, bot)
        RB.set_restlen!(bot.st,μ)
        RB.update!(bot.st)
        # κ = RB.get_cables_stiffness(bot.st)
        # ℓ = RB.get_cables_len(bot.st)
        # push!(ret,μ)
        if i > 1
            push!(bots[1].traj,bot.st.state.system)
        end
    end
    bots
end

function plotsave_inv_restlen(bots,xs,figname=nothing;
        xlabel = L"x",
        groups = [
            i for i = 1:8
        ]
    )
    friction_coefficients = RB.get_cables_restlen.(bots) |> VectorOfArray
    fs = RB.get_cables_tension.(bots) |> VectorOfArray

    @show fs[round(Int,length(fs)/2)]
    with_theme(theme_pub;
        figure_padding = (0,fontsize,0,0),
        resolution = (0.9tw,0.23tw),
        markersize = 0.8fontsize,
        palette = (
            marker = [
                :xcross,:cross,
                :utriangle,:dtriangle,
                :ltriangle,:rtriangle,
                :diamond,:hexagon,
                :star8,:star5
            ],
            color = cgrad(:seaborn_bright, 10, categorical=true).colors.colors,
        ),
        Scatter = (
            cycle = Cycle(
                [:color, :marker], 
                covary = true
            ),
        ),
        Lines = (cycle = [:color],),
    ) do 
        fig = Figure()
        ax1 = Axis(fig[1,1];xlabel, ylabel = L"\mu~(\mathrm{m})")
        Label(fig[1,1,TopLeft()],"($(alphabet[1]))")
        for j in eachindex(groups)
            lines!(ax1,xs,friction_coefficients[groups[j][1],:];)
            scatter!(ax1,xs,friction_coefficients[groups[j][1],:];label = "No. $(string(groups[j]))")
        end
        ax2 = Axis(fig[1,2];xlabel, ylabel = L"f~(\mathrm{N})")
        Label(fig[1,2,TopLeft()],"($(alphabet[2]))")
        for j in eachindex(groups)
            lines!(ax2,xs,fs[groups[j][1],:];)
            scatter!(ax2,xs,fs[groups[j][1],:];label = "No. $(string(groups[j]))")
        end
        ylims!(ax2,0,(maximum.(fs) |> maximum)+1)
        Legend(fig[1,3],ax2,)
        savefig(fig,figname)
        fig
    end
end

c1bot = class1(;)

GM.activate!(); with_theme(theme_pub;
        fontsize = 8 |> pt2px,
        figure_padding = (fontsize,0,fontsize,0),
    ) do
    plot_traj!(c1bot;
        figsize = (0.45tw,0.32tw),
        AxisType=Axis3,
        doslide=false,
        showinfo=false,
        showground=false,
        showlabels=false,
        showpoints=false,
        xlims=(-0.1, 0.15),
        ylims=(-0.1, 0.15),
        zlims=(-0.01, 0.18),
        (sup!)=(ax, _) -> begin
            ax.azimuth = 6.655530633326979
            ax.elevation = 0.2526990816987239
            ax.title = ""
            ax.xlabeloffset = 2fontsize
            ax.ylabeloffset = 1fontsize
            ax.zlabeloffset = 2.5fontsize
        end,
        figname = "class1"
    )
end

# pitch
xs = 0:1.5:30 |> collect
xs = 0:15:30 |> collect # vis
ψθϕs = [
  [0,x,0].|> deg2rad
  for x in xs
]
xlabel = L"\theta~(\degree)"
figsuffix = "theta"

# roll
xs = 0:1.5:30 |> collect
xs = 0:15:30 |> collect # vis
ψθϕs = [
  [0,0,x].|> deg2rad
  for x in xs
]
xlabel = L"\phi~(\degree)"
figsuffix = "phi"

class1s = inv!([
    class1(; R1=RotZYX(ψθϕ...))
    for ψθϕ in ψθϕs
];)

GM.activate!(); with_theme(theme_pub;
        figure_padding = (fontsize,0,fontsize,0),
    ) do
    plot_traj!(class1s[1];
        figsize = (0.50tw,0.3tw),
        AxisType=Axis3,
        atsteps=collect(1:length(class1s)) |> reverse,
        # atsteps = [1,4],
        doslide=false,
        showinfo=false,
        showground=false,
        showpoints=false,
        showlabels=false,
        xlims=(-0.1, 0.15),
        ylims=(-0.1, 0.15),
        zlims=(-0.01, 0.18),
        colorbar=(
            colormap=cgrad(:RdYlGn, length(class1s), categorical=true),
            label="Deg.",
            limits=(xs[begin]-xs[begin+1]/2, xs[end]+xs[begin+1]/2),
            ticks=(xs,string.(xs)|> reverse),
            # flipaxis = true,
        ),
        (sup!)=(ax, _) -> begin
            ax.azimuth = 6.655530633326979
            ax.elevation = 0.2526990816987239
            ax.title = ""
            ax.xlabeloffset = 2fontsize
            ax.ylabeloffset = 1fontsize
            ax.zlabeloffset = 3fontsize
        end,
        figname = "class1_$(figsuffix)_vis"
    )
end

GM.activate!(); plotsave_inv_restlen(class1s,xs; 
    xlabel,
    groups = [
        i for i = 1:6
    ]
)
CM.activate!(); plotsave_inv_restlen(class1s,xs,"class1_$figsuffix"; 
    xlabel,
    groups = [
        i for i = 1:6
    ],
)

spinebot = spine3d(2; c=0.0)

GM.activate!(); with_theme(theme_pub;
        fontsize = 8 |> pt2px,
        figure_padding = (fontsize,0,fontsize,0),
    ) do
    plot_traj!(spinebot;
        figsize = (0.4tw,0.36tw),
        AxisType=Axis3,
        doslide=false,
        showinfo=false,
        showground=false,
        showlabels=false,
        xlims=(-0.05, 0.05),
        ylims=(-0.05, 0.05),
        zlims=(-0.03, 0.07),
        (sup!)=(ax,_,_) -> begin
            ax.azimuth = 3.4555306333269855
            ax.elevation = 0.4426990816987242
            ax.title = ""
            ax.xlabeloffset = 2fontsize
            ax.ylabeloffset = 1fontsize
            ax.zlabeloffset = 2.5fontsize
        end,
        figname = "ultra"
        # figname = "ultra_rot_vis"
    )
end

# yaw
xs = 0:3:60 |> collect
xs = 0:30:60 |> collect # vis
ψθϕs = [
  [x,0,0] .|> deg2rad
  for x in xs
]
xlabel = L"\psi~(\degree)"
figsuffix = "psi"

# roll
xs = 0:3:60 |> collect
xs = 0:30:60 |> collect # vis
ψθϕs = [
  [0,0,x].|> deg2rad
  for x in xs
]
xlabel = L"\phi~(\degree)"
figsuffix = "phi"

# yaw, roll
xs = 0:2.5:40 |> collect
xs = 0:20:40 |> collect # vis
ψθϕs = [
  [x,0,-x].|> deg2rad
  for x in xs
]
xlabel = L"\psi,~-\phi~(\degree)"
figsuffix = "psi,-phi"

spine3ds = inv!([
    spine3d(2; RR=RotZYX(ψθϕ...))
    for ψθϕ in ψθϕs
];)


GM.activate!(); with_theme(theme_pub;
        figure_padding = (fontsize,0,fontsize,0),
    ) do
    plot_traj!(spine3ds[1];
        figsize = (0.42tw,0.3tw),
        AxisType=Axis3,
        atsteps=collect(1:length(spine3ds)) |> reverse,
        # atsteps = [1,4],
        doslide=false,
        showinfo=false,
        showground=false,
        showlabels=false,
        xlims=(-0.05, 0.05),
        ylims=(-0.05, 0.05),
        zlims=(-0.03, 0.07),
        colorbar=(
            colormap=cgrad(:RdYlGn, length(spine3ds), categorical=true),
            label="Deg.",
            limits=(xs[begin]-xs[begin+1]/2, xs[end]+xs[begin+1]/2),
            ticks=(xs,string.(xs)|> reverse),
            # flipaxis = true,
        ),
        (sup!)=(ax,_, _) -> begin
            ax.azimuth = 3.4555306333269855
            ax.elevation = 0.4426990816987242
            ax.title = ""
            ax.xlabeloffset = 2fontsize
            ax.ylabeloffset = 1fontsize
            ax.zlabeloffset = 3fontsize
        end,
        figname = "ultra_$(figsuffix)_vis"
    )
end

groups = 1:8

GM.activate!(); plotsave_inv_restlen(spine3ds,xs; xlabel,groups)
CM.activate!(); plotsave_inv_restlen(spine3ds,xs,"ultra_$figsuffix";xlabel,groups)


