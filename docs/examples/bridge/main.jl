# # 三维纯直杆可展开张拉整体桥动力学仿真
# 本代码文件位置
print(@__FILE__)
# 加载所需程序包
using LinearAlgebra
using SparseArrays, StaticArrays, OffsetArrays, CircularArrays
using RecursiveArrayTools, StructArrays
using TypeSortedCollections
using Rotations, CoordinateTransformations
using Interpolations
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using LaTeXStrings, Unitful, Printf
using GeometryBasics, Meshing
using EponymTuples
using Match
using Revise
import Meshes
import TensegrityRobots as TR

const CA = CircularArray

# 加载所需源代码
includet("bridge.jl")
includet("../vis.jl")
includet("../analysis.jl")
includet("../dyn.jl")

## deployable bridge
n = 5; m = 2

# 初始构型
tgbridge1 = nrailbridge(n,m;θ=0.0,c=0.0,h=1.0)
# 初始构型可视化
plot_traj!(tgbridge1;AxisType=LScene,showground=false)

ℓ = TR.get_cables_len(tgbridge1.tg)
Y = TR.build_Y(tgbridge1)
ax = float.(Y)\ℓ
B,F̃ = TR.build_inverse_statics_for_actuation(tgbridge1,tgbridge1)
ap,an = TR.get_solution_set(B,F̃)
a1 = ap

# 目标构型
tgbridge0 = nrailbridge(n,m;r=0.7714704084729095,c=25.0,h=0.45,right=false)
# 目标构型可视化
plot_traj!(tgbridge0;AxisType=LScene,showground=false)

ℓ = TR.get_cables_len(tgbridge0.tg)
Y = TR.build_Y(tgbridge0)
ax = float.(Y)\ℓ
B,F̃ = TR.build_inverse_statics_for_actuation(tgbridge0,tgbridge0)
ap,an = TR.get_solution_set(B,F̃)
a0 = ap

# 计算驱动量
μ0,μ1 = Y*a0, Y*a1
TR.actuate!(tgbridge0,a0)
TR.update!(tgbridge0.tg)
TR.get_cables_tension(tgbridge0.tg)
_ = TR.check_static_equilibrium_output_multipliers(tgbridge0.tg)

# 预设驱动的桥
l = 0
actor = make_pres_actor(repeat(μ0,l+1),repeat(μ1,l+1),1.0,4.0)
twobridges = nrailbridge(n,m;r=0.7714704084729095,c=25.0,h=0.45,right=false)

dpbridge = TR.TensegrityRobot(twobridges.tg,(actuators=[actor],))
TR.actuate!(dpbridge,[0.0])
TR.update!(dpbridge.tg)
_ = TR.check_static_equilibrium_output_multipliers(dpbridge.tg)

# 动力学仿真
dt = 1e-3

TR.solve!(TR.SimProblem(dpbridge,
        (x)->dynfuncs(x;actuate=true,gravity=false)
        ),
        TR.Zhong06(),
        dt=dt,
        tspan=(0.0,5.0),
        ftol=1e-14
    )

# 仿真结果可视化
plot_traj!(dpbridge;AxisType=LScene,showground=false,actuate=true)
