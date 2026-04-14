# # Rible.jl 入门

# Rible.jl 是面向刚柔机器人的多体动力学仿真平台。
# 本页通过一个滚动 superball 张拉整体机器人示例，演示从建模到求解的完整流程。

# ## 安装

# 在 Julia REPL 中执行以下命令安装 Rible.jl：

# ```julia
# import Pkg
# Pkg.activate(; temp=true)
# Pkg.add(url="https://github.com/Rible-Sim/Rible.jl")
# Pkg.add(url="https://github.com/Rible-Sim/Rible.jl", subdir="RibleTensegrity")
# ```

# > **注意**：Rible.jl 需要 Julia ≥ 1.11。

# ## 导入

using Rible
import Rible as RB
using RibleTensegrity
import RibleTensegrity as RT
using CircularArrays
using StaticArrays
using Rotations
using LinearAlgebra
using TypeSortedCollections

# ## 构建模型

# Rible.jl 在 `examples/` 目录下提供了多个可复用示例定义。
# 以下加载 superball 张拉整体机器人及其刚杆辅助定义：

include(joinpath(pathof(Rible), "../../examples/bodies/rigidbar_nonsmooth_repro.jl"));
include(joinpath(pathof(Rible), "../../examples/robots/superball.jl"));

# 创建一个具有适中滚动初始条件的 superball，并关闭外部网格加载，使示例保持自包含：

l = 1.7 / 2
d = l / 2

bot = superball(0.0;
    origin_velocity = SVector(2.0, 1.0, 0.0),
    ω = SVector(0.0, 0.0, 1.0),
    μ = 0.05,
    e = 0.0,
    l,
    d,
    z0 = l^2 / (sqrt(5) * d) - 1e-3,
    loadmesh = false,
    visible = true,
    constrained = false,
)

# ## 定义仿真问题

# `DynamicsProblem` 将机器人模型、物理环境和接触模型封装为完整的仿真问题。

# 定义地面碰撞面——一个 z 轴正方向的半空间：

ground_plane = RB.StaticContactSurfaces([
    RB.HalfSpace([0.0, 0.0, 1.0], zeros(3))
]);

# 选择接触模型（牛顿恢复系数 + 库仑摩擦）：

contact_model = RB.RestitutionFrictionCombined(
    RB.NewtonRestitution(),
    RB.CoulombFriction(),
);

# 组装 `DynamicsProblem`：

prob = RB.DynamicsProblem(bot; env=ground_plane, contact_model=contact_model)

# ## 求解

# 使用 `Zhong06` 积分器搭配内层接触求解器，运行一段简短仿真：

solver = RB.DynamicsSolver(
    RB.Zhong06(),
    RB.InnerLayerContactSolver(RB.InteriorPointMethod()),
)

result = RB.solve!(prob, solver; tspan=(0.0, 1.0), dt=1e-2, ftol=1e-12, maxiters=200, exception=false);

# 仿真轨迹存储在 `bot` 的内部数据结构中，可通过 `result` 访问完整的仿真状态。

# ## 可视化

# 使用轨迹绘图接口查看滚动过程：

using CairoMakie
plot_traj!(bot; do_slide = false, show_loci = false, show_mesh = true, show_ground = false)

# ## 下一步

# 至此，你已完成第一个 Rible.jl 动力学仿真。继续浏览文档以了解更多：

# - **建模基础**：[体与器](body_and_apparatus.md)、[结构](structure.md) 和 [控制枢纽](hub.md)
# - **动力学计算**：[机器人动力学计算](dynamics.md) 和 [误差、代价与灵敏度](adjoint/index.md)
# - **可视化**：[可视化](vis.md)
