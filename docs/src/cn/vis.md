# 可视化

```@meta
CurrentModule = Rible
ShareDefaultModule = true
```

```@setup
using Rible
import Rible as RB
using StaticArrays
using LinearAlgebra
using Makie
import CairoMakie
CairoMakie.activate!()
Makie.inline!(true)
using Rotations
using TypeSortedCollections
using FileIO
using MeshIO
import GeometryBasics as GB
```

Rible.jl 基于 Makie.jl 绘图生态提供了一套可视化工具链，读取由求解器生成的轨迹数据，支持二维与三维渲染，并提供开箱即用的交互式回放功能。

> **核心架构**：可视化模块作为独立的后处理环节，依赖于 `Robot` 对象的静态结构以及记录其中的数据轨迹。通过底层的绘图配方 (Recipe)，模块能自动识别并渲染各类物理组件。

## 核心配方 Recipe

可视化层建立在一个名为 `Vis` 的 Makie 配方 **`Makie.@recipe Vis (object,)`** 上。用户侧看到的 `vis` 与 `vis!`，本质上都是把对象送入这套 recipe 系统，再由其中保存的对象`object`类型决定后续具体走哪条渲染路径。

结构层的入口方法为
```julia
Makie.plot!(::Vis{Tuple{<:Structure}})
```

Structure 对应的方法并不会直接生成图元。它会先把输入包装成 `VizStructure`，再把这个内部对象交回 `vis!`。真正执行系统遍历的是 `VizStructure` 对应的方法：它会依次遍历其中的器件与物体，并把每个对象重新分派到更具体的绘图方法。

刚体层的可视化方法为
```julia
Makie.plot!(::Vis{Tuple{<:AbstractRigidBody}})
```

当分派到刚体方法后，它会先计算质心状态，再把质心和每个 `LocusState` 分别重新交给 `vis!` 处理。如果启用了网格渲染，它还会调用 `build_mesh`，并把结果交给 `Makie.poly!`。因此，刚体方法本身更像一个协调层：它负责把一个力学对象拆成多个更细的可视化目标。


器件层的可视化方法为
```julia
Makie.plot!(::Vis{Tuple{<:Apparatus}})
```

Rible 中对 `Apparatus` 的默认方法实际上不会画出内容，它只是一个占位入口。真正的器件渲染通常由扩展包补充。例如 `RibleTensegrity` 为 `Apparatus{<:CableJoint}` 和 `Apparatus{<:ClusterJoint}` 增加了专门的 `Makie.plot!` 方法，用线段和标签来表达缆索结构。这种设计让顶层 API 保持稳定，而具体力学系统仍然可以按包扩展自己的可视化分支。

位点层的可视化负责点、标签和可选的坐标轴箭头
```julia
Makie.plot!(::Vis{Tuple{<:LocusState}})
```

到了最底层，recipe 不再围绕系统组件继续分派，而是直接调用标准 Makie 图元接口绘图命令，如`scatter!`、`text!`、`arrows3d!`、`linesegments!`、`poly!`等等。

---

## 其他核心组件

除了 `recipe.jl` 提供的一整套对象渲染法则外，可视化模块还包含以下核心组件以提升其实用性：

- **交互式轨迹渲染 `plot_traj.jl`**：提供了顶层接口 `plot_traj!`。它能够直接读取 `Simulator` 的执行结果，并生成一个带有时间滑块的交互式图形界面，方便您观察系统随时间的动态演变。
- **数据分析与对比 `plot.jl`**：封装了轨迹抽取与误差分析函数，例如抽取速度 `get_velocity!` 或计算对比误差 `get_err_avg`。这些工具适合在算法性能测试或基准对比时快速绘制收敛曲线。
- **三维网格处理 `mesh.jl`**：借助 GeometryBasics.jl 生成并更新三维刚体的表面网格映射。

---

## 典型可视化流程

下面先构建一个旋转陀螺并执行短时积分。这样页面内后续所有可视化示例都可以复用同一个 `top` 与 `sim_result`。

```@example
include(joinpath(pathof(Rible), "../../examples/robots/spinningtop.jl"))

origin_position = [0.0, 0.0, 0.3]
R = RotX(0.0)
origin_velocity = [0.5, 0.0, 0.0]
Ω = [0.0, 0.0, 100.0]

top = make_top(origin_position, R, origin_velocity, Ω, :NCF; μ = 0.5, e = 0.5, loadmesh = true)
top.structure.connectivity.num_of_bodies
```


### 1. 局部结构静态渲染

如果在建模阶段需要检查初始位姿或拼装拓扑是否正确，可以直接调用通用的 `vis!` 指令。借助强大的多重分派，它可以直接接收最顶层的 `Structure`，然后逐层向下渲染其全貌：

```@example
fig_static = Figure()
ax_static = Axis3(fig_static[1, 1], aspect = :data)
Rible.vis!(ax_static, top.structure; show_loci = false, show_labels = false)
fig_static
```

### 2. 交互式轨迹回放

```@example
planes = StaticContactSurfaces([
    HalfSpace([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
])

contact_model = RestitutionFrictionCombined(
    NewtonRestitution(),
    CoulombFriction(),
)

prob = DynamicsProblem(top; env=planes, contact_model)
solver = DynamicsSolver(
    Zhong06(),
    InnerLayerContactSolver(InteriorPointMethod()),
)
sim_result = solve!(prob, solver; tspan = (0.0, 0.2), dt = 5e-4, ftol = 1e-14, maxiters = 50, exception = false)
```

完成动力学积分后，最直观的检查方式是调用 `plot_traj!`。它会自动提取状态序列，并构建一个可交互的时间流逝界面：


有了轨迹之后，`plot_traj!` 也可以同时展示多个时间截面。对静态文档页面而言，小型网格布局往往比滑块更适合直接比较不同时间步的姿态。

```@example
key_steps = round.(Int, range(1, length(top.traj.t), length = 4))

fig_traj = plot_traj!(
    top;
    do_slide = false,
    show_info = false,
    show_background = false,
    show_loci = false,
    show_labels = false,
    show_ground = false,
    gridsize = (2, 2),
    at_steps = key_steps,
)
```

---

## 性能与定制

Rible 通过更新 Makie 的 `Compute Pipeline` 来保持回放和多视图渲染的流畅性，而不是在每个时间步重复创建图形对象。视觉样式方面，它也尽量保持与标准 Makie 工作流一致，因此可以直接配合自定义 `Theme`、坐标轴设置和图布局使用。
