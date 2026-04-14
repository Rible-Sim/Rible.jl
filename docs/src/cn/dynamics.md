# 机器人、动力学与仿真

```@meta
CurrentModule = Rible
ShareDefaultModule = true
```

```@setup
using Rible
import Rible as RB
using StaticArrays
using LinearAlgebra
using Rotations
using TypeSortedCollections
using Makie
import CairoMakie
import GeometryBasics as GB
import RungeKutta

CairoMakie.activate!()
Makie.inline!(true)
const Point3 = GB.Point3
```

## 机器人

机器人 `Robot` 结构体是物理系统的数字孪生。它负责存储机器人的所有静态属性与动态历史，是整个仿真框架的数据锚点。

它的核心构成包括：
- **`structure` (物理结构)**：维护多体系统的拓扑连接（详见 [结构](@ref "结构")），包含所有的刚柔体 (`Body`)、关节 (`Joint`) 与外力装置 (`Apparatus`)（详见 [体与器](@ref "体与器")）。
- **`hub` (控制枢纽)**：统一管理控制所需的执行器 (`Actuators`) 和量器 (`Gauges`)（详见 [控制枢纽](@ref "控制枢纽")）。
- **轨迹容器**：持久化存储仿真推演过程中的所有数据，包括 `traj`（状态轨迹）、`contacts_traj`（接触轨迹）和 `control_traj`（控制变量轨迹）。

---

## 动力学与仿真

Rible 提供了一个灵活的动力学仿真框架，支持多种积分算法、接触力学模型和控制策略。

> **核心架构**：在 Rible 中，**物理模型** (`DynamicsProblem`) 与 **计算策略** (`DynamicsSolver`) 分离。这意味着您可以对同一个物理系统无缝切换不同的求解算法，而无需修改任何模型定义。

## 核心组件

### 动力学问题 DynamicsProblem

`DynamicsProblem` 是仿真任务的数据容器。它是一个不可变结构体，用于绑定物理系统及其环境的所有静态属性，为求解器提供明确的计算边界。

它主要包含以下要素：
- **`Robot` (机器人)**：机械系统的拓扑结构、质量分布和运行时状态（详见 [机器人](@ref "机器人")）。
- **`Objective` (优化目标)**：用于伴随分析代价函数（详见 [误差、代价与灵敏度分析](@ref "误差、代价与灵敏度分析")）。
- **`Environment` (环境)**：提供几何空间与连续的被动外力。
  - **力场 (`Field`)**：提供全局外力，如重力场 (`Gravity`)。
  - **几何体 (`Geometry`)**：提供用于碰撞检测的表面或网格地形。
- **`Policy` (控制策略)**：将当前系统状态 ``(t, q, v)`` 映射为主动控制力 `u`。
- **`Models` (交互模型)**：约束物理交互的数学形式，如接触模型（`AbstractContactModel`）。

### 动力学求解器 DynamicsSolver

如果说 `DynamicsProblem` 负责定义“是什么”，那么 `DynamicsSolver` 则负责“怎么算”。它本身不存储物理状态，而是将不同的数值积分器与子求解器组合起来，用于计算系统的演化。

- **积分器 (`Integrator`)**：指明积分格式，即计算离散时间步上的下一个系统状态 ``(q_{k+1}, v_{k+1})``。
- **子求解器**：包含物体求解器、装置求解器以及接触求解器。

---

## 积分器选择

积分器是 `DynamicsSolver` 的数值引擎。选择合适的积分器对仿真的稳定性与精度至关重要。

### Zhong 06
对于大多数带约束的多体系统，这是默认推荐的积分器。最早由钟万勰提出，后改进并用于接触动力学 [luoNonsmoothModifiedSymplectic2024](@cite)。它是一种二阶、保约束且类辛的格式，能够长时间保持能量守恒与约束流形。

```@example
# 使用默认设置构造
integrator = Zhong06()
solver = DynamicsSolver(integrator)
```

### 龙格-库塔法
通过 `RungeKutta.jl` 支持各种隐式 Runge-Kutta 格式，如 Gauss-Legendre。隐式 RK 方法非常适合对平滑系统有高精度要求的场景。

```@example
using RungeKutta
# 构造一个二阶 Gauss-Legendre (隐式) 积分器
integrator = RKIntegrator(RungeKutta.TableauGauss(2))
solver = DynamicsSolver(integrator)
```

### 额外积分器 (Generalized-α, Moreau-Jean)
为满足对比验证的需要，Rible 在 [RibleExtraIntegrators](@ref "额外积分器 RibleExtraIntegrators") 扩展包中提供了更多积分器。

---

## 典型仿真流程

### 1. 便捷 API：`solve!`

对于大多数标准场景，直接使用 `solve!` 方法即可。它会在底层自动完成内存分配并推演主循环。

下面的示例复用了仓库中的旋转陀螺模型。该陀螺具有较高的初始自旋速度，会与平面地面发生接触，并生成可供后续检查的轨迹。

```@example
include(joinpath(pathof(Rible), "../../examples/robots/spinningtop.jl"))

origin_position = [0.0, 0.0, 0.5]
R = RotX(0.0)
origin_velocity = [1.0, 0.0, 0.0]
Ω = [0.0, 0.0, 200.0]

top = make_top(origin_position, R, origin_velocity, Ω, :NCF; μ = 0.95, e = 0.5, loadmesh = false)
top.structure.connectivity.num_of_bodies
```

这个模型已经内置了 `Structure` 与 `ControlHub`，下一步只需定义接触表面与交互模型。

```@example
planes = StaticContactSurfaces([
    HalfSpace([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
])

contact_model = RestitutionFrictionCombined(
    NewtonRestitution(),
    CoulombFriction(),
)

prob = DynamicsProblem(top; env=planes, contact_model)
```

准备完成后，再定义带接触求解器的 `DynamicsSolver` 并执行一次短时积分。

```@example
solver = DynamicsSolver(
    Zhong06(),
    InnerLayerContactSolver(InteriorPointMethod()),
)

sim_result = solve!(
    prob,
    solver;
    tspan = (0.0, 0.2),
    dt = 5e-4,
    ftol = 1e-14,
    maxiters = 50,
    exception = false,
)

tip_position = get_trajectory!(top, 1, 1).u[end]
tip_velocity = get_velocity!(top, 1, 1).u[end]
tip_position, tip_velocity
```

### 2. 仿真容器 Simulator

`Simulator` 是负责维护运行状态的执行容器。`DynamicsProblem` 和 `DynamicsSolver` 本身是无状态的。

当你需要精细控制主循环（如结合外部强化学习框架、添加运行时控制器或实现异常热重启）时，应该显式创建 `Simulator` 并逐步执行：

下面先创建 `Simulator`，再调用针对容器的 `solve!`。

```@example
top_restart = make_top(origin_position, R, origin_velocity, Ω, :NCF; μ = 0.95, e = 0.5, loadmesh = false)
prob_restart = DynamicsProblem(top_restart; env=planes, contact_model)

sim = Simulator(prob_restart, solver; tspan = (0.0, 0.02), dt = 5e-4, restart = true)
solve!(sim, solver; dt = 5e-4)
traj_len_sim = length(top_restart.traj)
t_end_sim = top_restart.traj.t[end]
traj_len_sim, t_end_sim
```

完成积分后，可以直接从陀螺模型的轨迹容器提取后处理指标。

`mechanical_energy!` 用于检查积分过程中的能量行为，常用于快速诊断动力学是否符合预期。

```@example
E = mechanical_energy!(top).E
E[begin], E[end]
```

如果希望直接查看路径，可以先用 `get_trajectory!` 提取一个点的空间轨迹，再用 `Makie.lines` 绘制三维曲线。

```@example
tip_traj = get_trajectory!(top, 1, 1).u
tip_points = Point3.(tip_traj)

fig_traj = Figure()
ax_traj = Axis3(fig_traj[1, 1], aspect = :data)
Makie.lines!(ax_traj, tip_points, linewidth = 2)
fig_traj
```

此外，轨迹数据也可以继续对接到[可视化](vis.md)页面中的工作流，生成更丰富的动画与图形。

---

## 核心机制与性能优化

Rible 利用多重分派 (Multiple Dispatch) 与特性系统 (Traits) 来优化计算性能。

1. **轨迹内存预分配 (`prepare_traj!`)**
   初始化 `Simulator` 时自动调用。它会为整个时间段内的状态序列预先分配连续的内存块。

2. **工作空间自动生成 (`generate_cache`)**
   仿真开始前，根据问题与求解器的具体类型，在编译期生成专属的工作空间 (`Workspace`)，用于缓存矩阵（如雅可比矩阵）。

3. **基于特性 (Traits) 的分派优化**
    通过 `has_constant_mass_matrix(bot)` 判断质量矩阵是否为常数。若是，将自动跳过所有时间步的质量矩阵重复计算。

## 参考文献

```@bibliography
Pages = ["dynamics.md"]
```

