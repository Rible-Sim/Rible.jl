# 误差、代价与目标

```@meta
CurrentModule = Rible
ShareDefaultModule = true
```

```@setup
using Rible
import Rible as RB
using StaticArrays
using LinearAlgebra
import TypeSortedCollections as TSC
```

在机器人的控制设计、轨迹优化（如 iLQR）和强化学习训练中，我们需要一个定量的标量指标来评价系统的表现。Rible.jl 通过 `AbstractObjective` 体系将散落在各处的物理观测（Gauges）与执行消耗（Actions）聚合为统一的数学目标。

## 1. 从误差到代价的链条

在 Rible 中，**误差 (Error)** 是在量器层面定义的局部标量，而 **代价 (Cost)** 是在系统层级定义的加权全局指标。

### 1.1 ErrorGauge：评价的原子单元

正如在 [控制枢纽](../hub.md) 中所提到的，`ErrorGauge` 负责计算单一维度的偏差。例如：

- “末端执行器距离目标的欧式距离”
- “极杆偏离垂直线的角度”

每一个 `ErrorGauge` 产生的都是一个可微分的**标量误差**。

### 1.2 Objective：权重的指挥棒

`Objective` 的核心作用是“加权求和”。它定义了哪些误差更重要，以及需要付出多少控制代价。

总代价 ``J`` 通常由两部分组成：

```math
J = \underbrace{\sum_{k=0}^{N-1} L(x_k, u_k)}_{\text{运行代价 (Trajectory Cost)}} + \underbrace{\Phi(x_N)}_{\text{终端代价 (Terminal Cost)}} 
```

其中瞬时代价 ``L`` 是所有激活的 `ErrorGauge` 和 `Actuator` 的加权和。

## 2. 目标函数 (Objective) 的结构

一个典型的 `Objective` 对象包含了以下关键信息：

- **运行代价权重 (Trajectory Weights)**：每个时间步都会计算的权重。
- **终端代价权重 (Terminal Weights)**：仅在轨迹最后一个时间步计算。通常用于强化对“最终目标”的达成（如“必须停在特定位置”）。
- **致动器权重 (Actuator Weights)**：对控制输入 ``u`` 的大小进行惩罚（通常对应能量消耗）。
- **时间缩放因子 (``\eta``)**：用于处理变步长或特定积分方案下的权重调整。

### 核心接口：代价计算

基于 `Objective` 的评价通过以下函数实现：

- **`cost!(bot, objt, inst_state, u; mode=:trajectory)`**: 计算单步瞬时代价。
- **`cost!(bot, objt, solver, dt)`**: 计算整条轨迹的总代价。

## 3. 灵敏度分析：微分接口

对于轨迹优化和基于梯度的算法（如伴随灵敏度分析），仅有代价数值是不够的。Rible 提供了高效的缓存无关（Buffer-free）微分接口：

### 梯度 (Gradient)

**`cost_gradient!`**：计算代价关于以下变量的梯度：

- 坐标向量 ``q`` 和动量向量 ``p``。
- 控制输入向量 ``u``。
- 系统参数 ``c`` 或策略参数 ``\theta``。

### 黑塞矩阵 (Hessian)

**`cost_hessian!`**：计算代价关于状态和控制的二阶偏导数。这在 iLQR 等算法中用于构建局部二次近似。

> [!TIP]
> 这些微分接口底层自动利用了 `measure_jacobian!`。由于量器层级已经实现了高效的计算逻辑，`Objective` 层级的微分几乎不产生额外的计算开销。

## 4. 示例：Documenter 可运行的 Objective 与 cost! 流程

下面的示例全部使用 `@example`，可由 Documenter 直接执行并内联输出。

### 步骤 1：先构建一个最小可运行的 `Structure`

这一步的目标是准备一个可被 `cost!` 评估的最小物理系统。我们会依次创建 `Locus`、`RigidBodyProperty`、`RigidBodyState` 与 `NCF.NC1P3V`，然后组装成 `Structure` 并执行一次 `update!`。

这样做的目的是让后续量器与目标函数示例只聚焦在“代价计算逻辑”，而不被复杂模型细节干扰。

```@example
mass_locus = RB.Locus(SVector(0.0, 0.0, 0.0))
loci = [RB.Locus(SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0))]
inertia = SMatrix{3,3,Float64}(I)

prop = RB.RigidBodyProperty(1, true, 1.0, inertia, mass_locus, loci; visible = false)
r0 = SVector(0.0, 0.0, 0.0)
R0 = SMatrix{3,3,Float64}(I)
state = RB.RigidBodyState(prop, r0, R0, zero(r0), zero(r0))
coords = RB.NCF.NC1P3V(r0, r0, R0)
body = RB.RigidBody(prop, state, coords, nothing)
body.prop.id, length(body.prop.loci)
```

上一个代码块完成了单体建模；下面把它放进容器并初始化结构状态。

```@example
bodies = TSC.TypeSortedCollection([body])
apparatuses = TSC.TypeSortedCollection(Int[])
st = RB.Structure(bodies, apparatuses, RB.Connectivity(bodies, apparatuses))
RB.update!(st)
n_q = RB.get_num_of_coords(st)
n_c = RB.get_num_of_cstr(st)
n_q, n_c
```

### 步骤 2：定义一个误差规并组装 `ControlHub`

`ErrorGauge` 负责把“当前状态与参考值的偏差”映射成一个可微分标量；`Signifier` 指定测量作用在哪个 body 的哪个 locus。随后通过 `Coalition` 与 `ControlHub` 建立误差/动作在系统向量中的索引关系，再构造 `Robot` 作为统一入口。

```@example
error_gauge = RB.ErrorGauge(
    1,
    RB.Signifier(body, 1),
    RB.PositionCaptum(),
    [0.2, 0.0, 0.0],
)

capta_gauges = TSC.TypeSortedCollection(Int[])
error_gauges = TSC.TypeSortedCollection([error_gauge])
actuators = TSC.TypeSortedCollection(Int[])
num_error_gauges = length(error_gauges)
error_gauge.id, num_error_gauges
```

接着把量器、致动器和结构绑定为控制枢纽，并生成可直接用于 `cost!` 的 `Robot`。

```@example
coalition = RB.Coalition(st, capta_gauges, error_gauges, actuators)
hub = RB.ControlHub(st, capta_gauges, error_gauges, actuators, coalition)
bot = RB.Robot(st, hub)
num_errors = length(error_gauges)
num_actions = RB.get_num_of_actions(bot)
num_errors, num_actions
```

### 步骤 3：构造 `Objective` 并检查时间权重

这里的四个权重向量分别对应：运行误差项、运行控制项、终端误差项、终端控制项。由于这四个参数在构造器中的位置语义并不自解释，示例里保留了最小必要注释来标明每个向量的角色。

```@example
objt = RB.Objective(
    [10.0],      # trajectory_error_gauges_weights
    Float64[],   # trajectory_actuators_weights
    [100.0],     # terminal_error_gauges_weights
    Float64[],   # terminal_actuators_weights
    t -> 1 / (1 + t),
)

ηs = RB.get_trajectory_cost_weights(objt, [0.0, 0.1, 0.2], 0.1)
ηs
```

### 步骤 4：计算单步 trajectory / terminal cost

下面直接调用单步接口 `cost!(bot, objt, inst_state, u; mode=...)`。`mode=:trajectory` 用运行权重，`mode=:terminal` 用终端权重；并排展示两者可以直观看到同一状态在两套权重下的代价值差异。

```@example
inst = bot.structure.state.system
u = Float64[]

ϕ_traj = RB.cost!(bot, objt, inst, u; mode = :trajectory)
ϕ_term = RB.cost!(bot, objt, inst, u; mode = :terminal)
ϕ_traj, ϕ_term
```

## 5. 总结

通过将“测量什么”（Gauge）与“占比多少”（Objective）解耦，Rible 允许开发者：

1. **快速调参**：只需修改 `Objective` 里的权重数组，即可改变机器人的行为风格。
2. **多目标切换**：在同一个物理模型上轻松切换“速度优先”或“平稳优先”的任务目标。
3. **支持高级算法**：通过内置的 `cost_gradient!` 和 `cost_hessian!`，无缝对接高性能优化器。
