# 误差、代价与灵敏度分析

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

> [!WARNING]
> 灵敏度分析相关接口仍在开发中。本页描述的相关API仍可能于后续或针对不同求解器而调整。

在轨迹优化（如 iLQR）和强化学习中，我们需要一个标量指标来评价系统表现。Rible 通过 `AbstractObjective` 体系将物理观测（Gauges）与执行消耗（Actions）聚合为统一的数学目标。灵敏度分析则计算该目标对参数的梯度，支撑梯度类优化算法。

## 从误差到代价

在 Rible 中，**误差**是在量器层面定义的局部标量，而**代价**是在系统层级定义的加权全局指标。

### 误差量器

如 [控制枢纽](../hub.md) 所述，`ErrorGauge` 计算单一维度的偏差——例如末端执行器到目标的欧氏距离，或极杆偏离垂直线的角度。每个 `ErrorGauge` 产生一个可微分的标量误差。

### 目标函数

`Objective` 对所有激活的量器和致动器执行加权求和。总代价 ``J`` 由两部分组成：

```math
J = \underbrace{\sum_{k=0}^{N-1} L(x_k, u_k)}_{\text{轨迹代价}} + \underbrace{\Phi(x_N)}_{\text{终端代价}}
```

一个典型的 `Objective` 对象包含：

- **轨迹误差权重** — 每个时间步施加。
- **终端误差权重** — 仅在最终时间步施加，强调目标达成。
- **致动器权重** — 对控制输入 ``u`` 大小的惩罚（对应能量消耗）。
- **时间缩放因子** (``\eta``) — 处理变步长或特定积分方案下的权重调整。

核心评价接口：

- **`cost!(bot, objt, inst_state, u; mode=:trajectory)`** — 计算单步瞬时代价。
- **`cost!(bot, objt, solver, dt)`** — 计算整条轨迹的总代价。

### 微分接口

为梯度类算法提供无缓存微分：

- **`cost_gradient!`** — 代价关于 ``q``、``p``、``u`` 和参数 ``\theta`` 的梯度。
- **`cost_hessian!`** — 二阶导数，用于局部二次近似（如 iLQR）。

> [!TIP]
> 这些接口底层自动调用 `measure_jacobian!`。

## 灵敏度分析概述

灵敏度分析度量系统输出对输入与参数变化的响应程度。Rible 提供两类方法，均兼容 Zhong06 离散积分框架：

| 特性 | 伴随法 | 直接法 |
| :--- | :--- | :--- |
| 计算方向 | 反向模式 | 前向模式 |
| 适用场景 | 参数多、输出少 | 参数少、状态维度高 |
| 计算复杂度 | ``O(N_{obj})``，与参数数目弱相关 | ``O(N_{param})``，与参数数目线性 |
| 内存消耗 | 较高（存储前向轨迹） | 较低（在线推进） |
| 典型应用 | 大规模优化、策略训练 | 实时控制、在线监控 |

- 参数维度高（如神经网络）且目标为标量时，优先用**伴随法**。
- 参数少且需要在线灵敏度轨迹时，优先用**直接法**。

## 伴随灵敏度分析

伴随法通过反向递推变量，一次反向扫描即可得到目标函数对大规模参数的梯度。

### 理论基础

考虑离散系统：

```math
x_{k+1} = f(x_k, u_k, \theta), \quad k = 0, \dots, N-1
```

目标函数：

```math
J = \sum_{k=0}^{N-1} L_k(x_k, u_k) + \Phi(x_N)
```

引入伴随变量 ``\lambda_k``，反向递推：

```math
\lambda_k = \frac{\partial L_k}{\partial x_k} + \left(\frac{\partial f}{\partial x_k}\right)^T \lambda_{k+1}
```

参数梯度为：

```math
\frac{dJ}{d\theta} = \sum_{k=0}^{N-1} \left(\frac{\partial L_k}{\partial \theta} + \lambda_{k+1}^T \frac{\partial f}{\partial \theta}\right)
```

## 直接灵敏度分析

直接法将状态方程与灵敏度方程联立，前向推进状态与敏感度。

### 理论基础

定义灵敏度矩阵 ``S = \partial x / \partial \theta``，灵敏度方程为：

```math
\dot{S} = \frac{\partial f}{\partial x} S + \frac{\partial f}{\partial \theta}
```

可与原方程同时前向积分，得到 ``x(t)`` 与 ``S(t)`` 的全时域轨迹。

## 灵敏度分析接口

本页优先介绍公共工作流：先计算目标函数，再构造灵敏度求解器，随后调用 `solve!`，最后检查 `sim.solver_cache`。具体缓存类型名属于求解器族内部实现，只有在扩展底层算法时才需要关心。

灵敏度求解器本质上是在已有动力学求解器和目标函数之外再包一层。

- **伴随法** — `adj_solver = DiscreteAdjointDynamicsSolver(dyn_solver, objt)`，随后 `adj_sen_solver = AdjointDynamicsSensitivitySolver(dyn_solver, adj_solver)`。
- **直接法** — `drc_sen_solver = DirectDynamicsSensitivitySolver(dyn_solver, objt)`。
- **执行入口** — `solve!(prob, sensitivity_solver; ...)` 会返回一个仿真结果，灵敏度输出存放在 `sim.solver_cache`。

这样做的好处是它延续了动力学模块的公共接口：`DynamicsProblem` 负责物理模型，`DynamicsSolver` 负责前向积分，而灵敏度分析只是额外叠加的一层求解策略。

求解完成后，优先从 `sim.solver_cache` 读取结果，而不是直接依赖某个求解器族的缓存类型名。

- **伴随梯度** — `∂J∂x₀ᵀ`、`∂J∂θᵀ`、`∂J∂cᵀ`
- **直接法雅可比** — `Jac_state`、`Jac_action`、`Jac_control_params`

伴随法字段描述标量目标函数对初始状态和参数组的变化率。直接法字段描述传播后状态对初始状态、驱动量和控制参数的变化率。


## 示例：代价与灵敏度工作流

### 步骤 1：构建最小系统

创建一个含误差量器和控制枢纽的单体结构。

```@example
using Rible
using StaticArrays
using LinearAlgebra
import TypeSortedCollections as TSC

mass_locus = Locus(SVector(0.0, 0.0, 0.0))
loci = [Locus(SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0))]
inertia = SMatrix{3,3,Float64}(I)

prop = RigidBodyProperty(1, true, 1.0, inertia, mass_locus, loci; visible = false)
r0 = SVector(0.0, 0.0, 0.0)
R0 = SMatrix{3,3,Float64}(I)
state = RigidBodyState(prop, r0, R0, zero(r0), zero(r0))
coords = NCF.NC1P3V(r0, r0, R0)
body = RigidBody(prop, state, coords, nothing)

bodies = TSC.TypeSortedCollection([body])
apparatuses = TSC.TypeSortedCollection(Int[])
st = Structure(bodies, apparatuses, Connectivity(bodies, apparatuses))
update!(st)
```

### 步骤 2：定义量器、枢纽和目标函数

```@example
error_gauge = ErrorGauge(
    1,
    Signifier(body, 1),
    PositionCaptum(),
    [0.1, 0.0, 0.0],
)

capta_gauges = TSC.TypeSortedCollection(Int[])
error_gauges = TSC.TypeSortedCollection([error_gauge])
force_actuator = ExternalForceActuator(
    1,
    Signifier(body, 1),
    NaiveOperator(2),
    [1.0 0.0; 0.0 1.0; 0.0 0.0],
    [0.0, 0.0],
)
actuators = TSC.TypeSortedCollection([force_actuator])

coalition = Coalition(st, capta_gauges, error_gauges, actuators)
hub = ControlHub(st, capta_gauges, error_gauges, actuators, coalition)
bot = Robot(st, hub)
prob = DynamicsProblem(bot; env=GravityEnv())

objt = Objective(
    [10.0],      # 轨迹误差权重
    [1.0],       # 轨迹控制权重
    [100.0],     # 终端误差权重
    [5.0],       # 终端控制权重
    t -> 1 / (1 + t),
)
```

### 步骤 3：计算代价

```@example
inst = bot.structure.state.system
u = [0.2, -0.1]

ϕ_traj = cost!(bot, objt, inst, u; mode = :trajectory)
ϕ_term = cost!(bot, objt, inst, u; mode = :terminal)
ϕ_traj, ϕ_term
```

### 步骤 4：伴随灵敏度求解

```@example
dyn_solver = DynamicsSolver(Zhong06())
adj_solver = DiscreteAdjointDynamicsSolver(dyn_solver, objt)
adj_sen_solver = AdjointDynamicsSensitivitySolver(dyn_solver, adj_solver)

adj_sim = solve!(prob, adj_sen_solver; tspan = (0.0, 0.01), dt = 1e-3)
```

### 步骤 5：直接灵敏度求解

直接法复用同一个 `dyn_solver` 和 `objt`，只需更换外层求解器类型：

```@example
drc_sen_solver = DirectDynamicsSensitivitySolver(dyn_solver, objt)
drc_sim = solve!(prob, drc_sen_solver; tspan = (0.0, 0.01), dt = 1e-3)
```

## 提取梯度与雅可比

下面的示例直接从 `sim.solver_cache` 读取结果。

### 伴随法结果

```@example
∂J∂x0 = adj_sim.solver_cache.∂J∂x₀ᵀ
∂J∂θ = adj_sim.solver_cache.∂J∂θᵀ
∂J∂c = adj_sim.solver_cache.∂J∂cᵀ

∂J∂x0, ∂J∂θ, ∂J∂c
```

`∂J∂x0` 表示目标函数对初始状态的梯度，`∂J∂θ` 与 `∂J∂c` 则分别汇集了对控制参数和结构参数的灵敏度轨迹。

### 直接法结果

直接法会把状态雅可比直接存放在 `drc_sim.solver_cache` 中：

```@example
Jac_state = drc_sim.solver_cache.Jac_state
Jac_action = drc_sim.solver_cache.Jac_action
Jac_control_params = drc_sim.solver_cache.Jac_control_params

state_sample = Jac_state[end][1:3, 1:3]
action_sample = Jac_action[end][1:3, :]
control_param_total = sum(sum(abs, block) for block in Jac_control_params)

(; state_sample, action_sample, control_param_total)
```

这些雅可比分别描述状态对初始状态、驱动量和控制参数的灵敏度。这里 `control_param_total` 仍为零，是因为这里未引入参数化控制策略。
