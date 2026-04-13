# 控制枢纽

```@meta
CurrentModule = Rible
ShareDefaultModule = true
```

如果说 `Structure` 描述了机器人的由肌肉骨骼组成的“身体”，那么`ControlHub`则定义了它的“神经中枢”：它规定了机器人能感知什么以及能执行什么。

1. **通过量器（Gauges）感知**：`Gauges` 定义如何从复杂的物理状态（坐标、速度）中提取反馈控制策略或价值函数所需的特征向量。
2. **通过致动器（Actuators）驱动**：`Actuators` 定义如何将驱动指令转化为物理系统中的广义力。

这种架构使得相同的物理模型可以无缝配合不同的控制界面，便于控制算法的实现。

## 2. 致动器 (Actuators)

所有的致动器都继承自抽象类型 `AbstractActuator`。致动器的职责是将控制器输出的标量序列映射为作用在物理系统中的力。

- `ExternalForceActuator`: 最常用的具体实现，用于在物体的特定位点（Locus）施加力或力矩。

致动器的行为通过以下函数定义：

- **`actuate!(bot, policy, state)`**: 顶层入口。它根据当前的策略计算整个机器人的驱动量，并触发所有子致动器。
- **`execute!(structure, actuator, u)`**: 执行单个致动器，将驱动量 `u` 转化为广义力并累加到系统的力向量中。
- **`gen_force_actu_jacobian!(∂F∂u, structure, actuator, u)`**: 计算产生的广义力相对于驱动量 `u` 的雅可比矩阵。这是伴随矩阵计算和轨迹优化的基础。

## 3. 量器 (Gauges)

量器定义了机器人对自身或环境的观测方式。所有量器都继承自 `AbstractGauge`。

### 采集器与误差量器的区分

为了适应反馈控制和优化求解的不同需求，Rible 将量器分为两类：

#### 3.1 采集器 (CaptumGauge)

`CaptumGauge` 用于获取**原始反馈信号**。它通常输入给策略函数（Policy）或用于状态监控。

- **组成**：由一个 `Signifier`（指向物体位点）和一个 `Captum`（如 `PositionCaptum`）组成。
- **输出**：通常是一个向量（如 3D 位置向量）。

#### 3.2 误差量器 (ErrorGauge)

`ErrorGauge` 在观测的基础上引入了**参考值 (Reference)**。它直接计算“现状”与“目标”之间的偏差。

- **数学本质**：通常计算加权平方和误差

```math
e = \frac{1}{2}\|a - a_{ref}\|^2
```

- **输出**：是一个标量。
- **应用**：`ErrorGauge` 是构建 Objective（轨迹优化代价函数）的原子单元。

### 核心接口：测量

量器的计算围绕以下接口展开：

- **`measure!(out, st, sig, cap)`**: 根据具体的标识符（Signifier）和采集器（Captum）类型，将物理属性（如特定位点的位置、速度等）直接写入输出缓冲区 `out`。
- **`measure_jacobian!(Jq, Jv, Js, st, sig, cap)`**: 用于计算测量值相对于系统广义坐标 ``q``、广义速度 ``v`` 以及辅助变量 ``s`` 的雅可比矩阵。

```@setup
using Rible
using StaticArrays
using LinearAlgebra
import TypeSortedCollections as TSC
```

## 4. 示例：构建一个简易倒立摆

以下展示了如何为一个简化的倒立摆（Cart-pole）系统构建 `ControlHub`。

### 步骤 1：物理结构准备

首先，我们需要为机器人创建一个物理结构（`Structure`）。在这个例子中，我们定义了一个带有杆子的小车。

```@example cartpole
using Rible
using StaticArrays
using LinearAlgebra
import TypeSortedCollections as TSC

# 物理结构准备 (小车 L=1.0)
l = 1.0
mass_locus = Locus(@SVector zeros(3))
loci = [Locus(@SVector zeros(3)), Locus(SVector(0.0, l, 0.0))]
prop_cart = RigidBodyProperty(1, true, 1.0, SMatrix{3,3,Float64}(I), mass_locus, loci)
prop_pole = RigidBodyProperty(2, true, 0.5, SMatrix{3,3,Float64}(I), mass_locus, loci)

cart = RigidBody(prop_cart, RigidBodyState(prop_cart, zeros(3), SMatrix{3,3,Float64}(I), zeros(3), zeros(3)), NC1P3V(zeros(3), zeros(3), SMatrix{3,3,Float64}(I)), nothing)
pole = RigidBody(prop_pole, RigidBodyState(prop_pole, SVector(0.0, l/2, 0.0), SMatrix{3,3,Float64}(I), zeros(3), zeros(3)), NC1P3V(SVector(0.0, l/2, 0.0), SVector(0.0, l/2, 0.0), SMatrix{3,3,Float64}(I)), nothing)

bodies = TSC.TypeSortedCollection([cart, pole])
apparatuses = TSC.TypeSortedCollection(Int[])
st_obj = Structure(bodies, apparatuses, Connectivity(bodies, apparatuses))
```

### 步骤 2：定义量器

接下来，我们定义两种类型的量器：一个用于观测小车的状态，另一个用于测量杆子的位置误差。

```@example cartpole
# 定义量器：观测 (Captum) 与 误差 (Error)
cart_captum = CaptumGauge(1, Signifier(cart, 1), PosVelCaptum())
pole_pos_error = ErrorGauge(1, Signifier(pole, 2), PositionCaptum(), [0.0, l, 0.0])

capta_gauges = TSC.TypeSortedCollection([cart_captum])
error_gauges = TSC.TypeSortedCollection([pole_pos_error])
```

### 步骤 3：定义致动器

致动器定义了驱动量如何作用于物理系统。这里我们定义了一个在小车上施加水平推力的致动器。

```@example cartpole
# 定义致动器：水平推力
force_actuator = ExternalForceActuator(
    1, Signifier(cart, 1), NaiveOperator(1), [0; 1.0; 0;;], [0.0]
)
actuators = TSC.TypeSortedCollection([force_actuator])
```

### 步骤 4：组装 Hub 并进行测量

最后，我们将所有组件组装成一个 `ControlHub`，并演示如何执行测量。

```@example cartpole
# 组装 Hub 并进行测量
coalition = Coalition(st_obj, capta_gauges, error_gauges, actuators)
hub = ControlHub(st_obj, capta_gauges, error_gauges, actuators, coalition)

# 执行测量
c = measure(st_obj, cart_captum)
e = measure(st_obj, pole_pos_error)
(c, e)
```
