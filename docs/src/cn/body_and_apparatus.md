# 体与器

```@meta
CurrentModule = Rible
ShareDefaultModule = true
```

```@setup
using Rible
import Rible as RB
using StaticArrays
using LinearAlgebra
using FileIO, MeshIO
using Makie
import CairoMakie
CairoMakie.activate!()
Makie.inline!(true) # Make sure to inline plots into Documenter output!
```

## 位点（Locus）

!!! note "The place where something is situated or occurs. — Cambridge Dictionary"

在多体动力学中，准确描述物体上各个物质点的位置和方位是建立连接、施加外力和处理接触碰撞的基础。此外，在处理接触和碰撞时，这些点表面的法线方向、摩擦系数和回弹系数等物理属性也至关重要。

为了统一且高效地管理这些信息，Rible 引入了**位点（`Locus`）**和**位点状态（`LocusState`）**类型。其中`Locus` 是一个immutable类型，用于在物体的**局部坐标系**中定义一个物质点的不变属性：

- `position`：该点在局部坐标系中的位置向量。
- `axes`：该点处的局部参考轴系（`Axes`），通常用于定义接触面的法线方向或关节的旋转/平移轴。默认情况下，第一根轴（x轴）即为接触法向。
- `friction_coefficient` 和 `restitution_coefficient`：分别定义该点处的摩擦系数和碰撞回弹系数。

### 构造位点

在构建刚体时，我们需要提供质心的 `mass_locus` 以及一系列用于后续连接或接触的位点 `loci` 向量。

```@example
# 定义局部坐标系下的位点 (Locus)
mass_locus = Locus(SVector(-0.2, 0.0, 0.0)) # 质心坐标，默认无摩擦/回弹

# 定义两个位点。由于构造函数的灵活分派，我们可以按需提供参数
loci = [
    # 位点1：指定位置和参考轴法线
    Locus(SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0)),
    # 位点2：指定位置、参考轴法线、摩擦系数(0.5)和回弹系数(0.1)
    Locus(SVector(0.5, 0.0, 0.0), SVector(1.0, 0.0, 0.0), 0.5, 0.1)
]
```

`LocusState` 则是mutable类型，记录了该点在仿真过程中的**全局（世界）坐标系**状态：

- `frame`：该点在全局坐标系下的完整运动学框架（`CartesianFrame`），包括它当前的全局位置、全局姿态（旋转矩阵）、全局线速度和全局角速度。
- `force` 和 `torque`：当前作用在该点上的全局力和力矩。
- `contact_state`：维护当前点的接触状态，例如是否正处于激活的接触中、接触间隙（gap）以及接触力等。

## 体（Body）

!!! note "体用一源，显微无间。 — 程颐《易传序》"

体（Body）是多体系统的基础构建块。在 Rible 中，所有的体都是抽象类型 `AbstractBody{N,T}` 的子类，其中 `N` 代表空间维度（2或3），`T` 代表数值类型（通常为 `Float64`）。这种参数化设计确保了系统既能处理平面问题也能处理空间问题，同时保持类型安全和高性能。

### 刚体 （RigidBody）

刚体（`RigidBody`）是多体系统中最基本的单元之一，在Rible中它由`RigidBodyProperty`、`RigidBodyState`、坐标描述（Coordinates）和用于可视化的 Mesh 共同构成。

#### 刚体属性

`RigidBodyProperty` 定义刚体的不变量。这包括质量属性（质量、惯量）和几何特征（质心位置、位点 Loci）。
一个刚体的物理特性由 `RigidBodyProperty` 定义，它不仅包含动力学参数，还包含仿真所需的元数据。动力学参数主要涉及质量（`mass`）和惯性张量（`inertia`），其中惯性张量通常存储为 `SMatrix` 以优化内存和计算速度。
此外，每个刚体都有唯一的 `id` 和可选的 `type`（符号标记）。
布尔字段 `contactable` 和 `visible` 分别控制刚体是否参与接触检测以及是否在可视化中渲染。

```@example
id = 1
mass = 2.5
inertia = SMatrix{3,3,Float64}(I)  # 单位惯性矩阵

prop = RigidBodyProperty(
    id,
    true,           # contactable
    mass,
    inertia,
    mass_locus,
    loci;
    visible=true
)
```

#### 刚体状态

`RigidBodyState` 维护刚体的瞬时状态。这包含了在全局坐标系中的位置矢量和表示方向的旋转矩阵，以及线速度和角速度。这些状态量是仿真器在每个时间步更新的对象。

```@example
r0 = SVector(0.0, 0.0, 0.0)      # 初始位置
R = SMatrix{3,3}(I)              # 初始姿态 (无旋转)
r_dot = zero(r0)                 # 初始线速度
omega = zero(r0)                 # 初始角速度
state = RigidBodyState(prop, r0, R, r_dot, omega)
```

#### 坐标描述

系统使用**节点坐标法 (NCF)** 来描述运动。根据物体的类型（如自由浮动的块、连接两点的连杆），选择合适的坐标公式。例如，`NC1P3V`（1点3矢量）适用于标准的6自由度刚体，`NC3D2P`（3D 2点）适用于由两个端点定义的连杆，而 `NC2D2P` 则用于平面内的连杆。

```@example
# 对于一般的空间刚体，使用一个基准点和三个方向矢量 (NC1P3V)
coords = NCF.NC1P3V(r0, r0, R)
```

#### 可视化网格

可以为刚体加载可视化用的网格。

```@example
# 加载 STL 文件并进行缩放/平移 (示例)
# mesh = load(assetpath("link.STL")) |> make_patch(; scale=0.1)
mesh = nothing # 此处省略详细网格
```

#### 构造

最后，将上述所有实例组装为刚体。

```@example
# 最终组装实体
body = RigidBody(prop, state, coords, mesh)
```

在内部，会自动构建一个 `RigidBodyCache` 用于存储中间计算缓存。

## 器 （Apparatus）

!!! note "大器免成。    ――《老子》帛书本"

器件（`Apparatus`）由 **铰接（Joint）** 和**力（Force）**组成，是系统中代表物理上的关节或力（如弹簧）的核心组件，也是多体动力学方程中**约束项**和**控制力**的载体。

### 力（Force）

力（Force）负责计算作用在体上的广义力。所有力模型均继承自抽象类 `AbstractForce`。

基础力元包括：

- **距离弹簧阻尼器** (`DistanceSpringDamper`)，这是最基础的力元，连接两个点（Loci），维护着 `DistanceSpringDamperState` 以存储跟踪当前的长度、伸长率和张力。

- **旋转/扭转弹簧** (`Torsional/RotationalSpringDamper`) 则用于产生恢复力矩，通常作用于关节的转轴上，抵抗相对旋转。

### 铰接（Joint）

铰接（Joint）定义了体与体之间的相对连接关系。所有铰接模型均继承自 `AbstractJoint`。

铰接通常（但不总是）会对多体动力学方程引入**约束力**。标准铰接及其变体通过 `ProtoJoint` 类型实现：

| 关节类型 | 自由度 | 说明 |
| :--- | :---: | :--- |
| `:Fixed` | 0 | **固定关节**。将两个刚体完全锁定，或将刚体上的某点锁定在空间固定位置（`FixedPointJoint`）。常用于基座固定。 |
| `:Revolute` | 1 | **转动关节**。允许绕共用轴旋转。需要指定两个物体上的点重合以及轴线共线。 |
| `:Prismatic` | 1 | **移动关节**。允许沿轴线平移，限制旋转。 |
| `:Cylindrical` | 2 | **圆柱关节**。允许绕轴旋转和沿轴滑动。 |
| `:Universal` | 2 | **万向节**。由两个正交的转动轴组成，允许两个方向的旋转。 |
| `:Spherical` | 3 | **球关节**。仅约束平动，允许绕连接点任意旋转。 |

#### 定义连接

构建一个关节的核心在于清楚描述"连接关系"。这通过以下两个类型实现：

**锚点 (`Anchor`)** — 精确定位连接点的几何信息（位置与约束轴）。

- `Anchor(body, pid)` — 从同一个位点读取位置与轴。
- `Anchor(body, position_pid, axes_pid)` — 位置与轴分别来自不同位点（常用于转动关节等需要指定轴向的场景）。

**连接对 (`Hen2Egg`)** — 定义"母体"（Hen）与"子体"（Egg）的配对。

- `Hen2Egg(hen_anchor, egg_anchor)` — 建立两个锚点之间的对应关系。

构建关节时，先用 `Anchor` 明确两端几何，再用 `Hen2Egg` 描述拓扑方向。

#### 铰接构造示例

```@example
prop1 = RigidBodyProperty(1, true, 1.0, inertia, mass_locus, loci)
state1 = RigidBodyState(prop1, r0, R, r_dot, omega)
nmcs1 = NCF.NC1P3V(r0, r0, R)
body1 = RigidBody(prop1, state1, nmcs1, nothing)

prop2 = RigidBodyProperty(2, true, 1.0, inertia, mass_locus, loci)
r0_2 = r0 + SVector(1.0, 0.0, 0.0)
state2 = RigidBodyState(prop2, r0_2, R, r_dot, omega)
nmcs2 = NCF.NC1P3V(r0_2, r0_2, R)
body2 = RigidBody(prop2, state2, nmcs2, nothing)
```

```@example
# 1. 定义连接锚点：连接 body1 的第2个位点和 body2 的第1个位点
hen = Anchor(body1, 2)
egg = Anchor(body2, 1)

# 2. 建立连接关系
conn = Hen2Egg(hen, egg)
```

如果关节需要特定的旋转轴或平移轴（例如转动关节），应使用
`Anchor(body, position_pid, axes_pid)` 分别指定连接位置与约束轴：

```@example
# 指定 body1 上第1个位点作为位置/轴，body2 上第2个位点作为位置/轴
hen_axis = Anchor(body1, 1, 1)
egg_axis = Anchor(body2, 2, 2)
conn_axis = Hen2Egg(hen_axis, egg_axis)
```

创建关节对象时，使用 `ProtoJoint` 函数和关节类型符号。

```@example
# 创建一个球关节
joint_spherical = ProtoJoint(conn, :Spherical)

# 创建一个转动关节
joint_revolute = ProtoJoint(conn_axis, :Revolute)
```
