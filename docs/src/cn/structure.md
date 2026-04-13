# 结构

```@meta
CurrentModule = Rible
ShareDefaultModule = true
```

```@setup
using Rible
import Rible as RB
using StaticArrays
using LinearAlgebra
```

`AbstractStructure`/`Structure` 是管理多体系统所有力学属性、拓扑关系及运动状态的核心容器。它通过整合物体（body）和器件（apparatus），利用半自动生成的 `Connectivity` 进行编号与索引，并维护 `StructureState` 和 `StructureCache` 对外提供内部动力学接口。

本项目核心采用基于泛型函数与多重分派的接口规范。围绕 `AbstractStructure` 定义的系列函数（如 `update!`、`cstr_function!`）确立了行为标准，开发者只需为新子类型实现对应方法，即可无缝接入仿真框架。

- `Structure` 是标准具体实现，适用于大多数涉及刚体及复杂连接的场景。
- 用户可通过继承 `AbstractStructure` 并实现特定方法来扩展功能。例如，`RibleTensegrity` 中的 `TensegrityStructure` 针对张拉整体结构提供了专门的分派。

## 结构的组成

一个典型的结构由三大部分组成：

### 物体与器件

系统中所有运动单元的集合存储在 `bodies` 字段中。由于不同物体的数学描述可能不同，它们通常是异构的集合。系统中所有连接元件（如弹簧、阻尼器、约束关节等）则存储在 `apparatuses` 字段中。

为了在处理异构集合时保持类型稳定并获得高性能，Rible 使用了 `TypeSortedCollection` 进行存储。这确保了在执行大规模系统更新时，相同类型的物体或器件能够触发特定方法的特化，从而避免动态分派的开销。

### 结构状态

`StructureState` 作为系统的规范状态，通过多维度视角管理变量：

- **系统视角**: 在 `system` 字段中以大数组形式存储全局坐标 `q`、速度 `q̇`、加速度 `q̈`、约束力 `λ` 等。
- **成员视角**: 利用 `view` 机制将全局大数组的切片映射给各个物体（`members`），使组件在更新时能直接操作全局内存，消除重复数据的开销。

### 坐标状态

`AbstractCoordinatesState` 是所有坐标级状态容器的抽象超类型，供动力学求解器使用。其标准具体实现为 `CoordinatesState`，存储以下字段：

- **`q`**、**`q̇`**、**`q̈`** — 广义坐标、速度和加速度的平坦向量。
- **`p`** — 共轭动量。
- **`F`** — 由物体力、器件力和外场组装而成的广义力向量。
- **`λ`** — 拉格朗日乘子（约束力）。
- **`s`** — 器件内部状态的辅助变量。
- **`c`** — 结构参数。

`StructureState.system` 即为 一个`CoordinatesState`，求解器在调用结构级函数（如 `cstr_function!`、`assemble_forces!`）时直接操作 `CoordinatesState`。这种分离使求解器处理轻量状态对象，而结构在内部维护完整的成员级映射。

## 连接关系与索引

`Connectivity` 通常由系统根据物体和器件的拓扑布局自动生成。它负责管理复杂的索引映射逻辑。

### 自动化生成

通过调用 `Connectivity(bodies, apparatuses)`，Rible 会遍历所有组件，自动推导出系统的全局自由度、索引偏移以及连接拓扑。

### 内部索引映射

`Connectivity` 内部维护了关键的索引信息：

- **坐标索引**: 精确管理物体坐标在系统坐标向量中的位置。
- **约束索引**: 区分物体的固有约束（如自然坐标下的物体约束）和器件引入的外在约束（如运动副）。
- **位置映射**: 建立物体作用点（Loci）与系统全局状态之间的关联。

## 典型操作与核心功能

结构提供了一系列用于动力学仿真的核心函数。

### 更新管线

`update!` 是完整结构更新的顶层入口，按顺序执行以下步骤：

```julia
function update!(st::AbstractStructure, field=NoField(); )
    clear_forces!(st)
    stretch!(st)
    update_bodies!(st)
    update_apparatuses!(st)
    apply_field!(st, field)
    assemble_forces!(st)
end
```

每一步各有其职责：

- **`clear_forces!(st)`** — 将系统级状态和各个物体缓存中的所有力向量清零。
- **`stretch!(st)`** — 将结构参数值传播到各个物体的位点和各个器件，为后续计算做准备。
- **`update_bodies!(st)`** — 从系统级坐标填充各物体的局部状态，重新计算惯量缓存、变换和位点状态。
- **`update_apparatuses!(st)`** — 利用更新后的物体状态计算所有器件（关节、力元）。
- **`apply_field!(st, field)`** — 对每个物体施加外场（如重力）。
- **`assemble_forces!(st)`** — 将所有物体级和器件级的力汇总到系统级力向量 `F` 中。

较轻量的变体 `lazy_update!` 用 `lazy_update_bodies!` 替代 `update_bodies!`，跳过惯量缓存的重新计算。这在不需要质量矩阵时（如迭代约束求解过程中）很有用。

### 约束与雅可比

- **`cstr_function!(Φ, st, inst_state)`**: 计算系统中所有活动约束的残差向量。
- **`cstr_jacobian!(A, st, inst_state)`**: 计算系统约束对广义坐标的雅可比矩阵。

### 质量矩阵与能量

- **`assemble_M(st)`**: 组装并返回系统的全局质量矩阵。
- **`kinetic_energy(st)`**: 基于当前状态分派计算系统的总动能。

通过这些高度抽象的泛型接口，动力学求解器能够以统一的逻辑处理复杂的物理系统。

## 示例：最小可运行结构

以下示例仅使用包内 API 构建一个最小系统。

### 步骤 1：定义一个刚体

```@example
using Rible
using StaticArrays
using LinearAlgebra

mass_locus = Locus(SVector(0.0, 0.0, 0.0))
loci = [Locus(SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0))]
inertia = SMatrix{3,3,Float64}(I)

prop = RigidBodyProperty(1, true, 1.0, inertia, mass_locus, loci; visible = false)
r0 = SVector(0.0, 0.0, 0.0)
R0 = SMatrix{3,3,Float64}(I)
state = RigidBodyState(prop, r0, R0, zero(r0), zero(r0))
coords = NCF.NC1P3V(r0, r0, R0)
body = RigidBody(prop, state, coords, nothing)
```

### 步骤 2：组装 `Connectivity` 与 `Structure`

```@example
bodies = [body]
apparatuses = Int[]
cnt = Connectivity(bodies, apparatuses)
st = Structure(bodies, apparatuses, cnt)
inst = st.state.system
(; n_bodies = length(st.bodies), n_apparatuses = length(st.apparatuses))
```

### 步骤 3：访问系统状态向量

```@example
q = inst.q
qdot = inst.q̇
(; q_len = length(q), qdot_len = length(qdot))
```

### 步骤 4：执行一次结构更新并读取规模信息

```@example
update!(st)
n_q = get_num_of_coords(st)
n_c = get_num_of_cstr(st)
(; n_q, n_c, dof = n_q - n_c)
```

### 步骤 5：计算约束残差与雅可比

```@example
Φ = zeros(n_c)
cstr_function!(Φ, st, inst)
A = cstr_jacobian(st, inst)
(; residual_norm = norm(Φ), jacobian_size = size(A))
```

### 步骤 6：组装质量矩阵

```@example
M = assemble_M(st)
size(M)
```

### 步骤 7：构建 `Robot` 并进行短时积分

```@example
bot = Robot(st)
prob = DynamicsProblem(bot; env=GravityEnv())
solver = DynamicsSolver(Zhong06())
solve!(prob, solver; dt = 1e-3, tspan = (0.0, 0.01))
length(bot.traj)
```

### 步骤 8：后处理 — 机械能与中点速度

```@example
E_total = mechanical_energy!(bot).E[end]
v_mid = first(get_mid_velocity!(bot, 1, 1))
(; E_total, v_mid)
```
