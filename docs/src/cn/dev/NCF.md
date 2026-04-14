# Natural Coordinates 自然坐标系统

`Natural Coordinates` 是 Rible 自带的用于刚体建模的坐标系统。子模块 `NCF` 位于 `src/coordinates/natural_coordinates`，并通过多重分派接入质量矩阵、约束方程和关节求解流程。

> **核心特点**：该方法直接使用点坐标与方向向量描述姿态。它避免了姿态奇异，并且能得到常数质量矩阵。

## 基本表示

在自然坐标框架中，刚体状态由参考点和方向向量构成：

- 2D 刚杆常用表示：1 个基点 + 1 个基向量
- 2D 刚体常用表示：1 个基点 + 2 个基向量
- 3D 刚体常用表示：1 个基点 + 3 个基向量

为满足特殊建模需要，其他变体也完全支持，如2个基点+2个基向量。

对应到实现层，`NC` 是统一坐标结构，具体维度与自由度由构造函数决定。

## 模块结构与代码入口

`src/coordinates/natural_coordinates` 的主要文件如下：

- `types.jl`：定义 `NC` 和底层数据结构。
- `constructors.jl`：定义常用构造函数，如 `NC2D1P1V`、`NC1P2V`、`NC1P3V`。
- `functions.jl`：坐标映射与雅可比计算，如 `to_position`、`to_position_jacobian`。
- `constraints.jl`：内禀约束与雅可比实现，核心接口为 `cstr_function!`、`cstr_jacobian!`。
- `mass_matrix.jl`：质量矩阵构造，如 `make_M`。
- `joints.jl`：关节缓存、约束违约量与雅可比装配。

## 常见用法

### 1 构造坐标对象

自然坐标构造函数至少需要有一个基本点 `ri`，也可传其他基本点和基向量如 `u`，所有这些都在全局坐标系中表达。

```julia
# 2D 刚杆：1 点 + 1 向量
nc_bar_2d = NC2D1P1V(ri, u)

# 2D 刚体：1 点 + 2 向量
nc_body_2d = NC1P2V(ri)

# 3D 刚体：1 点 + 3 向量
nc_body_3d = NC1P3V(ri)
```

### 2 位置映射与雅可比

```julia
# 局部点 c 在全局下的位置
pos = to_position(nc_body_3d, q, c)

# 对坐标的雅可比
Jpos = to_position_jacobian(nc_body_3d, q, c)
```

### 3 约束值与约束雅可比

`constraints.jl` 的核心接口是就地更新版本：

```julia
ret_c = zeros(nc_num_constraints)
ret_J = zeros(nc_num_constraints, length(q))

cstr_function!(ret_c, nc_body_3d, q)
cstr_jacobian!(ret_J, nc_body_3d, similar(ret_J), q)
```

### 4 质量矩阵

```julia
M = make_M(nc_body_3d, mass, inertia, mass_center)
```

## 关节开发指南

本节面向需要新增关节类型的开发者。Rible 的关节实现采用“基础约束组合”策略：先定义基本约束类型，再在 `src/joint/joints.jl` 的 `get_joint_info` 中通过 `mask` 组合构造具体关节。

### 理论基础 基本约束

在关节模块中，所有约束都先写成统一的“联合坐标”形式。设

- ``q_h``：母体（hen body）坐标
- ``q_e``：子体（egg body）坐标
- ``q_j = [q_h; q_e]``：拼接后的联合坐标

`build_joint_cache` 的核心目标，是把几何约束预编译成两类矩阵对象：

- `transformations`：一次项矩阵（线性项）
- `halves`：二次型矩阵列表（每个约束对应一个 ``H_i``）

因此单个约束都可写为：

```math
C_i(q_j) = q_j^T H_i q_j + T_i q_j - c_i
```

其中 ``T_i`` 是 `transformations` 的第 ``i`` 行，``H_i`` 是 `halves[i]`。`violations` 保存参考构型下的常量偏置 ``c_i``，在运行时通过 `ret .-= violations` 归一到“零违约”基线。

> **线性 vs 二次 的具体含义**：
> - **线性约束**：``C_i(q_j)`` 对 ``q_j`` 是一次函数，雅可比是常量。
> - **二次约束**：``C_i(q_j)`` 对 ``q_j`` 是二次多项式，雅可比随 ``q_j`` 线性变化，但海森矩阵是常量。

### Mask 1st 重合约束 Type 1

约束两点重合：

```math
r_i - r_j = 0
```

在实现中，这类约束主要来自点位映射雅可比拼接矩阵 ``J=[-C_h\ C_e]`` 的线性部分：

- **代码索引**：`mask_1st`
- **矩阵贡献**：进入 `transformations = J[mask_1st, :]`
- **约束阶次**：线性（对应 ``H_i=0``，仅保留 ``T_i q_j``）

### Mask 2nd 距离约束 Type 4

约束两点距离为常数 ``d``：

```math
\|r_i - r_j\|^2 - d^2 = 0
```

这类约束在代码里由二次型构建：

- **代码索引**：`mask_2nd`
- **矩阵贡献**：`half_2nd[1] = J' * J`
- **约束阶次**：二次（仅有 ``q_j^T H_i q_j``，常量项由 `violations` 提供）

### Mask 3rd 正交约束 I Type 3

约束方向向量与相对位移正交：

```math
v_h^T (r_i - r_j) = 0
```

其二次型矩阵由方向向量选择矩阵与 Kronecker 积组合得到：

- **代码索引**：`mask_3rd`
- **矩阵构建**：通过 `kron(...)` 生成 `half_3rd`
- **约束阶次**：二次
- **常见用途**：平面内运动约束、滑块方向约束

### Mask 4th 正交约束 II Type 2

约束两个方向向量正交：

```math
v_h^T v_l = 0
```

这类约束完全由方向坐标构成，属于纯姿态耦合项：

- **代码索引**：`mask_4th`
- **矩阵构建**：通过方向基与 `kron(...)` 组合生成 `half_4th`
- **约束阶次**：二次
- **常见用途**：转轴对齐、万向副约束组合

## 实现细节与性能

### 缓存机制 ApparatusCache

`build_joint_cache` 在初始化阶段返回 `cache` 与 `violations`。其中 `cache` 内含运行时所需的全部稀疏结构：

- `transformations`：线性约束矩阵
- `halves`：二次型矩阵列表
- `hessians`：每个约束的对称海森矩阵 ``H_i + H_i^T``
- `joint_q`、`joint_work`、`trf_work`：就地计算缓冲，避免分配

这种设计让运行阶段只做数值代入，不再重建结构矩阵。

### joints.jl 关键函数

- `build_joint_cache`：根据 mask 组合预计算 `transformations`、`halves`、`hessians`。
- `get_joint_violations!`：先算 ``q_j^T H_i q_j``，再加线性项 ``T_i q_j``，最后减去 `violations`。
- `get_joint_jacobian!`：使用预计算的 `hessians` 与当前 `q_j` 装配 `C_q`。

## 扩展新关节

扩展新关节时，通常不需要新增底层数学核函数。推荐流程：

1. 在 `src/joint/joints.jl` 的 `get_joint_info` 中定义目标关节的 mask 组合。
2. 用现有 `mask_1st` 到 `mask_4th` 组合目标自由度约束。
3. 复用 `build_joint_cache`、`get_joint_violations!`、`get_joint_jacobian!` 完成求解接入。

这样可以保持新关节与现有求解流程完全兼容，并继承缓存机制带来的性能优势。
