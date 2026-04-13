# RibleTensegrity

```@meta
ShareDefaultModule = true
```

```@setup tensegrity
using Rible
import Rible as RB
using RibleTensegrity
import RibleTensegrity as RT
using StaticArrays
using LinearAlgebra
using SparseArrays
using TypeSortedCollections
using CircularArrays
import CircularArrays as CA
using Rotations
using Unitful
using ElasticArrays
using DataStructures
using EponymTuples
using StructArrays
using FileIO, MeshIO
import GeometryBasics as GB
using LoggingExtras
using Polyhedra
import CDDLib
using Meshes
using Makie
import CairoMakie

CairoMakie.activate!()
Makie.inline!(true)
global_logger(ConsoleLogger(stdout, Logging.Warn; show_limited=false))

# Body definitions (from main Rible repo)
include(joinpath(pathof(Rible), "../../examples/bodies/rigidbar.jl"))
include(joinpath(pathof(Rible), "../../examples/bodies/rigidbar_nonsmooth_repro.jl"))
include(joinpath(pathof(Rible), "../../examples/bodies/make_3d_bar.jl"))
include(joinpath(pathof(Rible), "../../examples/bodies/make_3d_plate.jl"))
include(joinpath(pathof(Rible), "../../examples/bodies/make_3d_tri.jl"))
include(joinpath(pathof(Rible), "../../examples/bodies/new_deck.jl"))

# Robot definitions
include(joinpath(pathof(Rible), "../../examples/robots/superball.jl"))
include(joinpath(pathof(RibleTensegrity), "../../examples/robots/bridge3d.jl"))

const nofield = RB.NoField()
const gravity = RB.Gravity()
nothing
```

**RibleTensegrity** 是专门用于张拉整体（Tensegrity）结构建模与分析的扩展包 [luoUnifiedApproachDynamic2024a](@cite)。张拉整体结构起源于由 Kenneth Snelson 和 Buckminster Fuller 等人创造的建筑艺术，其核心特征是由张力维持的整体结构。

在机器人领域，张拉整体结构因其极高的强度重量比、可变形性以及天然的抗碰撞冲击性，被视为**刚柔耦合机器人**的理想框架。本软件包提供了从建模设计、逆静力找形到动力学仿真与刚度优化等工具链。

## 核心类型

### TensegrityStructure

`TensegrityStructure` 是管理所有刚体和张拉索的结构容器。它继承自 `AbstractStructure`，持有系统状态、连接性和接触相关数据：

```julia
struct TensegrityStructure{BodyType,TenType,CntType,StateType,CRType,CacheType} <:
    AbstractStructure{BodyType,TenType,CntType}
```

构造 `TensegrityStructure` 需要三个参数：刚体、装置和连接性对象：

```julia
st = RT.TensegrityStructure(bodies, apparatuses, connectivity)
```

与一般`Structure`不同，`TensegrityStructure` 在处理力组装时会特别考虑索的**单边受力**特性：索只能拉、不能推。

### CableJoint

`CableJoint` 是实现"拉索"的关键类型。与普通铰链不同，它**不引入**运动学硬约束，而仅定义力作用的方向向量：

```julia
struct CableJoint{hen2eggType} <: AbstractJoint
```

每个 `CableJoint` 通过 `Hen2Egg` 配对连接两个刚体上的一对节点。控制索力的单边弹性规律为：

```math
f_{\text{tension}} = \max(0, k(l - l_0) + c \dot{l})
```

其中 ``k`` 为刚度，``l`` 为当前长度，``l_0`` 为原长，``c`` 为阻尼系数。这体现了索的松弛特性（Slackness）：当索短于原长时，力为零。

## 结构构建

本包提供辅助函数来简化构建流程。典型工作流程为：

1. **定义刚体**（如杆件、板件）
2. **指定连接矩阵**，指明哪些刚体的哪些节点由索连接
3. **格式转换**和**创建索**：使用 `change_connecting_format` + `connect_spring`
4. **构造** `TensegrityStructure`

### 示例：Superball（6 杆张拉整体）

Superball 是经典的 6 杆张拉整体机器人，包含 24 根索。内置的 `superball()` 构造函数演示了完整的构建流程：

```@example tensegrity
# 使用默认参数构造 Superball
bot = superball(0.0;
    l = 1.7 / 2,        # 半杆长
    d = 1.7 / 4,        # 半杆距
    z0 = 2.0,
    visible = false,
    loadmesh = false,
    constrained = false,
)

st = bot.structure

println("刚体数量: $(length(st.bodies))")
println("索数量: $(length(RT.get_cables(st)))")
println("全局坐标数: $(st.connectivity.num_of_full_coords)")
```

`superball()` 内部遵循以下模式：

```julia
# 1. 使用 rigidbar 辅助函数定义刚体
rbs = [rigidbar(i, p₁, p₂; m=5.0, μ, e) for i in 1:6]
rigidbodies = TypeSortedCollection(rbs)

# 2. 为每根索定义弹簧-阻尼器
spring_dampers = [RB.DistanceSpringDamper3D(restlen, k, c; slack=false) for i in 1:24]

# 3. 构建连接矩阵并转换格式
#    格式：每行为 [bid₁, pid₁, bid₂, pid₂]
#    正值 = Hen 端点，负值 = Egg 端点
cm = RT.change_connecting_format(rigidbodies, connecting_matrix)

# 4. 创建 CableJoint + DistanceSpringDamper 对
cables = RT.connect_spring(rigidbodies, spring_dampers; cm, istart=0)

# 5. 选择连接类型并组装
cnt = RB.Connectivity(rigidbodies, cables)           # 自由结构
# 或：cnt = RB.PresFreeConnectivity(rigidbodies, cables)  # 预设-自由分解
st = RT.TensegrityStructure(rigidbodies, cables, cnt)
```

### 检视 API

构建完成后，可使用以下检视函数检查索的属性：

```@example tensegrity
# 当前构型下的索张力、长度、刚度
tensions  = RT.get_cables_tension(st)
lengths   = RT.get_cables_len(st)
stiffness = RT.get_cables_stiffness(st)
restlens  = RT.get_cables_restlen(st)

println("张力范围: [$(minimum(tensions)), $(maximum(tensions))]")
println("长度范围: [$(minimum(lengths)), $(maximum(lengths))]")
```

其他检视函数包括 `get_cables_deform`（变形量）、`get_cables_len_dot`（长度变化率）、`get_cables_force_density`（力密度）。

## 逆静力分析与找形

给定目标构型 ``q``，系统的静力平衡条件为受拉索力向量与外力的平衡：

```math
B \gamma = \tilde{F}_{ext}
```

其中 ``B \in \mathbb{R}^{n \times m}`` 是平衡矩阵，映射 ``m`` 根索的力密度 ``\gamma`` 到 ``n`` 个系统自由度上的合力。``\tilde{F}_{ext}`` 是当前位姿下的重力等广义外力。

### 逆向求解原长

给定目标构型，`inverse_for_restlength` 求解使系统在该构型下达到静力平衡的索原长：

```julia
μ = RT.inverse_for_restlength(bot, botref, field; fmin=0.0, eps_rel=1e-6, verbose=false)
```

- `bot`：目标构型下的机器人/结构
- `botref`：参考结构（通常与 `bot` 相同）
- `field`：外力场（如 `RB.NoField()`、`RB.Gravity()`）
- `fmin`：最小预张力，确保索不会松弛

对于静定系统，直接求解。对于超静定或存在单边约束（``γ ≥ 0``）的系统，自动使用 COSMO 进行二次规划求解：

```math
\min_{\mu} \frac{1}{2} \mu^T \text{diag}(\kappa)\,\mu + h^T\mu \quad \text{s.t.} \quad B\mu = \tilde{F}_{\text{ext}},\;\mu \geq 0,\;\kappa \cdot (\ell - \mu) \geq f_{\min}
```

- **`inverse_for_actuation(bot, target_q, field)`**：若结构包含驱动器（如卷扬机或电机），求解实现目标态所需的驱动变量序列。

### 示例：张拉整体桥梁

```@example tensegrity
# 构建一个 2 模块张拉整体桥梁
bridge = bridge3d(; n=2)

# 求解满足最小预张力的原长
μ = RT.inverse_for_restlength(bridge, bridge, nofield; fmin=1e4, eps_rel=1e-10)

# 应用求解得到的原长
RT.set_restlen!(bridge, μ)
RB.update!(bridge.structure, nofield)

# 验证张力范围
tension_min, tension_max = RT.get_cables_tension(bridge.structure) |> extrema
println("张力范围: [$tension_min, $tension_max]")
```

### GDR：广义动力学松弛

对于复杂结构，当解析逆求解困难时，可使用核心 `Rible` 包提供的 `GDR!` 进行数值松弛以寻找平衡。它引入虚构动能耗散，驱动结构趋向最低能量平衡点：

```@example tensegrity
# 在重力场下运行 GDR，寻找重力平衡构型
bridge_gdr = deepcopy(bridge)
RB.GDR!(bridge_gdr, gravity; β=1e-4, res=1e-9)

# 提取平衡构型
q_eq = bridge_gdr.traj.q[end]
RB.set_initial!(bridge_gdr, q_eq, zero(q_eq))

# 验证静力平衡
isequilibrium, _ = RT.check_static_equilibrium_output_multipliers(
    bridge_gdr.structure, gravity
)
println("达到平衡: $isequilibrium")
```

### 特征值分析

达到平衡后，可计算结构的固有频率：

```@example tensegrity
eigenvals, _ = RB.undamped_eigen(bridge_gdr.structure, gravity)
println("固有频率（前 5 阶）: $(eigenvals[1:min(5, length(eigenvals))])")
```

## 动力学仿真

张拉整体动力学由标准的 Rible 求解器管线处理。关键考虑因素是接触建模——张拉整体机器人通常在地面上滚动或与之交互。

完整的动力学仿真示例（含 superball 滚动仿真）请参阅[入门指南](@ref "Rible.jl 入门")。

## 刚度与稳定性分析

### 静力-运动学确定性

`static_kinematic_determine` 对转置平衡矩阵进行 SVD 分解，将结构分解为**自应力模态**和**机构位移模态**：

```@example tensegrity
# 构建用于刚度分析的 Superball
ballbot_s = superball(0.0;
    θ = 0.0,
    l = 2.0 / 2,
    d = 2.0 / 4,
    z0 = 2.0 / 2,
    visible = false,
    loadmesh = false,
    constrained = false,
)
st_s = ballbot_s.structure

# 验证静力平衡并获取索张力
RT.check_static_equilibrium_output_multipliers(st_s, nofield)
RB.update!(st_s)
f = RT.get_cables_tension(ballbot_s)

# 构建平衡矩阵 Bᵀ
Q̃ = RT.build_Q(st_s)
L̂ = RT.build_L̂(st_s)
Bᵀ = -Q̃ * L̂

# 投影到自由度零空间
Ǎ = RB.cstr_jacobian(st_s, st_s.state.system)
Ň = RB.nullspace(Ǎ)
ℬᵀ = transpose(Ň) * Bᵀ

# 分解为自应力模态和机构位移模态
S, D = RT.static_kinematic_determine(ℬᵀ)
println("自应力模态数: $(size(S, 2))")
println("机构位移模态数: $(size(D, 2))")
```

- **自应力模态**（`S` 的列）：不依赖外力、仅靠内部预拉力维持平衡的内力向量。
- **机构位移模态**（`D` 的列）：结构在当前预应力水平下可能失稳的运动方向（无穷小机构）。

### 材料刚度与几何刚度

切线刚度矩阵 ``K_T`` 是判断结构稳定性的关键：

```math
K_T = K_{\text{material}} + K_{\text{geometric}} + K_{\text{constraint}}
```

材料刚度 ``K_{\text{material}}`` 反映了缆索受拉伸时的弹性阻力。几何刚度 ``K_{\text{geometric}}`` 与约束刚度 ``K_{\text{constraint}}`` 均源于几何扰动，且受缆索预应力水平的影响；然而，研究发现仅 ``K_{\text{constraint}}`` 是唯一导致失稳的根源。

```@example tensegrity
q = RB.get_coords(st_s)
k = RT.get_cables_stiffness(st_s)

# 材料刚度矩阵
Ǩm = RT.build_material_stiffness_matrix!(st_s, q, k)

# 几何刚度矩阵
Ǩg = RT.build_geometric_stiffness_matrix!(st_s, q, f)

# 投影到自由坐标空间
𝒦m = transpose(Ň) * Ǩm * Ň |> Symmetric
𝒦g = transpose(Ň) * Ǩg * Ň |> Symmetric

vals_m = sort(eigvals(𝒦m))
println("材料刚度最小特征值: $(vals_m[1])")
```

结构稳定的充分条件是总切线刚度 ``K_T`` 在机构位移方向上为正定。

## 预应力优化

借助凸优化求解器（COSMO、Clarabel），RibleTensegrity 可在自应力空间内搜索最优性能点 [luoStabilityConditionsStiffness2024](@cite)：

- **`optimize_maximum_stiffness`**：最大化刚度矩阵的最小特征值，提升结构抵御外载荷的能力。
- **`optimize_zero_stiffness`**：寻找最小特征值过零的临界预应力，识别失稳起点。

### 示例：最大化刚度

```@example tensegrity
# 为每个自应力模态构建预应力刚度贡献
ns = size(S, 2)
vec𝒦ps = map(1:ns) do i
    si = S[:, i]
    λi = inv(Ǎ * transpose(Ǎ)) * Ǎ * Bᵀ * si
    Ǩai = -RB.cstr_forces_jacobian(st_s, q, λi)
    Ǩgi = RT.build_geometric_stiffness_matrix!(st_s, q, si)
    𝒦pi = transpose(Ň) * (Ǩgi .+ Ǩai) * Ň |> Symmetric
    vec(𝒦pi)
end
mat𝒦ps = reduce(hcat, vec𝒦ps)

# 设置优化问题
vec𝒦m = vec(𝒦m)
vecI = vec(Matrix(1.0I, size(𝒦m)))
ᾱ = ones(ns)                         # 自应力组合权重
nx = ns + 2                           # 决策变量：[α..., σ, ρ]

A = hcat(-Matrix(1.0I, ns, ns), ᾱ, zero(ᾱ))
b = zeros(ns)

# 最大化最小刚度 (ρ)
result = RT.optimize_maximum_stiffness(mat𝒦ps, vec𝒦m, vecI, A, b, nx)
σ_opt = result.x[end-1]               # 最优预应力缩放
ρ_opt = result.x[end]                 # 最大最小特征值

println("最优预应力缩放: $σ_opt")
println("最大最小特征值: $ρ_opt")
```

### 寻找零刚度点

```@example tensegrity
# 搜索最小刚度趋于零的临界预应力
result_zero = RT.optimize_zero_stiffness(
    mat𝒦ps, vec𝒦m, vecI,
    hcat(-Matrix(1.0I, ns, ns), ᾱ),   # 等式约束
    zeros(ns),                          # 约束右端项
    ns + 1,                             # 决策变量数
    result.x[1:end-1],                  # 用最大刚度结果热启动
)

σ_zero = result_zero.x[end]
println("零刚度预应力缩放: $σ_zero")

# 验证：该缩放下刚度矩阵应有近零特征值
𝒦_zero = 𝒦m + σ_zero * reshape(mat𝒦ps * ᾱ, size(𝒦m))
ρ_zero = minimum(eigvals(Symmetric(𝒦_zero)))
println("零刚度点最小特征值: $ρ_zero")
```

!!! tip "求解器与热启动"
    `optimize_maximum_stiffness` 和 `optimize_zero_stiffness` 内部使用 COSMO。零刚度搜索接受热启动向量（`x_0`）——用 `optimize_maximum_stiffness` 的解作为起点可显著提升收敛速度。

## 参考文献

```@bibliography
Pages = ["RibleTensegrity.md"]
```
