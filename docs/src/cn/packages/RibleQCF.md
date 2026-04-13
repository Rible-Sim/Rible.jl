# 四元数坐标 — RibleQCF

```@meta
CurrentModule = Rible
ShareDefaultModule = true
```

```@setup qcf
using Rible
import Rible as RB
using RibleQCF
import RibleQCF as QCF
using StaticArrays
using LinearAlgebra
using Rotations
using TypeSortedCollections
using Makie
import CairoMakie
using FileIO
using MeshIO
import GeometryBasics as GB

CairoMakie.activate!()
Makie.inline!(true)
```

RibleQCF 是一个扩展包，在默认的自然坐标格式（NCF）之外提供基于四元数的坐标类型 `QC`。单个刚体由 **7 个广义坐标**描述——3 个平移分量和 4 个单位四元数分量——在单位模长约束被处理之后得到 6 个自由度。其惯量表示及数值影响遵循 [xuNumericalInfluencesInertia2016](@cite) 和 [xuNumericalInfluenceAdditional2020](@cite) 所建立的框架。


!!! tip 
    **QCF** 会得到非常数质量矩阵，需要在状态更新时重新计算，并涉及到复杂的雅可比计算。

## 坐标类型 `QC`

核心类型 `QCF.QC` 存储质量、对角化的转动惯量张量以及由它们推导出的若干预计算标量：

```julia
struct QC{T,JT} <: AbstractCoordinates{3,T}
    m      # 质量
    m⁻¹    # 1/质量
    γ      # 惯量对角元素最大值（正则化参数）
    Jγ     # J - γ·I
    γ⁻¹    # 1/γ
    J⁻¹γ   # J⁻¹ - (1/γ)·I
end
```

通过传入质量和 3×3 惯量矩阵来构造 `QC` 实例：

```julia
m = 0.5
J = Diagonal(SA[0.001, 0.001, 0.002])

qc = QCF.QC(m, J)
```

惯量矩阵仅有对角元素被使用，非对角项会在内部被丢弃。可选关键字参数 `γ` 默认为 `maximum(diag(J))`，作为正则化参数保持旋转质量矩阵的良好条件数。

!!! tip "查询坐标特性"
    可以在任意层级查询刚体是否使用恒定质量坐标：
    ```julia
    has_constant_mass_matrix(body.coords)   # QC 返回 Val(false)
    has_constant_mass_matrix(body)          # 从 coords 向上传播
    has_constant_mass_matrix(structure)     # 从 bodies 向上传播
    ```

## 非恒定质量矩阵

与 NCF（`has_constant_mass_matrix` 返回 `Val(true)`）不同，QCF 的广义质量矩阵**依赖于当前四元数状态**：

```julia
has_constant_mass_matrix(::QCF.QC) = Val(false)
```

7×7 质量矩阵为块对角形式：

```math
M(q) = \begin{bmatrix} mI_3 & 0 \\ 0 & M_\gamma(q) \end{bmatrix}, \quad M_\gamma(q) = 4\bigl(L^\top(q)\,J_\gamma\,L(q) + \gamma I_4\bigr)
```

其中 ``L(q)`` 是 4×4 四元数左乘矩阵。由于 ``M`` 每步都变化，求解器必须在每次迭代中重新组装和分解它。Rible 通过求解器分派自动处理这一点：当检测到 `Val(false)` 时，会选择 `Zhong06_Nonconstant_Mass` 求解器族而非恒定质量族。

### 实际影响

- **每步计算代价**高于 NCF，因为需要逐步组装并分解质量矩阵。对于单体仿真差异可忽略；对于包含大量 QCF 刚体的多体系统，需权衡 NCF 是否足够。
- **紧容差有助于稳定性**。四元数约束 ``\|q\|^2 = 1`` 在求解器容差较松时会漂移。推荐默认使用 `ftol = 1e-14`、`maxiters = 50`。

## 内禀约束

`QC` 携带一个内禀代数约束——四元数单位模长：

```math
c(q) = \frac{1}{2}(q^\top q - 1) = 0
```

此约束由求解器与关节约束一起自动执行。框架在 `constraints.jl` 中提供了约束函数、其雅可比矩阵以及力雅可比矩阵。

## 完整示例：基于 QCF 的旋转陀螺

内置的旋转陀螺模型通过 `cT` 参数同时支持 NCF 和 QCF。下面我们构造一个 QCF 陀螺，验证其非恒定质量特性，并执行完整的接触仿真。

### 构造

```@example qcf
include(joinpath(pathof(Rible), "../../examples/robots/spinningtop.jl"))

origin_position = [0.0, 0.0, 0.5]
R = RotX(0.0)
origin_velocity = [1.0, 0.0, 0.0]
Ω = [0.0, 0.0, 200.0]

top = make_top(origin_position, R, origin_velocity, Ω, :QCF; μ = 0.95, e = 0.5, loadmesh = true)
```

### 验证非恒定质量特性

```@example qcf
body1 = get_bodies(top.structure)[1]
has_constant_mass_matrix(body1.coords) === Val(false)
```

### 接触设置与仿真

求解器对象与 NCF 相同——一个带有内层接触求解器的 `Zhong06` 积分器。框架在内部将 QCF 刚体路由到非恒定质量代码路径。

```@example qcf
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

sim_result = solve!(
    prob,
    solver;
    tspan = (0.0, 0.2),
    dt = 1e-3,
    ftol = 1e-14,
    maxiters = 50,
    exception = false,
)

tip_pos = get_trajectory!(top, 1, 1).u[end]
tip_vel = get_velocity!(top, 1, 1).u[end]
tip_pos, tip_vel
```

### 可视化

由于陀螺在构造时使用了 `loadmesh = true`，`vis!` 会直接渲染 STL 几何体：

```@example qcf
key_steps = round.(Int, range(1, length(top.traj.t), length = 4))

fig_traj = plot_traj!(
    top;
    do_slide = false,
    show_info = false,
    show_loci = false,
    show_background = false,
    show_ground = false,
    gridsize = (2, 2),
    at_steps = key_steps,
)
```

### 对比 QCF 与 NCF 结果

可以用相同的初始条件和 `:NCF` 重建模型进行对比：

```@example qcf
top_ncf = make_top(origin_position, R, origin_velocity, Ω, :NCF; μ = 0.95, e = 0.5, loadmesh = false)

has_constant_mass_matrix(get_bodies(top_ncf.structure)[1].coords) === Val(true)
```

## 参考文献

```@bibliography
Pages = ["RibleQCF.md"]
```
