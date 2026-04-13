# 额外积分器 RibleExtraIntegrators

```@meta
ShareDefaultModule = true
```

`RibleExtraIntegrators` 为核心框架补充了额外的时间积分算法，适用于默认 Zhong06 家族之外的特殊场景。当前最直接的用途，是为非光滑动力学提供 `Moreau` 这类替代积分器。

## 广义α方法 GeneralizedAlpha

`GeneralizedAlpha` 是一种面向结构动力学的高阶隐式积分算法。它可在高频数值阻尼与低频响应精度之间进行可控折中，常用于刚柔耦合系统中的长时间积分。

它的主要特点包括：

- 面向目标问题保持二阶精度，
- 通过高频谱半径 ``\rho_\infty`` 控制数值耗散，
- 对刚性较强的结构系统具有良好的长时间稳定性。

### 参数配置

数值耗散由 ``\rho_\infty`` 控制。当 ``\rho_\infty = 1`` 时，算法接近无阻尼极限；当 ``\rho_\infty = 0`` 时，则会施加最强的高频阻尼。

```julia
using Rible
using RibleExtraIntegrators

solver = DynamicsSolver(
    RibleExtraIntegrators.GeneralizedAlpha(0.8)
)

solve!(prob, solver; tspan=(0.0, 10.0), dt=1e-3)
```

## Moreau-Jean 方法 Moreau

`Moreau` 实现了经典的非光滑时间步进法，适用于碰撞、单边接触和摩擦等会引入速度跳变的问题。对于这类系统，它通常比逐事件的平滑积分更直接，也更稳健。

### 理论基础

在 Moreau-Jean 框架下，动力学可写成如下测度形式：

```math
M \, \mathrm{d}v = f(t, q, v) \, \mathrm{d}t + \mathrm{d}r
```

其中 ``M`` 是质量矩阵，``v`` 是广义速度，``f`` 是平滑外力/内力项，``\mathrm{d}r`` 是由接触引起的非光滑冲量项。

在给定 ``\theta`` 之后，连续问题会被离散为一系列互补型子问题。

### 构造方式

`Moreau` 积分器接受离散化参数 ``\theta``：

```julia
using Rible
using RibleExtraIntegrators

solver = DynamicsSolver(
    RibleExtraIntegrators.Moreau(0.5)
)

solve!(prob, solver; tspan=(0.0, 0.1), dt=1e-3)
```

### 测试用例

包的测试套件中包含一个滑块-曲柄基准测试，用于比较 `Moreau(0.5)` 与默认 `Zhong06()` 在同一模型上的表现。该测试验证了两种方法在轨迹和能量方面的吻合程度。

完整基准测试代码见 [`RibleExtraIntegrators/test/slider_crank.jl`](https://github.com/Rible-Sim/Rible.jl/tree/main/RibleExtraIntegrators/test/slider_crank.jl)。
