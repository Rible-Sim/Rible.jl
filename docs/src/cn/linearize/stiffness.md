# 刚度

stiffness.jl 文件实现了一系列用于分析和优化张拉整体结构刚度的函数。主要功能包括:

1. 静力学和运动学分析
2. 刚度矩阵的构建和优化
3. 自应力状态和刚度方向的计算

## 主要函数

### static_kinematic_determine

此函数执行张拉整体结构的静力学和运动学分析。

参数:

- `ℬᵀ`: 平衡矩阵的转置
- `atol`: 绝对容差(可选)

返回:

- 自应力状态
- 刚度方向

### optimize_maximum_stiffness 和 optimize_zero_stiffness

这两个函数使用凸优化来最大化或最小化结构刚度。

参数:

- `mat𝒦ps`: 预应力刚度矩阵
- `vec𝒦m`: 材料刚度矩阵
- `vecI`: 单位矩阵向量
- `A`, `b`: 等式约束
- `nx`: 变量数量

返回:

- 优化结果

## 常见用法

1. 分析结构的自应力状态:

```julia
S, D = static_kinematic_determine(ℬᵀ)
```

2. 优化最大刚度:

```julia
result = optimize_maximum_stiffness(mat𝒦ps, vec𝒦m, vecI, A, b, nx)
```

3. 寻找零刚度构型:

```julia
result = optimize_zero_stiffness(mat𝒦ps, vec𝒦m, vecI, A, b, nx, x_0)
```

## 内部辅助函数

- `classical_gram_schmidt` 和 `modified_gram_schmidt`: 用于正交化矩阵列
- `optimize_zero_stiffness_Clarabel`: 使用 Clarabel 求解器的零刚度优化版本

这些功能共同提供了一个强大的工具集,用于分析和优化张拉整体结构的刚度特性。
