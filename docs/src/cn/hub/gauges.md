# 量器器

定义了测量和误差评估的核心组件。这些组件用于监控和控制机器人的状态。

## 涉及类型

1. Captum 类型用于选择要测量的量, 例如 `PositionCaptum`, `VelocityCaptum`, `PosVelCaptum`, `AngularPositionCaptum`， `ActionCaptum` 等具体子类，分别用于如位置、速度、角度等。

3. Gauge 类型将Captum和Signifier组合在一起,形成一个完整的测量单元。
`Signifier` 用于标识被测量的特定部件或关节，可以是*体*、*器*，也可以是致动器。
`CaptumGauge` 为单纯测量量器；`ErrorGauge` 增加了参考值,用于计算误差。

## 涉及函数

- `measure()`: 执行实际的测量

- `measure_jacobian()`: 计算测量量关于Signifier状态的雅克比矩阵

- `measure_hessians()`: 计算测量量关于Signifier状态的黑塞矩阵

## 使用示例

```julia
# 创建一个位置Captum
pos_captum = PositionCaptum()

# 创建一个Signifier
sig = Signifier(body=some_body, pid=1)

# 创建一个CaptumGauge
gauge = CaptumGauge(1, sig, pos_captum)

# 使用Gauge进行测量
measurement = measure(structure, gauge)

# 创建一个ErrorGauge
error_gauge = ErrorGauge(1, sig, pos_captum, reference_position)

# 计算误差
error = measure(structure, error_gauge)
```

通过组合不同的Captum、Signifier和Gauge,可以灵活地定义各种测量和误差评估方案,为机器人的状态监控和控制提供基础。