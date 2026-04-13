# Built-in Policies

Rible.jl comes with several ready-to-use policies for common control and trajectory optimization tasks.

## Open-Loop Policies

### `NoPolicy`
Represents the absence of a control policy. When this policy is active, no control actions are applied to the robot (actions are zeroed or remain unchanged depending on initialization, but effectively "do nothing" logic).

```julia
struct NoPolicy <: AbstractPolicy end
```

### `TimeFunctionPolicy`
A policy where control actions are a function of time only: $u = \pi(t)$.

```julia
struct TimeFunctionPolicy{F <: Function} <: AbstractTimePolicy
    f::F # Typically holds a function f(t)
end
```

### `TrajectoryPolicy`
Represents a pre-computed open-loop trajectory of control actions. It interpolates actions between specified time nodes.

```julia
struct TrajectoryPolicy{T,time_nodesT,actuations_trajT,itpT} <: AbstractPolicy
    num_of_actions::T
    time_nodes::time_nodesT
    actuations_traj::actuations_trajT
    itp::itpT # Interpolation objects
end
```
- **Constructor**: `TrajectoryPolicy(time_nodes, actuations_traj)`
- **Parameters**: The values in `actuations_traj` are treated as optimization parameters.

## Closed-Loop / Feedback Policies

### `StaticLinearPolicy`
Implements a simple linear feedback controller: $u = Kx + b$, where $x = [q; \dot{q}]$.

```julia
struct StaticLinearPolicy{mat_biasT,matT,biasT} <: AbstractPolicy
    mat_bias::mat_biasT
    mat::matT
    bias::biasT
end
```
- **Parameters**: Entries of the bias vector and feedback matrix.

### `LuxPolicy`
Wraps a neural network actor (using Lux.jl) for learning-based control.

```julia
struct LuxPolicy{A,P,S} <: AbstractPolicy
    actor::A
    ps::P
    st::S
end
```
- **Usage**: Used in Reinforcement Learning contexts (e.g., `RibleQCF`).
- **Parameters**: Neural network weights `ps`.

### `iLQRRolloutPolicy`
A policy specialized for Iterative Linear Quadratic Regulator (iLQR) rollouts. It combines a feedforward term, feedback gains, and a line-search step size.

```julia
struct iLQRRolloutPolicy{T, VK, Vd, Vx, Vu, Vss} <: AbstractPolicy
    time_nodes::T
    K::VK
    d::Vd
    x::Vx
    u::Vu
    step_size::Vss
end
```
- **Logic**: $u = u_{ref} + K(x - x_{ref}) + \alpha \cdot d$
- **Parameters**: `step_size` ($\alpha$).
