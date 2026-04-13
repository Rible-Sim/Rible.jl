# Policy System

The Policy system in Rible.jl provides a flexible interface for defining how a robot's actuators are controlled over time or based on state. 

At its core, a policy determines the control actions $u$ applied to the system.

## Abstract Interface

All policies must subtype `AbstractPolicy`.

```julia
abstract type AbstractPolicy end
```

The main entry point for using a policy is the `actuate!` function, which applies the policy's logic to update the robot's control inputs.

### Core Functions

- `actuate!(bot::Robot, policy::AbstractPolicy, t::Real)`: Update the robot's actuation state at time `t`.
- `actuate!(bot::Robot, policy::AbstractPolicy, state::AbstractCoordinatesState)`: Update actuation based on the full instantaneous state (time, configuration, velocity).
- `get_params(policy)`: Retrieve optimization parameters (if any).
- `set_params!(policy, params)`: Update optimization parameters.
- `get_num_of_params(policy)`: Get the number of optimization parameters.

See [Custom Policies](custom.md) for details on implementing your own policy.
