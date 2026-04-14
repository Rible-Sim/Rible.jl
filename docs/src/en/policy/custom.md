# Custom Policies

To implement a custom control logic in Rible.jl, you need to define a new struct inheriting from `AbstractPolicy` and implement the `actuate!` method.

## Basic Implementation

1. **Define the Struct**
   ```julia
   struct MyCustomPolicy <: AbstractPolicy
       # Add fields as needed (e.g., gains, target positions)
   end
   ```

2. **Implement `actuate!`**
   This function writes the computed control actions into the robot's state.
   
   ```julia
   import Rible: actuate!
   
   function actuate!(bot::Robot, policy::MyCustomPolicy, inst_state::AbstractCoordinatesState)
       (;structure,hub) = bot
       (; q, q̇, t) = inst_state
       (; actuators, state, coalition) = hub
       (; actid2sys_actions_idx) = coalition
       
       # 1. Compute your controls
       my_u = ... # calculate based on q, q̇, t
       
       # 2. Assign to global control vector
       # Note: This example assumes a simple mapping. 
       # In practice, iterate over actuators to map correctly.
       state.u .= ... 
       
       # 3. Execute on actuators
       foreach(actuators) do actuator
           idx = actid2sys_actions_idx[actuator.id]
           execute!(
               structure,
               actuator,
               (@view state.u[idx])
           )
       end
   end
   ```

## Optimization Interface

If your policy has tunable parameters (e.g., for trajectory optimization or learning), implement the parameter interface:

```julia
import Rible: get_params, get_num_of_params, set_params!

# Return a flat vector of parameters
function get_params(policy::MyCustomPolicy)
    return ...
end

# Return the total number of parameters
function get_num_of_params(policy::MyCustomPolicy)
    return ...
end

# Update internal parameters from a flat vector
function set_params!(policy::MyCustomPolicy, params)
    # update fields...
end
```

## Differentiable Control

For gradient-based methods (e.g., trajectory optimization, reinforcement learning), Rible.jl requires derivatives of the policy actions.

### 1. Continuous-Time / Force Jacobian

When simulating dynamics, the system need to know how generalized forces change with state.

```julia
import Rible: gen_force_state_jacobian!

function gen_force_state_jacobian!(∂F∂q̌, ∂F∂q̌̇, bot::Robot, policy::MyCustomPolicy, inst_state::AbstractCoordinatesState)
    (;structure, hub) = bot
    # Populate ∂F∂q̌ and ∂F∂q̌̇ (matrices)
    # This usually involves:
    # 1. Computing ∂u/∂x (Jacobian of control action w.r.t state)
    # 2. Computing ∂F/∂u (Jacobian of generalized force w.r.t control action)
    # 3. Applying chain rule: ∂F/∂x = ∂F/∂u * ∂u/∂x
end
```

### 2. Discrete-Time Optimization Interface

For solvers that operate on discrete time steps (like iLQR), you must provide the Jacobian of the *action* $u_k$ with respect to the discrete state variables ($q_k, q_{k+1}, p_k, p_{k+1}, \dots$) and the policy parameters $\theta$.

> [!NOTE]
> These functions dispatch on the `solver` type. While the interface is generic, the specific arguments (like `qₖ`, `qₖ₊₁`) often depend on the solver's discretization scheme. Currently, this interface is primarily designed for the **Zhong06** family of symplectic integrators.

#### Reverse Mode (VJP)
For efficiency in reverse-mode differentiation (e.g., adjoint sensitivity), implementing `vjp_wrt_state` and `accumulate_param_grad!` is preferred.

```julia
import Rible: vjp_wrt_state, accumulate_param_grad!

function vjp_wrt_state(v, policy::MyCustomPolicy, bot::Robot, num_of_actions, solver::Zhong06, solver_state)
    # v is the cotangent vector (adjoint) of size (num_of_actions,)
    
    ∂ϕ∂qₖᵀ   = ... # v' * ∂u/∂qₖ
    ∂ϕ∂qₖ₊₁ᵀ = ...
    # ... etc
    
    return ∂ϕ∂qₖᵀ, ∂ϕ∂qₖ₊₁ᵀ, ∂ϕ∂pₖᵀ, ∂ϕ∂pₖ₊₁ᵀ, ∂ϕ∂λᵀ, ∂ϕ∂sₖᵀ, ∂ϕ∂sₖ₊₁ᵀ
end

function accumulate_param_grad!(grad_storage, policy::MyCustomPolicy, v_total, solver_state, bot)
    # Add parameter gradients: grad_storage .+= transpose(∂u/∂θ) * v_total
end
```
