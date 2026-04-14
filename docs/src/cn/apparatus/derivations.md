# Generalized Force Derivation for Torsional Spring (NCF)

## Goal
To derive the generalized force vector $\mathbf{Q}(q, s)$ acting on the natural coordinates $q$, where the force *magnitude* (torque) is parameterized by the auxiliary variable $s$ (the angle $\theta$), while the *direction* is determined by the system kinematics $q$.

## 1. Virtual Work Formulation
The virtual work $\delta W$ done by the torsional spring torque $\tau$ for a virtual rotation $\delta \theta$ is:
$$ \delta W = \tau \delta \theta $$

We define the torque $\tau$ as a function of the auxiliary variable $s$:
$$ \tau(s) = -k(s - s_0) $$
where $s = \theta$ is the state variable representing the angle.

The system angle $\theta$ is also related kinematically to the natural coordinates $q$.
Let the body orientation be defined by the unit vector $\mathbf{u}(q)$ (e.g., in NCF, $\mathbf{u}$ is a basis vector derived from $q$).
The kinematic variation $\delta \theta$ is related to the variation in coordinates $\delta q$ by the Jacobian $\mathbf{J}_\theta(q) = \nabla_q \theta$:
$$ \delta \theta = \mathbf{J}_\theta(q) \cdot \delta q $$

Substituting this into the virtual work equation:
$$ \delta W = \tau(s) (\mathbf{J}_\theta(q) \cdot \delta q) $$

The generalized force vector $\mathbf{Q}(q, s)$ is the coefficient of $\delta q$ in the virtual work expression:
$$ \mathbf{Q}(q, s) = \tau(s) \mathbf{J}_\theta(q)^T $$

This explicitly shows the structure requested:
*   **Magnitude** $\tau(s)$ depends on the auxiliary variable $s$.
*   **Direction** $\mathbf{J}_\theta(q)^T$ depends on the generalized coordinates $q$.

## 2. Sensitivity Jacobian w.r.t. Auxiliary Variable ($\nabla_s \mathbf{Q}$)
To implement the Adjoint method or implicit integration, we need the Jacobian of the generalized force with respect to the auxiliary variable $s$.

Taking the partial derivative of $\mathbf{Q}(q, s)$ with respect to $s$:
$$ \nabla_s \mathbf{Q} = \frac{\partial}{\partial s} \left( \tau(s) \mathbf{J}_\theta(q)^T \right) $$

Since $\mathbf{J}_\theta(q)$ does not depend on $s$:
$$ \nabla_s \mathbf{Q} = \frac{\partial \tau}{\partial s} \mathbf{J}_\theta(q)^T $$

Given $\tau(s) = -k(s - s_0)$, the derivative is:
$$ \frac{\partial \tau}{\partial s} = -k $$

Thus, the sensitivity Jacobian is:
$$ \nabla_s \mathbf{Q} = -k \mathbf{J}_\theta(q)^T $$

## 3. Explicit Form for NCF (Natural Coordinates)
In Natural Coordinates (NCF), the orientation is often defined by a vector $\mathbf{u} = [u_x, u_y]^T$ derived linearly from $q$.
Assuming $\mathbf{u}$ is a unit vector, the angle $\theta$ satisfies:
$$ \mathbf{u} = [\cos\theta, \sin\theta]^T $$

The variation $\delta \theta$ is related to $\delta \mathbf{u}$ by the cross-product relationship:
$$ \delta \theta = \mathbf{v} \cdot \delta \mathbf{u} $$
where $\mathbf{v} = [-u_y, u_x]^T$ is the vector orthogonal to $\mathbf{u}$.

Since $\mathbf{u}$ is linear in $q$ (for standard NCF elements), $\delta \mathbf{u} = \frac{\partial \mathbf{u}}{\partial q} \delta q$.
Let $\mathbf{C} = \frac{\partial \mathbf{u}}{\partial q}$ be the constant conversion matrix (or block thereof). Then:
$$ \delta \theta = \mathbf{v}^T \mathbf{C} \delta q $$
So the kinematic Jacobian is:
$$ \mathbf{J}_\theta(q) = \mathbf{v}^T \mathbf{C} $$

Substituting back into the force equations:

**Generalized Force:**
$$ \mathbf{Q}(q, s) = \tau(s) \mathbf{C}^T \mathbf{v} $$

**Sensitivity Jacobian:**
$$ \nabla_s \mathbf{Q} = -k \mathbf{C}^T \mathbf{v} $$

Where $\mathbf{v} = [-u_y, u_x]^T$ is computed from the current configuration $q$.

## 4. Implementation Correspondence
*   **Primal Force (`mutate.jl`)**: Accumulates $\mathbf{Q} = \tau(s) \mathbf{C}^T \mathbf{v}$.
*   **Adjoint Jacobian (`linearization.jl`)**: Accumulates $\nabla_s \mathbf{Q} = -k \mathbf{C}^T \mathbf{v}$.

This formulation is robust (no division by Jacobian determinant) and mathematically consistent with the Separation of Variables approach $\mathbf{Q}(q, s) = f(s) g(q)$.
