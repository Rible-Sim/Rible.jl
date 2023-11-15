# Quaternion

quaternion as a 4-components vector

$$\bm{q} = [q_0, q_1, q_2, q_3]^\mathrm{T} \in \mathbb{R}^4$$

unit quaternion encode axis-angle ($\bm{n}$, $\theta$) rotation as 

$$\bm{q} = [
    \cos\frac{\theta}{2}, 
    n_x\sin\frac{\theta}{2}, 
    n_y\sin\frac{\theta}{2}, 
    n_z\sin\frac{\theta}{2}
]$$

so it satisfies

$$ \varPhi(\bm{q}) = \frac{1}{2}(\bm{q}^\mathrm{T}\bm{q} - 1) = 0$$

Jacobian of the above constraint is 

$$ \frac{\partial \varPhi}{\partial \bm{q}} = \bm{q}^\mathrm{T}\dot{\bm{q}} = 0$$

which becomes a velocity-level constraint.

If the above position-level constraint is satisfied, conjugate is inverse  

$$ \bm{q}* = [q_0, -q_1, -q_2, -q_3]^\mathrm{T}$$

pure quaternion $\bm{v}$ is 

$$\bm{v} = [0, v_1, v_2, v_3]^\mathrm{T}$$

rotation matrix 

$$\bm{R}(\bm{q}) = \begin{bmatrix}
    
\end{bmatrix}$$

Let $\bar{\bm{\eta}}$ be either a 3-element vector or a pure quaternion.

rotate $\bar{\bm{\eta}}$ using either

$$\bm{\eta} = \bm{R}(\bm{q})\bar{\bm{\eta}}$$

or equivalently 

$$ \bm{\eta} = \bm{q}\bar{\bm{\eta}}\bm{q}*$$

Local angular velocity

$$\Omega = 2L(\bm{q})\dot{\bm{q}}$$

Velocity of the rotated vector

$$\dot{\bm{\eta}} = \bm{R}(\dot{\bm{q}})\bar{\bm{\eta}}$$

Jacobian of a rotated vector 

$$ \frac{\partial
    \left(
        \bm{R}(\bm{q})\bm{\eta}
    \right)
    }{\partial \bm{q}} 
=
$$

