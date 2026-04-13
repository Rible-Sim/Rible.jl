
rotation_matrix(θ::Real) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]
rotation_matrix(R::AbstractMatrix) = R
function angle2mode(θ::Real)
    modθ = mod2pi(θ+π/4) 
    if modθ < π/2
        return 1
    elseif modθ < π
        return 2
    elseif modθ < 3π/2
        return 3
    else
        return 4
    end
end

"""
    cartesian_to_angular_velocity(x, y, vx, vy)

Computes angular velocity ω = dθ/dt given Cartesian position and velocity.
Assumes θ = atan(y, x).
Formula: ω = (x*vy - y*vx) / (x^2 + y^2)
"""
function cartesian_to_angular_velocity(x, y, vx, vy)
    r2 = x^2 + y^2
    return (x*vy - y*vx) / r2
end