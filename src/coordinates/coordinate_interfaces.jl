

# Coordinate System Interface

# This module defines the interface for coordinate systems in Rible.
# All coordinate implementations should implement these methods.



"""
Return the indices of intrinsic constraints for the given coordinates.
$(TYPEDSIGNATURES)
"""
function get_cstr_idx(coords)
    collect(1:get_num_of_intrinsic_cstr(coords))
end


"""
Return the number of local dimensions.
$(TYPEDSIGNATURES)
"""
function get_num_of_local_dims(coords) end


# Transformation methods

"""
Transform global coordinates to local coordinates.
$(TYPEDSIGNATURES)
"""
function to_local_coords(coords, r̄) end

"""
Get position from coordinates and local parameters.
$(TYPEDSIGNATURES)
"""
function to_position(coords, q, c) end

"""
Get position Jacobian with respect to coordinates.
$(TYPEDSIGNATURES)
"""
function to_position_jacobian(coords, q, c) end

"""
Get velocity Jacobian with respect to coordinates and velocities.
$(TYPEDSIGNATURES)
"""
function to_velocity_jacobian(coords, q, q̇, c) end

# Rotation and angular velocity methods

"""
Find rotation matrix from coordinates.
$(TYPEDSIGNATURES)
"""
function find_rotation(coords, q) end

"""
Find angular velocity from coordinates and velocities.
$(TYPEDSIGNATURES)
"""
function find_angular_velocity(coords, q, q̇) end

"""
    find_local_angular_velocity(coords, q, q̇)

Find local angular velocity from coordinates and velocities.
"""
function find_local_angular_velocity(coords, q, q̇) end

# Conversion methods

"""
Convert Cartesian frame to coordinates.
$(TYPEDSIGNATURES)
"""
function cartesian_frame2coords(coords, frame) end

"""
Get deformation parameters.
$(TYPEDSIGNATURES)
"""
function get_deform(coords) end

# Joint methods

"""
    build_joint_cache(coords, ...)

Build joint cache for coordinate system.
"""
function build_joint_cache end

"""
    get_joint_violations!(violations, coords, ...)

Get joint constraint violations.
"""
function get_joint_violations! end

"""
    get_joint_jacobian!(jacobian, coords, ...)

Get joint constraint jacobian.
"""
function get_joint_jacobian! end

"""
    add_joint_forces_jacobian!(jacobian, coords, ...)

Get joint forces jacobian.
"""
function add_joint_forces_jacobian! end

"""
    get_joint_velocity_jacobian!(jacobian, coords, ...)

Get joint velocity jacobian.
"""
function get_joint_velocity_jacobian! end

# Utility methods

"""
    find_independent_free_idx(coords, ...)

Find independent free indices.
"""
function find_independent_free_idx(coords, q) end

"""
    nullspace_mat(coords, ...)

Get nullspace matrix.
"""
function nullspace_mat(coords, q) end

"""
    kinetic_energy_coords(coords, ...)

Calculate kinetic energy from coordinates.
"""
function kinetic_energy_coords end
