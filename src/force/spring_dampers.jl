abstract type AbstractForce end
struct NoForce <: AbstractForce end

get_num_of_aux_var(::AbstractForce) = 0


"""
$(TYPEDEF)
"""
mutable struct DistanceSpringDamperState{N,T}
    restlen::T
    length::T
    lengthdot::T
    tension::T
    direction::MArray{Tuple{N},T,1,N}
    force::MArray{Tuple{N},T,1,N}
    start::MArray{Tuple{N},T,1,N}
    stop::MArray{Tuple{N},T,1,N}
    start_vel::MArray{Tuple{N},T,1,N}
    stop_vel::MArray{Tuple{N},T,1,N}
end

"""
$(TYPEDSIGNATURES)
"""
function DistanceSpringDamperState(restlen,direction)
    DistanceSpringDamperState(
        restlen,
        restlen,
        zero(restlen),
        zero(restlen),
        direction,
        zero(direction),
        zero(direction),
        zero(direction),
        zero(direction),
        zero(direction)
        )
end

"""
$(TYPEDEF)
"""
struct DistanceSpringDamper{N,T} <: AbstractForce
    k::T
    c::T
	slack::Bool
    state::DistanceSpringDamperState{N,T}
end

"""
$(TYPEDSIGNATURES)
"""
function DistanceSpringDamper2D(restlen::T,k::T,c=zero(k);slack=true) where T
    direction = MVector{2}(one(T),zero(T))
    state = DistanceSpringDamperState(restlen,direction)
    DistanceSpringDamper(k,c,slack,state)
end

"""
$(TYPEDSIGNATURES)
"""
function DistanceSpringDamper3D(restlen::T,k::T,c=zero(k);slack=true) where T
    direction = MVector{3}(one(T),zero(T),zero(T))
    state = DistanceSpringDamperState(restlen,direction)
    DistanceSpringDamper(k,c,slack,state)
end

"""
$(TYPEDSIGNATURES)
"""
function update!(cab::DistanceSpringDamper,p1,p2,ṗ1,ṗ2)
	(;k,c,state) = cab
    state.start .= p1
    state.stop .= p2
    state.start_vel .= ṗ1
    state.stop_vel .= ṗ2
	@. state.direction = p2 - p1
	len2 = dot(state.direction, state.direction)
	len = sqrt(len2)
	state.length = len
	inv_len = inv(len)
	@. state.direction *= inv_len
	@. state.force = ṗ2 - ṗ1
	state.lengthdot = dot(state.direction, state.force)
	deformation = len - state.restlen
	f = k*deformation + c*state.lengthdot
	state.tension = cab.slack ? ((deformation < 0 || f < 0) ? zero(f) : f) : f
	@. state.force = state.tension * state.direction
end

"""
$(TYPEDSIGNATURES)
"""
function potential_energy(cab::DistanceSpringDamper)
	(;k,state,slack) = cab
	Δ⁺ = state.length - state.restlen
	Δ = ifelse(slack,max(0,Δ⁺),Δ⁺)
	pe = 1/2*k*Δ^2
end


"""
$(TYPEDEF)
"""
mutable struct TorsionalSpringDamperState{T}
    rest_angle::T
    angle::T
    angular_velocity::T
    torque::T
    mode_id::Int
end

"""
$(TYPEDSIGNATURES)
"""
function TorsionalSpringDamperState(rest_angle)
    TorsionalSpringDamperState(
        rest_angle,
        zero(rest_angle),
        zero(rest_angle),
        zero(rest_angle),
        angle2mode(rest_angle)
    )
end

"""
$(TYPEDEF)
"""
struct TorsionalSpringDamper{T} <: AbstractForce
    k::T
    c::T
	slack::Bool
    state::TorsionalSpringDamperState{T}
end

"""
$(TYPEDSIGNATURES)

Create a torsional spring damper with specified rest angle, stiffness k, and optional damping coefficient c.

# Arguments
- `rest_angle::T`: The rest angle of the spring in radians
- `k::T`: Spring stiffness coefficient 
- `c=zero(k)`: Optional damping coefficient, defaults to 0
- `slack=false`: If true, spring only generates torque when angle > rest_angle

# Returns
- `TorsionalSpringDamper`: A torsional spring damper force element
"""
function TorsionalSpringDamper(rest_angle::T,k::T,c=zero(k);slack=false) where T
    state = TorsionalSpringDamperState(rest_angle)
    TorsionalSpringDamper(k,c,slack,state)
end

get_num_of_aux_var(::TorsionalSpringDamper) = 1

function get_auxilary(force::TorsionalSpringDamper)
    [force.state.angle]
end

"""
$(TYPEDSIGNATURES)
"""
function update!(sd::TorsionalSpringDamper,angle,angular_velocity=zero(angle))
	(;k,c,state) = sd
	deformation = angle - rest_angle
	state.angle = angle
	state.angular_velocity = angular_velocity
	state.torque = k*deformation + c*state.angular_velocity
end

"""
$(TYPEDSIGNATURES)
"""
function potential_energy(sd::TorsionalSpringDamper)
	(;k,c,state,slack) = sd
	Δ⁺ = (state.angle - state.rest_angle)
	Δ = ifelse(slack,max(0,Δ⁺),Δ⁺)
    ## Δ = Δ⁺
	pe = 1/2*k*(Δ^2)
end
