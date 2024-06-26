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
struct DistanceSpringDamper{N,T}
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
	# state.direction .= p2 .- p1 # Δr
	# state.length = norm(state.direction)
	# state.force .= ṗ2 - ṗ1 # Δṙ
	# state.lengthdot = (transpose(state.direction)*state.force)/state.length
	# state.direction ./= state.length
    state.start,     state.stop = p1, p2
    state.start_vel, state.stop_vel = ṗ1, ṗ2
	l = p2 - p1
	l̇ = ṗ2 - ṗ1
	state.length = norm(l)
	state.direction = l/state.length
	state.lengthdot = (transpose(state.direction)*l̇)
	deformation = state.length - state.restlen
	f = k*deformation + c*state.lengthdot
	if cab.slack
		if deformation < 0
			state.tension = 0.0
		elseif f < 0
			state.tension = 0.0
		else
			state.tension = f
		end
	else
		state.tension = f
	end
	state.force .= state.tension.*state.direction
    
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
mutable struct RotationalSpringDamperState{M,T}
    rest_angles::MArray{Tuple{M},T,1,M}
    angles::MArray{Tuple{M},T,1,M}
    angular_velocities::MArray{Tuple{M},T,1,M}
    torques::MArray{Tuple{M},T,1,M}
end

"""
$(TYPEDSIGNATURES)
"""
function RotationalSpringDamperState(rest_angles)
    RotationalSpringDamperState(
        rest_angles,
        zero(rest_angles),
        zero(rest_angles),
        zero(rest_angles),
    )
end

"""
$(TYPEDEF)
"""
struct RotationalSpringDamper{N,T}
    k::T
    c::T
	slack::Bool
    mask::Vector{Int}
    state::RotationalSpringDamperState{N,T}
end

"""
$(TYPEDSIGNATURES)
"""
function RotationalSpringDamper2D(rest_angles::AbstractVector{T},mask,k::T,c=zero(k);slack=true) where T
    state = RotationalSpringDamperState(MVector{1}(rest_angles))
    RotationalSpringDamper(k,c,slack,mask,state)
end

"""
$(TYPEDSIGNATURES)
"""
function RotationalSpringDamper3D(rest_angles::AbstractVector{T},mask,k::T,c=zero(k);slack=true) where T
    state = RotationalSpringDamperState(MVector{3}(rest_angles))
    RotationalSpringDamper(k,c,slack,mask,state)
end

"""
$(TYPEDSIGNATURES)
"""
function update!(rf::RotationalSpringDamper,angles,angular_velocities=zero(angles))
	(;k,c,state) = rf
	state.angles .= angles
	state.angular_velocities .= angular_velocities
	deformations = state.angles - state.rest_angles
	state.torques .= k*deformations + c*state.angular_velocities
end

"""
$(TYPEDSIGNATURES)
"""
function potential_energy(rf::RotationalSpringDamper)
	(;k,state,mask,slack) = rf
	Δ⁺ = (state.angles - state.rest_angles)[mask]
	# Δ = ifelse(slack,max(0,Δ⁺),Δ⁺)
    Δ = Δ⁺
	pe = 1/2*k*sum(Δ.^2)
end

mutable struct LinearLaw{T}
    k::T
    F0::T
end

function (ll::LinearLaw)(Δl)
    (;F0, k) = ll
    F = F0 + k*Δl
end

"""
$(TYPEDSIGNATURES)
"""
mutable struct SMADistanceSpringDamperState{N,T}
    temp::T
    restlen::T
    length::T
    lengthdot::T
    tension::T
    direction::MArray{Tuple{N},T,1,N}
end

function SMADistanceSpringDamperState(temp,restlen,direction)
    SMADistanceSpringDamperState(temp,restlen,restlen,zero(restlen),zero(restlen),direction)
end


struct SMADistanceSpringDamper{N,T,F}
    id::Int
    law::F
    c::T
    state::SMADistanceSpringDamperState{N,T}
end

function SMADistanceSpringDamper3D(restlen::T,law::LawT,c::T,original_temp::T=0.0) where {LawT,T}
    direction = MVector{3}(one(T),zero(T),zero(T))
    state = SMADistanceSpringDamperState(original_temp,restlen,direction)
    SMADistanceSpringDamper(law,c,state)
end

mutable struct SlidingPoint{T}
    μ::T
    θ::T
    α::T
    s::T
    s⁺::T
    s⁻::T
end

struct DistanceSpringDamperSegment{N,T}
    k::T
    c::T
    prestress::T
    original_restlen::T
    state::DistanceSpringDamperState{N,T}
end

function DistanceSpringDamperSegment( original_restlen, k; c=zero(k), prestress=zero(k))
    direction = MVector{3}(one(k), zero(k), zero(k))
    state = DistanceSpringDamperState(original_restlen, direction)
    DistanceSpringDamperSegment( k, c, prestress, original_restlen, state)
end

"""
$(TYPEDSIGNATURES)
"""
function update!(cab::DistanceSpringDamperSegment, p1, p2, ṗ1, ṗ2, s1, s2)
	(;k, c, state, prestress) = cab
    state.start,     state.stop = p1, p2
    state.start_vel, state.stop_vel = ṗ1, ṗ2
	l = p2 - p1
	l̇ = ṗ2 - ṗ1
	state.length = norm(l)
	state.direction = l/state.length
	state.lengthdot = (transpose(state.direction)*l̇)
	deformation = state.length - state.restlen - s2 + s1
	f = k*deformation + c*state.lengthdot + prestress
    if f < 0
        state.tension = 0.0
    else
        state.tension = f
    end
	state.force .= state.tension.*state.direction
end

function calculate_α(μ,θ)
    exp(-μ*θ)
end

function SlidingPoint(μ)
    θ = one(μ) * 2pi
    α = calculate_α(μ, θ)
    s = zero(μ)
    s⁺, s⁻ = s2s̄(s)
    SlidingPoint(μ, θ, α, s, s⁺, s⁻)
end

struct ClusterDistanceSpringDampers{spsType,segsType}
    sps::spsType
    segs::segsType
end

function ClusterDistanceSpringDampers(segs; μ=0.0)
    @show 11
    nsp = length(segs) - 1
    sps = StructArray([SlidingPoint(μ) for i = 1:nsp])
    ClusterDistanceSpringDampers(sps, segs)
end

function s2s̄(s::Number)
    abss = abs(s)
    s⁺ = (abss + s)/2
    s⁻ = (abss - s)/2
    s⁺,s⁻
end