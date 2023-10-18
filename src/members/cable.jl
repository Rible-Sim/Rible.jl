"""
绳索信息类。
$(TYPEDEF)
"""
mutable struct CableState{N,T}
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
绳索信息构造子。
$(TYPEDSIGNATURES)
"""
function CableState(restlen,direction)
    CableState(
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
普通绳索信息类。
$(TYPEDEF)
"""
struct Cable{N,T}
    id::Int
    k::T
    c::T
	slack::Bool
    state::CableState{N,T}
end

"""
2D普通绳索构造子。
$(TYPEDSIGNATURES)
"""
function Cable2D(id,restlen::T,k::T,c=zero(k);slack=true) where T
    direction = MVector{2}(one(T),zero(T))
    state = CableState(restlen,direction)
    Cable(id,k,c,slack,state)
end

"""
3D普通绳索构造子。
$(TYPEDSIGNATURES)
"""
function Cable3D(id,restlen::T,k::T,c=zero(k);slack=true) where T
    direction = MVector{3}(one(T),zero(T),zero(T))
    state = CableState(restlen,direction)
    Cable(id,k,c,slack,state)
end

"""
更新绳索拉力。
$(TYPEDSIGNATURES)
"""
function update!(cab::Cable,p1,p2,ṗ1,ṗ2)
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
计算绳索势能。
$(TYPEDSIGNATURES)
"""
function potential_energy(cab::Cable)
	(;k,state,slack) = cab
	Δ⁺ = state.length - state.restlen
	Δ = ifelse(slack,max(0,Δ⁺),Δ⁺)
	pe = 1/2*k*Δ^2
end

mutable struct LinearLaw{T}
    k::T
    F0::T
end

function (ll::LinearLaw)(Δl)
    @unpack F0, k = ll
    F = F0 + k*Δl
end

"""
？记忆合金绳索（一种拉力与温度有关的绳索）信息类。
$(TYPEDSIGNATURES)
"""
mutable struct SMACableState{N,T}
    temp::T
    restlen::T
    length::T
    lengthdot::T
    tension::T
    direction::MArray{Tuple{N},T,1,N}
end

function SMACableState(temp,restlen,direction)
    SMACableState(temp,restlen,restlen,zero(restlen),zero(restlen),direction)
end


struct SMACable{N,T,F}
    id::Int
    law::F
    c::T
    state::SMACableState{N,T}
end

function SMACable3D(id,restlen::T,law::LawT,c::T,original_temp::T=0.0) where {LawT,T}
    direction = MVector{3}(one(T),zero(T),zero(T))
    state = SMACableState(original_temp,restlen,direction)
    SMACable(id,law,c,state)
end

mutable struct SlidingPoint{T}
    μ::T
    θ::T
    α::T
    s::T
    s⁺::T
    s⁻::T
end

struct CableSegment{N,T}
    id::Int
    k::T
    c::T
    prestress::T
    original_restlen::T
    state::CableState{N,T}
end

function CableSegment(id, original_restlen, k; c=zero(k), prestress=zero(k))
    direction = MVector{2}(one(k), zero(k))
    state = CableState(original_restlen, direction)
    CableSegment(id, k, c, prestress, original_restlen, state)
end

function calculate_α(μ,θ)
    exp(-μ*θ)
end

function SlidingPoint(μ)
    θ = one(μ) * 2pi
    #θ = zero(μ)
    α = calculate_α(μ,θ)
    s = zero(μ)
    s⁺,s⁻ = s2s̄(s)
    SlidingPoint(μ,θ,α,s,s⁺,s⁻)
end

struct ClusterCables{spsType,segsType}
    ID::Int
    sps::spsType
    segs::segsType
end

function ClusterCables(id, nsp, segs;μ=0.0)
    sps = StructArray([SlidingPoint(μ) for i = 1:nsp])
    ClusterCables(id, sps, segs)
end

function s2s̄(s::Number)
    abss = abs(s)
    s⁺ = (abss + s)/2
    s⁻ = (abss - s)/2
    s⁺,s⁻
end