

struct ClusterJoint{spType} <: AbstractJoint
    sliding_points::Vector{spType}
    num_of_cstr::Int
end

mutable struct SlidingPoint{T}
    μ::T
    θ::T
    α::T
    s::T
    s⁺::T
    s⁻::T
    ϵ⁺::T
    κ₁⁺::T
    κ₂⁺::T
    ϵ⁻::T
    κ₁⁻::T
    κ₂⁻::T
end

function SlidingPoint(μ::T; ϵ=1e-7one(T), κ₁=one(T), κ₂=one(T)) where T
    θ = one(μ) * 2pi
    α = exp(-μ*θ)
    s = zero(μ)
    s⁺ = (abs(s) + s)/2
    s⁻ = (abs(s) - s)/2
    SlidingPoint(μ, θ, α, s, s⁺, s⁻, ϵ, κ₁, κ₂, ϵ, κ₁, κ₂)
end

struct ClusterDistanceSpringDampers{seg} <: AbstractForce
    segments::Vector{seg}
end


struct DistanceSpringDamperSegment{N,T} <: AbstractForce
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

function prepare_cache!(appar::Apparatus{<:CableJoint,<:DistanceSpringDamperSegment}, cnt::AbstractConnectivity)
    T = get_numbertype(appar)
    nq = cnt.num_of_full_coords
    nd = get_num_of_dims(appar)
    (; cache, joint) = appar
    if (size(cache.K, 1) != nq) || (size(cache.K, 2) != nq)
        cache.K = spzeros(T, nq, nq)
    else
        zero!(cache.K)
    end
    if (size(cache.C, 1) != nq) || (size(cache.C, 2) != nq)
        cache.C = spzeros(T, nq, nq)
    else
        zero!(cache.C)
    end
    J = cache.named.J
    DJ = cache.named.DJ
    if (size(J, 1) != nd) || (size(J, 2) != nq)
        cache.named = @eponymtuple(
            J = zeros(T, nd, nq),
            DJ = zeros(T, nd, nq),
        )
    else
        zero!(J)
        zero!(DJ)
    end
    nothing
end

function ApparatusCache(joint::CableJoint, force::DistanceSpringDamperSegment, full_coords_idx, )
    T = get_numbertype(joint)
    nd = get_num_of_dims(joint)
    #note to be overwrite with prepare_cache!
    num_appar_full_idx = length(full_coords_idx)
    (; num_of_cstr) = joint

    φ = zeros(T, num_of_cstr)
    A = spzeros(T, num_of_cstr, num_appar_full_idx)
    ∂Aᵀλ∂q = zeros(T, num_appar_full_idx, num_appar_full_idx)
    ∂Aq̇∂q = zeros(T, num_of_cstr, num_appar_full_idx)
    Q = zeros(T, num_appar_full_idx)
    K = zeros(T, num_appar_full_idx, num_appar_full_idx)
    C = zeros(T, num_appar_full_idx, num_appar_full_idx)

    J = zeros(T, nd, num_appar_full_idx)
    DJ = similar(J)
    named = @eponymtuple(J, DJ)

    ApparatusCache(
        φ,
        A,
        ∂Aᵀλ∂q,
        ∂Aq̇∂q,
        Q,
        K,
        C,
        named,
    )
end


"""
$(TYPEDSIGNATURES)
"""
function update!(cab::DistanceSpringDamperSegment, p1, p2, ṗ1, ṗ2, s1, s2)
	(;k, c, state, prestress) = cab
    state.start .= p1     
    state.stop .= p2
    state.start_vel .= ṗ1
    state.stop_vel .= ṗ2
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


function ApparatusCache(joint::AbstractJoint,force::ClusterDistanceSpringDampers,full_coords_idx,)
    T = get_numbertype(joint)
    (;num_of_cstr,) = joint
    num_appar_full_idx = length(full_coords_idx)
    
    φ = zeros(T,num_of_cstr)
    ∂φ∂q = spzeros(T,num_of_cstr,num_appar_full_idx)
    ∂Aᵀλ∂q = spzeros(T,num_appar_full_idx,num_appar_full_idx)
    ∂Aq̇∂q = spzeros(T,num_of_cstr,num_appar_full_idx)
    Q = zeros(T,num_appar_full_idx)
    K = spzeros(T,num_appar_full_idx,num_appar_full_idx)
    C = spzeros(T,num_appar_full_idx,num_appar_full_idx)
    
    segs = force.segments
    sps = joint.sliding_points
    nsps = length(sps)
    ns = 2 * nsps
    nsegs = length(segs)
    @assert nsps + 1 == nsegs
    # n = length(sps)
    ap = [-segs[i+1].force.k for i in 1:nsps-1]
    ap0 = [segs[i].force.k + segs[i+1].force.k/sps[i].α for i in 1:nsps]
    ap1 = [-segs[i+1].force.k/sps[i].α for i in 1:nsps-1]
    an0 = [segs[i].force.k + segs[i+1].force.k*sps[i].α for i in 1:nsps]
    an1 = [-segs[i+1].force.k*sps[i].α for i in 1:nsps-1]
    A⁺ = Tridiagonal(ap, ap0, ap1)
    A⁻ = Tridiagonal(ap, an0, an1)
    A = [A⁺ -A⁺; -A⁻ A⁻]
    𝐓 = zeros(Int, 2nsps, 2nsps)
    for (j, i) in enumerate(1:2:2nsps)
        𝐓[i, j] = 1
    end
    for (j, i) in enumerate(2:2:2nsps)
        𝐓[i, j+nsps] = 1
    end
    N = zeros(T, nsegs, 2nsegs-2)
    N[1, 1:2] = [1 -1]
    N[end, end-1:end] = [-1 1]
    for i in 2:nsegs-1
        N[i, 2i-3:2i] = [-1 1 1 -1]
    end

    b⁺ = zeros(T, nsegs-1, nsegs)
    b⁻ = zeros(T, nsegs-1, nsegs)
    for i in 1:nsps
        k1 = segs[i  ].force.k
        k2 = segs[i+1].force.k
        α = sps[i].α
        b⁺[i, i:i+1] = [-k1  k2/α]
        b⁻[i, i:i+1] = [ k1 -k2*α]
    end
    b = vcat(b⁺, b⁻)
    Tb = 𝐓 * b
    ∂ζ∂s = 𝐓 * A * 𝐓'
    ζ = zeros(T, ns)
    named = @eponymtuple(
        ns,
        nsegs,
        nsps,
        ∂ζ∂s,
        Tb,
        N,
        ζ
    )
    
    cache = ApparatusCache(
        φ,
        ∂φ∂q,
        ∂Aᵀλ∂q,
        ∂Aq̇∂q,
        Q,
        K,
        C,
        named
    )
end


function execute!(structure::AbstractStructure,actuator::ExternalForceActuator{<:Apparatus{<:ClusterJoint},<:NaiveOperator},u) 
    (;id,signifier,operator,force) = actuator
    (;state) = structure
    segs = signifier.force.segments
    segs[1].force.state.restlen = segs[1].force.original_restlen - u[1]
end

function switch!(appar::Apparatus{<:ClusterJoint,},
        structure::AbstractStructure,
        q, s
    )
    (;ζ,) = appar.cache.named

    (;joint, force) = appar
    segs = force.segments
    sps = joint.sliding_points

    s⁺ = @view   s[begin:2:end]
    s⁻ = @view   s[begin+1:2:end]
    ζ⁺ = @view   ζ[begin:2:end]
    ζ⁻ = @view   ζ[begin+1:2:end]
    for (i, sp) in enumerate(sps)
        ζ⁺[i] = segs[i+1].force.state.tension/sp.α   - segs[i  ].force.state.tension
        ζ⁻[i] = segs[i  ].force.state.tension - sp.α*segs[i+1].force.state.tension
        # sp.κ₁⁺ = clamp(1/ζ⁺[i],1e-1,1e0)
        sp.κ₂⁺ = clamp(1/s⁺[i],1e-2,1e1)
        # sp.κ₁⁻ = clamp(1/ζ⁻[i],1e-1,1e0)
        sp.κ₂⁻ = clamp(1/s⁻[i],1e-2,1e1)
    end
    # @show ζ .* s
end


function get_appar_idx(appar::Apparatus{<:ClusterJoint},bodyid2sys_full_coords)
    segs = appar.force.segments
    first_hen_id = segs[1].joint.hen2egg.hen.body.prop.id
    first_egg_id = segs[1].joint.hen2egg.egg.body.prop.id
    first_hen_full = bodyid2sys_full_coords[first_hen_id]
    first_egg_full = bodyid2sys_full_coords[first_egg_id]
    appar_full_coords_idx = vcat(
        first_hen_full,
        first_egg_full
    )
    for seg in segs[2:end]
        (;id) = seg.joint.hen2egg.egg.body.prop
        full_idx = bodyid2sys_full_coords[id]
        appar_full_coords_idx = vcat(
            appar_full_coords_idx,
            full_idx
        )
    end
    unique(appar_full_coords_idx)
end

function gen_force_auxi_jacobian!(∂F∂s, appar::Apparatus{<:ClusterJoint}, st::AbstractStructure,q,q̇,t,s)
    cnt = st.connectivity
    (; apparid2sys_full_coords_idx, apparid2sys_aux_var_idx) = cnt
    (; N, nsegs) = appar.cache.named
    segs = appar.force.segments
    appar_full_coords_idx = apparid2sys_full_coords_idx[appar.id]
    kc = [seg.force.k for seg in segs]
    foreach(segs) do seg
        (; full_coords_idx, id) = seg
        (; hen, egg) = seg.joint.hen2egg
        if id == 1
            cluster_idx = collect(1:2)
        elseif id == nsegs
            cluster_idx = collect(2nsegs-3:2nsegs-2)
        else
            cluster_idx = collect(2id-3:2id)
        end
        num_of_cluster_idx = length(cluster_idx)
        num_full_coords_idx = length(full_coords_idx)
        seg_full_coords_idx = appar_full_coords_idx[full_coords_idx]
        body_hen = hen.body
        body_egg = egg.body
        T = get_numbertype(body_hen)
        ndim = get_num_of_dims(body_hen)
        J = zeros(T,ndim,num_full_coords_idx)
        lkn = zeros(T, num_of_cluster_idx, ndim)
        C_hen = body_hen.cache.Cps[hen.pid]
        C_egg = body_egg.cache.Cps[egg.pid]
        ncoords_hen = get_num_of_coords(body_hen.coords)
        ncoords_egg = get_num_of_coords(body_egg.coords)
        (; k, c, state) = seg.force
        (; direction, tension) = state
        if tension == 0
        else
            J .= 0
            J[:,1:ncoords_hen]              .-= C_hen
            J[:,ncoords_hen+1:ncoords_hen+ncoords_egg] .+= C_egg
            kN = kc[id] .* N[id, cluster_idx]
            for il = 1:ndim 
                for ik = 1:num_of_cluster_idx
                    lkn[ik, il] = direction[il] * kN[ik]
                end
            end
            ∂F∂s[seg_full_coords_idx, cluster_idx] .-= transpose(lkn * J)
        end
    end
end

get_numbertype(::ClusterJoint{SlidingPoint{T}}) where{T} = T

function auxi_function!(ret, appar::Apparatus{<:ClusterJoint,},
        cnt::AbstractConnectivity, q, s
    )
    (;ζ,) = appar.cache.named

    (;joint, force) = appar
    segs = force.segments
    sps = joint.sliding_points

    S⁺ = @view ret[begin:2:end]
    S⁻ = @view ret[begin+1:2:end]
    s⁺ = @view   s[begin:2:end]
    s⁻ = @view   s[begin+1:2:end]
    ζ⁺ = @view   ζ[begin:2:end]
    ζ⁻ = @view   ζ[begin+1:2:end]
    for (i, sp) in enumerate(sps)
        (;ϵ⁺, κ₁⁺, κ₂⁺, ϵ⁻, κ₁⁻, κ₂⁻) = sp
        ζ⁺[i] = segs[i+1].force.state.tension/sp.α   - segs[i  ].force.state.tension
        ζ⁻[i] = segs[i  ].force.state.tension - sp.α*segs[i+1].force.state.tension
        S⁺[begin-1+i] = Rible.FischerBurmeister(ϵ⁺, κ₁⁺, κ₂⁺)(ζ⁺[i], s⁺[begin-1+i])
        S⁻[begin-1+i] = Rible.FischerBurmeister(ϵ⁻, κ₁⁻, κ₂⁻)(ζ⁻[i], s⁻[begin-1+i])
    end
end

function add_tangent_stiffness_matrix!(∂Q̌∂q,appar::Apparatus{<:ClusterJoint},st::AbstractStructure,q,s=nothing)
    foreach(appar.force.segments) do seg
        add_tangent_stiffness_matrix!(∂Q̌∂q,seg,st,q)
    end
    nothing
end

function add_tangent_stiffness_matrix!(∂Q̌∂q,appar::Apparatus{<:CableJoint,<:DistanceSpringDamperSegment},st::AbstractStructure,q,s=nothing)
    
    cnt = st.connectivity
    (;
        bodyid2sys_full_coords,
        apparid2sys_full_coords_idx
    ) = cnt
    (;joint, force, cache) = appar
    T = get_numbertype(st)
    num_of_dim = get_num_of_dims(st)
    appar_full_coords_idx = get_appar_idx(appar,bodyid2sys_full_coords)
    D = @MMatrix zeros(T,num_of_dim,num_of_dim)
    Im = Symmetric(SMatrix{num_of_dim,num_of_dim}(one(T)*I))
    J = cache.named.J
    DJ = cache.named.DJ

    (;hen,egg) = joint.hen2egg
    body_hen = hen.body
    body_egg = egg.body
    C_hen = body_hen.cache.Cps[hen.pid]
    C_egg = body_egg.cache.Cps[egg.pid]
    ncoords_hen = get_num_of_coords(body_hen)
    ncoords_egg = get_num_of_coords(body_egg)
    (;k,c,state) = force
    (;direction,tension) = state
    l = state.length
    l̇ = state.lengthdot
    D .= direction*transpose(direction)
    density = tension/l
    β = c*l̇/l + density
    D .*= k-β
    D .+= β.*Im
    fill!(J, zero(T))
    J[:,1:ncoords_hen] .-= C_hen
    J[:,ncoords_hen.+1:ncoords_hen+ncoords_egg] .+= C_egg
    ## @show size(J),size(cache.K)
    mul!(DJ, D, J)
    mul!(cache.K, transpose(J), DJ, -one(T), zero(T))
    ## @tullio cache.K[i,j]  = -J[k,i]*D[k,l]*J[l,j]
    ## @show typeof(JDJ)
    ∂Q̌∂q[appar_full_coords_idx,appar_full_coords_idx] .+= cache.K
    nothing
end


function update_apparatus!(st::AbstractStructure, appar::Apparatus{<:ClusterJoint},s)
    (;id, joint, force) = appar
    segs = force.segments
    sps = joint.sliding_points
    (;nsegs) = appar.cache.named
    for (i, sp) in enumerate(sps)
        sp.s⁺ = s[2i-1]
        sp.s⁻ = s[2i  ]
        sp.s = sp.s⁺ - sp.s⁻
    end
    for (id, seg) in enumerate(segs)
        hen = seg.joint.hen2egg.hen
        egg = seg.joint.hen2egg.egg
        locus_state_hen = hen.body.state.loci_states[hen.pid]
        locus_state_egg = egg.body.state.loci_states[egg.pid]
        p_hen = locus_state_hen.frame.position
        ṗ_hen = locus_state_hen.frame.velocity
        f_hen = locus_state_hen.force
        p_egg = locus_state_egg.frame.position
        ṗ_egg = locus_state_egg.frame.velocity
        f_egg = locus_state_egg.force
        if (id==1)
            s1 = 0
            s2 = sps[id].s
        elseif (id==nsegs)
            s1 = sps[nsegs-1].s
            s2 = 0
        else
            s1 = sps[id-1].s
            s2 = sps[id].s
        end
        update!(seg.force, p_hen, p_egg, ṗ_hen, ṗ_egg, s1, s2)
        f_hen .+= seg.force.state.force
        f_egg .-= seg.force.state.force
    end
end


function auxi_jacobian!(∂S∂q,∂S∂s, appar::Apparatus{<:ClusterJoint, <:ClusterDistanceSpringDampers},
        cnt::AbstractConnectivity, q, s
    )
    (;num_of_full_coords, apparid2sys_aux_var_idx, bodyid2sys_full_coords) = cnt
    (;ns, nsegs, nsps, ∂ζ∂s, Tb, ζ,) = appar.cache.named
    (; force, joint, id) = appar
    segs = force.segments
    sps = joint.sliding_points

    T = eltype(q)
    ndim = get_num_of_dims(force.segments[begin])
        
    ζ⁺ = @view ζ[begin:2:end]
    ζ⁻ = @view ζ[begin+1:2:end]
    s⁺ = @view s[begin:2:end]
    s⁻ = @view s[begin+1:2:end]
    ∂FB∂x = deepcopy(s)
    ∂FB∂y = deepcopy(s)
    for (i, sp) in enumerate(sps)
        (;ϵ⁺, κ₁⁺, κ₂⁺, ϵ⁻, κ₁⁻, κ₂⁻) = sp
        ζ⁺[i] = segs[i+1].force.state.tension/sps[i].α - segs[i  ].force.state.tension
        ζ⁻[i] = segs[i].force.state.tension - sps[i].α * segs[i+1].force.state.tension
        ∂FB⁺ = Rible.FischerBurmeister_first_derivatives(ϵ⁺, κ₁⁺, κ₂⁺).(ζ⁺[i], s⁺[i])
        ∂FB⁻ = Rible.FischerBurmeister_first_derivatives(ϵ⁻, κ₁⁻, κ₂⁻).(ζ⁻[i], s⁻[i])
        ∂FB∂x[begin+2(i-1)]   = ∂FB⁺.∂x
        ∂FB∂x[begin+2(i-1)+1] = ∂FB⁻.∂x 
        ∂FB∂y[begin+2(i-1)]   = ∂FB⁺.∂y
        ∂FB∂y[begin+2(i-1)+1] = ∂FB⁻.∂y
    end
    ∂FB∂ζ = ∂FB∂x |> Diagonal
    ∂FB∂s = ∂FB∂y |> Diagonal

    ∂l∂q = zeros(T, nsegs, num_of_full_coords)
    J = zeros(T, ndim, num_of_full_coords)
    foreach(segs) do seg
        (; hen, egg) = seg.joint.hen2egg
        body_hen = hen.body
        body_egg = egg.body
        C_hen = body_hen.cache.Cps[hen.pid]
        C_egg = body_egg.cache.Cps[egg.pid]
        mfull_hen = bodyid2sys_full_coords[body_hen.prop.id]
        mfull_egg = bodyid2sys_full_coords[body_egg.prop.id]
        (; direction) = seg.force.state
        J .= 0
        J[:, mfull_egg] .+= C_egg
        J[:, mfull_hen] .-= C_hen
        ∂l∂q[seg.id, :] = direction'*J
    end

    ∂ζ∂q = Tb * ∂l∂q
    # @show size(∂S∂q) size(∂FB∂ζ) size(∂ζ∂q)
    ∂S∂q .= ∂FB∂ζ * ∂ζ∂q
    ∂S∂s .= ∂FB∂ζ * ∂ζ∂s + ∂FB∂s 
end
