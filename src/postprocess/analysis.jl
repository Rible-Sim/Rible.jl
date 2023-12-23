function get_trajectory!(bot::Robot,bodyid::Int,pid::Int,step_range=:)
    (; structure, traj)= bot
	T = get_numbertype(bot)
    N = get_num_of_dims(bot)
    position_traj = [
        zeros(T,N)
        for _ in eachindex(traj)
    ]
    bodies = get_bodies(structure)
    body = bodies[bodyid]
    for (i,state) in enumerate(traj)
        (;t,q,q̇) = state 
        update_bodies!(structure,q,q̇)
        if pid == 0
            position_traj[i] .= body.state.mass_locus_state.frame.position
        else
	        position_traj[i] .= body.state.loci_states[pid].frame.position
		end
    end
    position_traj[step_range] |> VectorOfArray
end

function get_velocity!(bot::Robot,bodyid::Int,pid::Int,step_range=:)
    (; structure, traj)= bot
	T = get_numbertype(bot)
    N = get_num_of_dims(bot)
    velocity_traj = [
        zeros(T,N)
        for _ in eachindex(traj)
    ]
    bodies = get_bodies(structure)
    body = bodies[bodyid]
    for (i,state) in enumerate(traj)
        (;t,q,q̇) = state 
        update_bodies!(structure,q,q̇)
        if pid == 0
            velocity_traj[i] .= body.state.mass_locus_state.frame.velocity
        else
	        velocity_traj[i] .= body.state.loci_states[pid].frame.velocity
		end
    end
    velocity_traj[step_range] |> VectorOfArray
end

function get_mid_velocity!(bot::Robot,bodyid::Int,pid::Int,step_range=:)
    (; structure, traj)= bot
	(; t, q) = traj
	h = t[begin+1] - t[begin]
	T = get_numbertype(bot)
    N = get_num_of_dims(bot)
    velocity_traj = [
        zeros(T,N)
        for _ in eachindex(traj)
    ]
    bodies = get_bodies(structure)
    body = bodies[bodyid]
    for (i,(qₖ,qₖ₋₁)) in enumerate(zip(q[begin+1:end], q[begin:end-1]))
        update_bodies!(structure,(qₖ.+qₖ₋₁)./2,(qₖ.-qₖ₋₁)./h)
        if pid == 0
            velocity_traj[i] .= body.state.mass_locus_state.frame.velocity
        else
	        velocity_traj[i] .= body.state.loci_states[pid].frame.velocity
		end
    end
    velocity_traj[begin:end-1][step_range] |> VectorOfArray
end

function get_mid_angular_velocity!(bot::Robot,bodyid::Int,step_range=:)
    (; structure, traj)= bot
	(; t, q) = traj
	h = t[begin+1] - t[begin]
	T = get_numbertype(bot)
    N = get_num_of_dims(bot)
    M = 2N-3
    angular_velocity_traj = [
        zeros(T,M)
        for _ in eachindex(traj)
    ]
    bodies = get_bodies(structure)
    body = bodies[bodyid]
    for (i,(qₖ,qₖ₋₁)) in enumerate(zip(q[begin+1:end], q[begin:end-1]))
        update_bodies!(structure,(qₖ.+qₖ₋₁)./2,(qₖ.-qₖ₋₁)./h)
        angular_velocity_traj[i] .= body.state.origin_frame.angular_velocity
    end
    angular_velocity_traj[begin:end-1][step_range] |> VectorOfArray
end

function get_mid_orientation!(bot::Robot,bodyid::Int,step_range=:)
    (; structure, traj)= bot
	(; t, q) = traj
	h = t[begin+1] - t[begin]
	T = get_numbertype(bot)
    N = get_num_of_dims(bot)
    M = 2N-3
    angles_traj = [
        zeros(T,M)
        for _ in eachindex(traj)
    ]
    bodies = get_bodies(structure)
    body = bodies[bodyid]
    for (i,(qₖ,qₖ₋₁)) in enumerate(zip(q[begin+1:end], q[begin:end-1]))
        update_bodies!(structure,(qₖ.+qₖ₋₁)./2,(qₖ.-qₖ₋₁)./h)
        angles_traj[i] .= body.state.origin_frame.axes.X |> rotation2angles
    end
    angles_traj[begin:end-1][step_range] |> VectorOfArray
end

function get_orientation!(bot::Robot,bodyid::Int,step_range=:)
    (; structure, traj)= bot
	(; t, q) = traj
	T = get_numbertype(bot)
    N = get_num_of_dims(bot)
    M = 2N-3
    angles_traj = [
        zeros(T,M)
        for _ in eachindex(traj)
    ]
    bodies = get_bodies(structure)
    body = bodies[bodyid]
    for (i,qₖ) in enumerate(q)
        update_bodies!(structure,qₖ)
        angles_traj[i] .= body.state.origin_frame.axes.X |> rotation2angles
    end
    angles_traj[step_range] |> VectorOfArray
end

function get_angular_velocity!(bot::Robot,bodyid::Int,step_range=:)
    (; structure, traj)= bot
	(; t, q, q̇) = traj
	T = get_numbertype(bot)
    N = get_num_of_dims(bot)
    M = 2N-3
    angular_velocity_traj = [
        zeros(T,M)
        for _ in eachindex(traj)
    ]
    bodies = get_bodies(structure)
    body = bodies[bodyid]
    for (i,(qₖ,q̇ₖ)) in enumerate(zip(q, q̇))
        update_bodies!(structure,qₖ,q̇ₖ)
        angular_velocity_traj[i] .= body.state.origin_frame.angular_velocity
    end
    angular_velocity_traj[step_range] |> VectorOfArray
end

function get_mid_times(bot::Robot)
    (;t) = bot.traj
    (t[begin+1:end] .+ t[begin:end-1])./2
end

function get_tension!(bot::Robot,cid::Int,step_range=:)
    (; structure, traj)= bot
    (; cables) = structure.apparatuses
    T = get_numbertype(structure)
    f = Vector{T}()
    h = traj.t[2] - traj.t[1]
    q_mids = [(traj.q[k] .+ traj.q[k+1])./2 for k = 1:length(traj)-1]
    q̇_mids = [(traj.q[k] .- traj.q[k+1])./h for k = 1:length(traj)-1]
    for (q,q̇) in zip(q_mids, q̇_mids)
        update!(structure,q)
        push!(f,cables[cid].state.tension)
    end
    f
end

function get_actuation_traj(tr::Robot)
    vecarr_to_vectors(VectorOfArray([ctrller.traj.us for ctrller in tr.hub.actuators]))
end

function analyse_energy(tr_input;actuation=false,gravity=false,elasticity=false,factor=1)
    tr = deepcopy(tr_input)
    (;structure, traj) = tr
    kes = [kinetic_energy_coords(structure,q,q̇) for (q,q̇) in zip(traj.q,traj.q̇)]
    epes = 0.0
    gpes = 0.0

    es = kes
    if elasticity
        if actuation
            as = get_actuation_traj(tr)
            epes = factor.*[elastic_potential_energy(tr,q,a) for (q,a) in zip(traj.qs,as)]
        else
            epes = factor.*[elastic_potential_energy(structure,q) for q in traj.qs]
        end

        es += epes
    end
    if gravity
        gpes = factor.*[gravity_potential_energy(structure,q) for q in traj.qs]
        es += gpes
    end
    es1 = es[1]
    es_err = abs.((es.-es1)./es1)

    kes,epes,gpes,es,es_err
end

function string_potential(tgstruct,sol)
    [
    [begin
        reset_forces!(tgstruct)
        distribute_q_to_rbs!(tgstruct,q)
        update_cables_apply_forces!(tgstruct)
        potential_energy(s)
    end
    for q in sol.qs] for s in tgstruct.cables]
end

function analyse_slackness(tgstruct,q)
    distribute_q_to_rbs!(tgstruct,q)
    for (i,s) in enumerate(tgstruct.cables)
        len = s.state.length
        restlen = s.state.restlen
        Δl = len - restlen
        if Δl < 0
            @info "String $i is slack."
        else
            @info "Slackness check pass. Δl_$i = $Δl."

        end
    end
end

function test_slackness(tr)
    (;structure, traj) = tr
     "Test Slackness" 
     begin
         for (j,q) in enumerate(traj.qs)
            "Solution $j" 
            distribute_q_to_rbs!(structure,q)
             for (i,s) in enumerate(structure.cables)
                "String $i" 
                len = s.state.length
                restlen = s.state.restlen
                @assert  len > restlen
            end
        end
    end
end

function compute_contraint_error(bot)
    (;structure, traj) = bot
    (;qs, q̇s) = traj
    Φ = build_Φ(structure)
    A = build_A(structure)
    position_cstr_errors = [norm(Φ(q)) for q in qs]
    velocity_cstr_errors = [norm(A(q)*q̇)/norm(q̇) for (q,q̇) in zip(qs,q̇s)]
    position_cstr_errors, velocity_cstr_errors
end

function get_tangent(x)
	@assert length(x) == 3
	mag = norm(x[2:3])
	dir = [x[2]/mag,x[3]/mag]
	mag, dir
end

function isactive(contact)
	contact.state.active
end
function isimpact(contact;Λtol=1e-8) # rule out false active
	cs = contact.state
	cs.active && !cs.persistent && norm(cs.force) > Λtol
end
function doespersist(contact;Λtol=1e-8) # rule out false active
	cs = contact.state
	cs.active && cs.persistent && norm(cs.force) > Λtol
end
function issliding(contact;vtol=1e-8)
	cs = contact.state
	cs.active && norm(cs.v[2:3]) > vtol
end

function get_contact_angle(contact;Λtol=1e-7,vtol=Λtol)
	(;id,μ,e,state) = contact
	(;active,persistent) = state
    Λ = state.force
    v = state.relative_velocity
	if active
		vₙ = v[1]
		Λₙ = Λ[1]
		vₜ, vₜ_dir = get_tangent(v)
		Λₜ, Λₜ_dir = get_tangent(Λ)
		α_v = atan(v[3],v[2])
		α_Λ = atan(Λ[3],Λ[2])
		# δΛ = Λₙ - Λₜ
		# δα = abs(Meshes.∠(vₜ_dir,Λₜ_dir)) - π
		δα = α_v - rem2pi(α_Λ-π,RoundNearest)
		# if !persistent
			if Λₜ > Λtol &&  vₜ > vtol
			# 	return δα
				return acos(clamp((v[2:3] ⋅ Λ[2:3])/(Λₜ*vₜ),-1,1))
			else
				return missing
			end
		# end
	end
	missing
end

function get_friction_direction(contact;Λtol=1e-8)
	(;id,μ,e,state) = contact
	(;active,persistent,v,Λ) = state
	if active
		Λₙ = Λ[1]
		Λₜ, Λₜ_dir = get_tangent(Λ)
		if Λₜ > Λtol
			θ = atan(Λₙ,Λₜ)
			return θ
		end
	end
	missing
end

function check_Coulomb(it,contact;tol=1e-3,report=false)
	(;id,μ,e,state) = contact
	(;active,persistent,v,Λ) = state
	if active
		vₙ = v[1]
		Λₙ = Λ[1]
		vₜ, vₜ_dir = get_tangent(v)
		Λₜ, Λₜ_dir = get_tangent(Λ)
		δΛ = Λₙ - Λₜ
		δα = (abs(Meshes.∠(vₜ_dir,Λₜ_dir)) - π) |> abs
		if vₙ > tol
			@info "Taking-off"
			@show δα, δΛ
		elseif vₙ > -tol && vₜ > tol
			@info "Sliding"
			@show δα, δΛ
		elseif vₜ < tol
			@info "Sticking"
			if δΛ < 0
				@warn "δΛ = $δΛ < 0"
			end
		else
			@warn "Penentrating!"
			# @assert isapprox(α,π,rtol=0.1)
		end
		vₙᵏ⁻¹ = -vₙ/e
		𝝂 = [vₙ+vₜ,v[2],v[3]]
		if !persistent
			𝝂[1] += e*vₙᵏ⁻¹
		end
		@info "it = $it, persistent = $persistent, Λₙ=$(Λₙ), Λₜ=$(Λₜ), vₙ=$(vₙ), vₜ=$(vₜ), 𝝂⋅Λ=$(𝝂⋅Λ)"
	end
	nothing
end
