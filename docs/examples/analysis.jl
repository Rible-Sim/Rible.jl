using Test

function get_trajectory!(bot::TR.TensegrityRobot,rbid::Int,pid::Int,step_range=:)
    (; tg, traj)= bot
	T = TR.get_numbertype(bot)
    rp = Vector{T}[]
    rbs = TR.get_bodies(tg)
    rb = rbs[rbid]
    for q in traj.q
        TR.update_rigids!(tg,q)
        if pid == 0
            push!(rp,rb.state.rg)
        else
            push!(rp,rb.state.rps[pid])
        end
    end
    rp[step_range] |> VectorOfArray
end

function get_velocity!(bot::TR.TensegrityRobot,rbid::Int,pid::Int,step_range=:)
    (; tg, traj)= bot
	T = TR.get_numbertype(bot)
    ṙp = Vector{T}[]
    rbs = TR.get_bodies(tg)
    rb = rbs[rbid]
    for (q,q̇) in zip(traj.q, traj.q̇)
        TR.update_rigids!(tg,q,q̇)
        if pid == 0
            push!(ṙp,rb.state.ṙg)
        else
	        push!(ṙp,rb.state.ṙps[pid])
		end
    end
    ṙp[step_range] |> VectorOfArray
end

function get_kinetic_energy!(bot::TR.TensegrityRobot,rbid::Int,step_range=:)
    (; tg, traj)= bot
	numberType = TR.get_numbertype(bot)
    T = numberType[]
    rbs = TR.get_bodies(tg)
    rb = rbs[rbid]
    for (q,q̇) in zip(traj.q, traj.q̇)
        TR.update_rigids!(tg,q,q̇)
        Tt = TR.kinetic_energy_translation(rb)
        Tr = TR.kinetic_energy_rotation(rb)
        push!(T,Tt+Tr)
    end
    T[step_range]
end

function get_mid_velocity!(bot::TR.TensegrityRobot,rbid::Int,pid::Int,step_range=:)
    (; tg, traj)= bot
	(; t, q) = traj
	T = TR.get_numbertype(bot)
	h = t[begin+1] - t[begin]
    ṙp = Vector{T}[]
    rbs = TR.get_bodies(tg)
    rb = rbs[rbid]
    for (qₖ,qₖ₋₁) in zip(traj.q[begin+1:end], traj.q[begin:end-1])
        TR.update_rigids!(tg,(qₖ.+qₖ₋₁)./2,(qₖ.-qₖ₋₁)./h)
        push!(ṙp,rb.state.ṙps[pid])
    end
    ṙp[step_range] |> VectorOfArray
end

function get_orientation!(bot::TR.TensegrityRobot,rbid::Int,step_range=:)
    (; tg, traj)= bot
	T = TR.get_numbertype(bot)
    R = VectorOfArray(Vector{Matrix{T}}())
    rbs = TR.get_bodies(tg)
    rb = rbs[rbid]
    for (q,q̇) in zip(traj.q, traj.q̇)
        TR.update_rigids!(tg,q,q̇)
        TR.update_orientations!(tg)
        push!(R,rb.state.R)
    end
    R[step_range] |> VectorOfArray
end

function get_angular_velocity!(bot::TR.TensegrityRobot,rbid::Int,step_range=:)
    (; tg, traj)= bot
	T = TR.get_numbertype(bot)
    ω = Vector{T}[]
    rbs = TR.get_bodies(tg)
    rb = rbs[rbid]
    for (q,q̇) in zip(traj.q, traj.q̇)
        TR.update_rigids!(tg,q,q̇)
        push!(ω,rb.state.ω)
    end
    ω[step_range] |> VectorOfArray
end

function get_time_mids(bot::TR.TensegrityRobot)
    (;t) = bot.traj
    (t[begin+1:end] .+ t[begin:end-1])./2
end

function get_tension!(bot::TR.TensegrityRobot,cid::Int,step_range=:)
    (; tg, traj)= bot
    (; cables) = tg.tensiles
    T = TR.get_numbertype(tg)
    f = Vector{T}()
    h = traj.t[2] - traj.t[1]
    q_mids = [(traj.q[k] .+ traj.q[k+1])./2 for k = 1:length(traj)-1]
    q̇_mids = [(traj.q[k] .- traj.q[k+1])./h for k = 1:length(traj)-1]
    for (q,q̇) in zip(q_mids, q̇_mids)
        TR.update!(tg,q)
        push!(f,cables[cid].state.tension)
    end
    f
end

function get_actuation_traj(tr::TR.TensegrityRobot)
    vecarr_to_vectors(VectorOfArray([ctrller.traj.us for ctrller in tr.hub.actuators]))
end

function analyse_energy(tr_input;actuation=false,gravity=false,elasticity=false,factor=1)
    tr = deepcopy(tr_input)
    (;tg, traj) = tr
    kes = [TR.kinetic_energy_coords(tg,q,q̇) for (q,q̇) in zip(traj.q,traj.q̇)]
    epes = 0.0
    gpes = 0.0

    es = kes
    if elasticity
        if actuation
            as = get_actuation_traj(tr)
            epes = factor.*[TR.elastic_potential_energy(tr,q,a) for (q,a) in zip(traj.qs,as)]
        else
            epes = factor.*[TR.elastic_potential_energy(tg,q) for q in traj.qs]
        end

        es += epes
    end
    if gravity
        gpes = factor.*[TR.gravity_potential_energy(tg,q) for q in traj.qs]
        es += gpes
    end
    es1 = es[1]
    es_err = abs.((es.-es1)./es1)

    kes,epes,gpes,es,es_err
end

function string_potential(tgstruct,sol)
    [
    [begin
        TR.reset_forces!(tgstruct)
        TR.distribute_q_to_rbs!(tgstruct,q)
        TR.update_cables_apply_forces!(tgstruct)
        TR.potential_energy(s)
    end
    for q in sol.qs] for s in tgstruct.cables]
end

function analyse_slackness(tgstruct,q)
    TR.distribute_q_to_rbs!(tgstruct,q)
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
    (;tg, traj) = tr
    @testset "Test Slackness" begin
        @testset "Solution $j" for (j,q) in enumerate(traj.qs)
            TR.distribute_q_to_rbs!(tg,q)
            @testset "String $i" for (i,s) in enumerate(tg.cables)
                len = s.state.length
                restlen = s.state.restlen
                @test len > restlen
            end
        end
    end
end

function compute_contraint_error(bot)
    (;tg, traj) = bot
    (;qs, q̇s) = traj
    Φ = TR.build_Φ(tg)
    A = TR.build_A(tg)
    position_constraint_errors = [norm(Φ(q)) for q in qs]
    velocity_constraint_errors = [norm(A(q)*q̇)/norm(q̇) for (q,q̇) in zip(qs,q̇s)]
    position_constraint_errors, velocity_constraint_errors
end

function get_tangent(x)
	@assert length(x) == 3
	mag = norm(x[2:3])
	dir = Meshes.Vec(x[2]/mag,x[3]/mag)
	mag, dir
end

function isactive(contact)
	contact.state.active
end
function isimpact(contact;Λtol=1e-8) # rule out false active
	cs = contact.state
	cs.active && !cs.persistent && norm(cs.Λ) > Λtol
end
function doespersist(contact;Λtol=1e-8) # rule out false active
	cs = contact.state
	cs.active && cs.persistent && norm(cs.Λ) > Λtol
end
function issliding(contact;vtol=1e-8)
	cs = contact.state
	cs.active && norm(cs.v[2:3]) > vtol
end

function get_contact_angle(contact;Λtol=1e-7,vtol=Λtol)
	(;id,μ,e,state) = contact
	(;active,persistent,v,Λ) = state
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
