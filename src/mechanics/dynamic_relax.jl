
function make_viscous_damper(c=0.95)
    γ(θ) = c
end

function make_kinetic_damper()
    function γ(θ)
        if θ>0
            return 1
        else
            return 0
        end
    end
end

function make_drift_damper(c=0.95,d=20)
    γ(θ) = c + θ/d
end


function get_pseudo_inverse(J)
    J⁺ = transpose(J)*inv(J*transpose(J))
end

function get_project(J,J⁺)
    I-J⁺*J
end

function get_project(J)
    J⁺ = get_pseudo_inverse(J)
    get_project(J,J⁺)
end

function project_gradient!(r,ω,proj)
    r .= -proj*ω
end

function parallel_transport!(q,proj)
    projected_q = proj*q
    q .= norm(q)./norm(projected_q).*projected_q
end

function compute_θ(ϕ,r)
    θ = transpose(ϕ)*r./(norm(ϕ)*norm(r))
end

function pull_back!(inst_state,structure;N=10,ξ=1e-14)
    for s = 1:N
        b = cstr_function(structure,inst_state)
        normb = norm(b,Inf)
        if normb < ξ
            return normb, s
        elseif s == N
            @warn "Max iternation reached for pull back"
            return normb, s
        end
        J = cstr_jacobian(structure,inst_state)
        J⁺ = get_pseudo_inverse(J)
        inst_state.q̌ .-= J⁺*b
    end
end

function GDR!(
        bot,
        field;
        β=1e-3,
        maxiters=Int(1e4),
        res=1e-7,
        N=10,
        ξ=1e-7,
        verbose=false,
        damper=Val(:kinetic),
        c=0.95,d=20
    )
    reset!(bot)
    (;structure,traj) = bot
    inst_state = deepcopy(structure.state.system)
    if damper isa Val{:kinetic}
        𝛄 = make_kinetic_damper()
    elseif damper isa Val{:viscous}
        𝛄 = make_viscous_damper(c)
    elseif damper isa Val{:drift}
        𝛄 = make_drift_damper(c,d)
    end
    t = zero(β)
    r = one.(inst_state.q̌)
    q = β.*r
    rs = Vector{eltype(r)}([Inf])
    bs = Vector{eltype(r)}([Inf])
    ss = Int[]
    for itr = 1:maxiters
        J = cstr_jacobian(structure,inst_state)
        gen_force!(inst_state, bot, field, NoPolicy())
        ω = - inst_state.F̌
        proj = get_project(J)
        project_gradient!(r,ω,proj)
        parallel_transport!(q,proj)
        θ = compute_θ(q,r)
        γ = 𝛄(θ)
        q .= γ.*q + β.*r
        inst_state.q̌ .+= β.*q
        normb, spull = pull_back!(inst_state,structure;N,ξ)
        push!(bs,normb)
        push!(ss,spull)

        t += β
        normr = norm(r,Inf)
        push!(rs,normr)
        push!(traj,deepcopy(traj[end]))
        traj.t[end] = t
        traj.q[end] .= inst_state.q
        # @show itr,normr
        if normr < res
            break
        elseif itr == maxiters
            @warn "Max iternation reached for GDR. Res = $normr"

        elseif rs[end] > rs[end-1]
            # @error("Diverging")
            # break
        end
        if verbose
            println("itr: $itr, normr: $normr")
        end
    end
    rs,bs,ss,bot
    bot
end
