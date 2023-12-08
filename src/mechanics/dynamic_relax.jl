
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

function initialize_GDR(st,F::Nothing;gravity=true)
    Q̃ = build_Q̃(st)
    Γ = build_Γ(st)
    # 𝛚(x) =
    function 𝛚(x)
        # Q = Q̃*Γ(x)
        clear_forces!(st)
        update_bodies!(st,x)
        update_tensiles!(st)
        if gravity
            apply_gravity!(st)
        end
        F = assemble_forces!(st)
        # @show abs.(F-Q) |> maximum
        -F
    end
    𝐛 = make_cstr_function(st)
    𝐉 = make_cstr_jacobian(st)
    x0 = st.state.system.q
    x̌0 = st.state.system.q̌
    x0,x̌0,𝛚,𝐛,𝐉
end

function initialize_GDR(st,F)
    x0,x̌0,𝛚_,𝐛,𝐉 = initialize_GDR(st,nothing)
    𝛚(x) = - 𝛚_(x) - F
    x0,x̌0,𝛚,𝐛,𝐉
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

function pull_back!(x,x̌,𝐛,𝐉;N=10,ξ=1e-14)
    for s = 1:N
        b = 𝐛(x)
        normb = norm(b,Inf)
        if normb < ξ
            return normb, s
        elseif s == N
            @warn "Max iternation reached for pull back"
            return normb, s
        end
        J = 𝐉(x)
        J⁺ = get_pseudo_inverse(J)
        x̌ .-= J⁺*b
    end
end

function GDR!(
        bot,
        F=nothing;
        gravity=true,
        β=1e-3,
        maxiters=Int(1e4),
        res=1e-7,
        N=10,
        ξ=1e-7,
        verbose=false,
    )
    reset!(bot)
    (;st,traj) = bot
    x,x̌,𝛚,𝐛,𝐉 = initialize_GDR(st,F;gravity)
    # 𝛄 = make_viscous_damper()
    𝛄 = make_kinetic_damper()
    # 𝛄 = make_drift_damper()
    t = zero(β)
    r = one.(x̌)
    q = β.*r
    rs = Vector{eltype(r)}([Inf])
    bs = Vector{eltype(r)}([Inf])
    ss = Int[]
    for itr = 1:maxiters
        J = 𝐉(x)
        ω = 𝛚(x)
        proj = get_project(J)
        project_gradient!(r,ω,proj)
        parallel_transport!(q,proj)
        θ = compute_θ(q,r)
        γ = 𝛄(θ)
        q .= γ.*q + β.*r
        x̌ .+= β.*q
        normb, spull = pull_back!(x,x̌,𝐛,𝐉;N,ξ)
        push!(bs,normb)
        push!(ss,spull)

        t += β
        normr = norm(r,Inf)
        push!(rs,normr)
        push!(traj,deepcopy(traj[end]))
        traj.t[end] = t
        traj.q[end] .= x
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
