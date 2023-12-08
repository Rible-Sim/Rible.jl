
function undamped_eigen!(bot::Robot;gravity=false,scaling=0.01)
    (;st,traj) = bot
    q̌ = get_free_coords(st)
    ω²,δq̌ = undamped_eigen(st;gravity)
    neg_idx = findall(ω².<=0)
    if !isempty(neg_idx)
        @warn "Negative ω² occurs, idx $neg_idx, zeroing."
        ω²[neg_idx] .= 0
    end
    ω = sqrt.(ω²)
    resize!(traj,1)
    nω = length(ω)
    for i = 1:nω
        push!(traj,deepcopy(traj[end]))
        traj.t[end] = ω[i]
        δq̌i = δq̌[i]
        ratio = norm(δq̌i)/norm(q̌)
        traj.q̌[end] .= q̌ .+ scaling.*δq̌i/ratio
    end
    bot
end