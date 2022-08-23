abstract type ForwardMode end
struct PrimalMode <: ForwardMode end
struct StiffMode <: ForwardMode end
struct DeformMode <: ForwardMode end
struct AllMode <: ForwardMode end

function initialize_sequence(nt,len=1)
    StructArray{typeof(nt)}(undef,len)
end

function split_by_lengths(x::AbstractVector, n::AbstractVector{<:Int})
     result = Vector{Vector{eltype(x)}}()
     start = firstindex(x)
     for len in n
       push!(result, x[start:(start + len - 1)])
       start += len
     end
     result
end

function split_by_lengths(x::AbstractVector, len::Int)
     result = Vector{Vector{eltype(x)}}()
     start = firstindex(x)
     last = lastindex(x)
     while start+len-1 <= last
       push!(result, x[start:(start+len-1)])
       start += len
     end
     result
end

function check_and_retrieve(result,lens::AbstractVector{<:Int})
    # real_path_results = results(result; only_real=true,
    #                             only_nonsingular=true,
    #                             only_finite=true)
    # @assert length(real_path_results) == 1
    path_results = results(result)
    if length(path_results) != 1
        @show failed(result)
        error("Tracking failed.")
    end
    path_result1 = path_results[1]
    nstep = steps(path_result1)
    res = residual(path_result1)
    sol = real(solution(path_result1))
    q̌,s,λ = split_by_lengths(sol,lens)
    (q̌=q̌,s=s,λ=λ,nstep=nstep,res=res)
end

function check_and_retrieve(result,Psys::HomotopyContinuation.System)
    check_and_retrieve(result,length.(Psys.variable_groups))
end

function check_slackness(𝐥,𝐮)
    @assert length(𝐥) == length(𝐮)
    for (i,(l,u)) in enumerate(zip(𝐥,𝐮))
        Δl = l - u
        if Δl < 0
            @warn "The $(i)th string is slack, Δl = $(Δl)"
        end
    end
end

function forward_system(tg,mode=PrimalMode();F̌=reshape(build_Ǧ(tg),:,1))
    (;nconstraints) = tg
    (;nfree) = tg.connectivity.indexed
    (;nc) = tg.connectivity.numbered
    ns = tg.tensiles.cables |> length
    nλ = nconstraints
    @var q̌[1:nfree]
    @var s[1:ns]
    @var λ[1:nλ]
    @var d[1:nλ]
    @var c[1:nc]
    @var k[1:ns]
    @var u[1:ns]
    @var g[1:size(F̌,2)]
    polyq̌ = 1.0q̌
    polys = 1.0s
    polyλ = 1.0λ
    polyd = 1.0d
    polyc = 1.0c
    polyk = 1.0k
    polyu = 1.0u
    polyg = 1.0g
    q0 = get_q(tg)
    Φ = make_Φ(tg,q0)
    A = make_A(tg,q0)
    Q̌ = make_Q̌(tg,q0)
    S = make_S(tg,q0)

    variable_groups = [q̌,s,λ]

    scaling = 1
    if mode isa PrimalMode
        P = [transpose(A(polyq̌))*polyλ - Q̌(polyq̌,polys,polyu) - F̌*g;
            S(polyq̌,polys);
            Φ(polyq̌)]
        parameters = [u;g]
    elseif mode isa StiffMode
        P = [transpose(A(polyq̌))*polyλ - Q̌(polyq̌,polys,polyu,polyk) - F̌*g;
            S(polyq̌,polys);
            Φ(polyq̌)]
        parameters = [k;u;g]
    elseif mode isa DeformMode
        P = [transpose(A(polyq̌))*polyλ - Q̌(polyq̌,polys,polyu) - F̌*g;
            S(polyq̌,polys);
            Φ(polyq̌,polyd)]
        parameters = [d;u;g]
    elseif mode isa AllMode
        P = [transpose(A(polyq̌,polyc))*polyλ -  Q̌(polyq̌,polys,polyu,polyk,polyc) - F̌*g;
            S(polyq̌,polys,polyc);
            Φ(polyq̌,polyd,polyc)]
        parameters= [d;c;k;u;g]
    else
        error("Invalid mode")
    end
    # PPP = [subs(f, u=>l,g=>0.0) for f in P]
    # vars = (q̌=q̌,s=s,λ=λ,d=d,u=u,k=k,g=g)
    P,variable_groups,parameters
end

function forward_once(Psys::HomotopyContinuation.System,
                        startsols,start_parameters,target_parameters)
    tracker_options = TrackerOptions(;extended_precision=true,parameters=:default)
    result = HomotopyContinuation.solve(Psys, startsols; start_parameters, target_parameters, tracker_options, threading = false)
    # result = HomotopyContinuation.solve(Fsys, startsols)
    check_and_retrieve(result,Psys)
end

function find_diff_system(
            P,variable_groups,parameters,
            parameter_points
        )
    diff_parameters = copy(parameters)
    diff_parameter_points = [
        deepcopy(pp)
        for pp in parameter_points
    ]
    first_parameter_point = first(parameter_points)
    npoints = length(parameter_points)
    plens = [length(p) for p in first_parameter_point]
    pidxs = split_by_lengths(1:length(parameters),plens)
    for (i,p) in enumerate(keys(first_parameter_point))
        p_points = [a[p] for a in parameter_points]
        counts = map(diff(p_points)) do d
                d.== 0
        end |> sum
        identicals = counts .== (npoints-1).*one.(counts)
        differents = .!(identicals)
        P = subs(P,parameters[pidxs[i]][identicals]=>first_parameter_point[p][identicals])
        pidxs[i] = pidxs[i][identicals]
        for a in diff_parameter_points
            deleteat!(a[p],collect(1:plens[i])[identicals])
        end
    end
    ide_pindx = reduce(vcat,pidxs)
    deleteat!(diff_parameters,ide_pindx)
    Psys = System(P;variable_groups,parameters=diff_parameters)
    Psys,ide_pindx
end

function forward_once(tg::TensegrityStructure,startsols_input,
                        start_parameters_input,
                        target_parameters_input,
                        mode=PrimalMode();F=reshape(build_G(tg),:,1))

    P_all,variable_groups,parameters = forward_system(tg,mode;F)
    Psys, diff_parameter_points  = find_diff_system(
                P_all,variable_groups,parameters,
                start_parameters_input,
                target_parameters_input
            )
    q0,s0,λ0 = startsols_input
    startsols = [[q0; s0; λ0]]


    # forward_once(Psys,startsols,start_parameters,target_parameters)
    target_sol = forward_once(Psys,startsols,diff_start_parameters,diff_target_parameters)
    check_slackness(inv.(target_sol.s),target_parameters_input.u)
    target_sol
end

function forward_sequence(Psys::HomotopyContinuation.System,
                        startsols_input,
                        start_parameters_input,
                        target_parameters_input,
                        ide_pindx;n=1)

    start_parameters = reduce(vcat,start_parameters_input)
    target_parameters = reduce(vcat,target_parameters_input)
    diff_parameters = (target_parameters.-start_parameters)./n

    q̌0,s0,λ0 = startsols_input
    plens = [length(p) for p in start_parameters_input]
    start_point = merge((q̌=q̌0,s=s0,λ=λ0,nstep=0,res=zero(eltype(q̌0))),start_parameters_input)
    # solseq = initialize_sequence(start_point)
    solseq = StructArray([start_point])

    for k = 1:n
        pᵏ⁻¹ = start_parameters .+ (k-1).*diff_parameters
        pᵏ = start_parameters .+ k.*diff_parameters
        xᵏ⁻¹ = [[solseq[k].q̌; solseq[k].s; solseq[k].λ]]
        sol = forward_once(Psys, xᵏ⁻¹,
                    deleteat!(copy(pᵏ⁻¹),ide_pindx),
                    deleteat!(copy(pᵏ),  ide_pindx)
            )
        parameters = deepcopy(start_parameters_input)
        foreach((x,y)-> x[:] = y, parameters, split_by_lengths(pᵏ,plens))
        # check_slackness(inv.(sol.s),parameters.u)
        StructArrays.append!!(solseq,[merge(sol,parameters)])
    end
    solseq
end

function forward_sequence(tg::TensegrityStructure,startsols,
                        start_parameters,
                        target_parameters,
                        mode=PrimalMode();F̌=reshape(build_Ǧ(tg),:,1))
    P,variable_groups,parameters = forward_system(tg,mode;F̌)
    Psys = System(P;variable_groups,parameters)
    forward_sequence(Psys,startsols,start_parameters,target_parameters)
end

function forward_multi_sequence(Psys::HomotopyContinuation.System,startsols_input,
                        parameter_points, ide_pindx, mode=PrimalMode();n=1)

    # parameter_point1 = parameter_points[1]
    startsols_inputs = [startsols_input]
    [
        begin
            @debug "Forwarding the $(i)th sequence."
            seq = forward_sequence(Psys,
                                startsols_inputs[i],
                                parameter_points[i],
                                parameter_points[i+1],
                                ide_pindx;n)
            push!(startsols_inputs,(q̌=seq[end].q̌,s=seq[end].s,λ=seq[end].λ))
            seq
        end
        for (i,parameter_point) in enumerate(parameter_points[begin:end-1])
    ]
end

function forward_multi_sequence(tg::TensegrityStructure,startsols,
                        parameter_points,mode=PrimalMode();
                        F̌=reshape(build_Ǧ(tg),:,1),n=1)
    P,variable_groups,parameters = forward_system(tg,mode;F̌)
    # q̌0,s0,λ0 = startsols
    # q̌,s,λ = variable_groups
    # Pz = map(P) do z
    #     z(q̌=>q̌0,s=>s0,λ=>λ0,parameters=>reduce(vcat,parameter_points[1]))
    # end
    # Psys = System(P;variable_groups,parameters)

    Psys, ide_pindx = find_diff_system(
                P,variable_groups,parameters,
                parameter_points
            )
    seqs = forward_multi_sequence(Psys,startsols,parameter_points, ide_pindx,mode;n)
    [recover.(seq,Ref(tg)) for seq in seqs]
end

function recover(state::NamedTuple,tg::AbstractTensegrityStructure)
    q = get_q(tg)
    (;sysfree) = tg.connectivity.indexed
    q[sysfree] = state.q̌
    merge((q=q,),state)
end

function apply!(bot,seq)
    (;traj) = bot
    resize!(traj,1)
    append!(traj.q, seq.q[1:end])
    # append!(traj.c, seq.c[1:end])
    append!(traj.t, collect(1:length(seq)))
end

function get_start_sol(bot)
    (;tg) = bot
    q̌ = get_q̌(tg)
    isequ, λ = check_static_equilibrium_output_multipliers(bot.tg)
    if isequ
        @info "Alreadly in static equilibrium, skipping inverse."
        # λ = inverse_for_multipliers(bot,bot);
        u = get_cables_restlen(bot)
    else
        λ, u = inverse_for_restlength(bot,bot)
    end
    ℓ = get_cables_len(bot)
    s = inv.(ℓ)
    @eponymtuple(q̌,s,λ),u
end

function get_start_system(bot,mode=PrimalMode();F=reshape(build_Ǧ(bot.tg),:,1))
    (;tg) = bot
    start_sol,u = get_start_sol(bot)
    g = zeros(get_numbertype(bot),size(F,2))
    if typeof(mode)<:PrimalMode
        start_parameters = @eponymtuple(u,g)
    elseif typeof(mode)<:StiffMode
        k = get_cables_stiffness(tg)
        start_parameters = @eponymtuple(k,u,g)
    elseif typeof(mode)<:DeformMode
        d = get_d(tg)
        start_parameters = @eponymtuple(d,u,g)
    elseif typeof(mode)<:AllMode
        d = get_d(tg)
        c = get_c(tg)
        k = get_cables_stiffness(tg)
        start_parameters = @eponymtuple(d,c,k,u,g)
    else
        error("Invalid mode")
    end
    start_sol, start_parameters
end
