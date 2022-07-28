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
    q,s,λ = split_by_lengths(sol,lens)
    (q=q,s=s,λ=λ,nstep=nstep,res=res)
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

function forward_system(tg,mode=PrimalMode();F=reshape(build_G(tg),:,1))
    @var q[1:tg.ncoords]
    @var s[1:tg.ncables]
    @var λ[1:tg.nconstraint]
    @var d[1:tg.nconstraint]
    @var k[1:tg.ncables]
    @var u[1:tg.ncables]
    @var g[1:size(F,2)]
    polyq = 1.0q
    polys = 1.0s
    polyλ = 1.0λ
    polyd = 1.0d
    polyk = 1.0k
    polyu = 1.0u
    polyg = 1.0g

    Φ = build_Φ(tg)
    A = build_A(tg)
    Q̃ = build_Q̃(tg)
    U = build_U(tg)
    S = build_S(tg)

    variable_groups = [q,s,λ]

    scaling = 1
    if mode isa PrimalMode
        P = [transpose(A(polyq))*polyλ - Q̃*U(polys,polyu)*polyq - F*g;
            S(polyq,polys);
            Φ(polyq)]
        parameters = [u;g]
        Psys = System(P; variable_groups, parameters)
    elseif mode isa StiffMode
        P = [transpose(A(polyq))*polyλ - Q̃*U(polys,polyu,polyk)*polyq - F*g;
            S(polyq,polys);
            Φ(polyq)]
        parameters = [k;u;g]
        Psys = System(P; variable_groups, parameters)
    elseif mode isa DeformMode
        P = [transpose(A(polyq))*polyλ - Q̃*U(polys,polyu)*polyq - F*g;
            S(polyq,polys);
            Φ(polyq,polyd)]
        parameters = [d;u;g]
        Psys = System(P; variable_groups, parameters)
    elseif mode isa AllMode
        P = [transpose(A(polyq))*polyλ - Q̃*U(polys,polyu,polyk)*polyq - F*g;
            S(polyq,polys);
            Φ(polyq,polyd)]
        parameters= [d;k;u;g]
        Psys = System(P; variable_groups, parameters)
    else
        error("Invalid mode")
    end
    # PPP = [subs(f, u=>l,g=>0.0) for f in P]
    # vars = (q=q,s=s,λ=λ,d=d,u=u,k=k,g=g)
    P,variable_groups,parameters
end

function forward_once(Psys::HomotopyContinuation.System,
                        startsols,start_parameters,target_parameters)
    tracker_options = TrackerOptions(;extended_precision=true,parameters=:default)
    result = HomotopyContinuation.solve(Psys, startsols; start_parameters, target_parameters, tracker_options, threading = false)
    # result = HomotopyContinuation.solve(Fsys, startsols)
    check_and_retrieve(result,Psys)
end

function forward_once(tg::TensegrityStructure,startsols_input,
                        start_parameters_input,
                        target_parameters_input,
                        mode=PrimalMode();F=reshape(build_G(tg),:,1))

    P_all,variable_groups,parameters = forward_system(tg,mode;F)
    # Psys = System(P;variable_groups,parameters)
    q0,s0,λ0 = startsols_input
    startsols = [[q0; s0; λ0]]
    start_parameters = reduce(vcat,start_parameters_input)
    target_parameters = reduce(vcat,target_parameters_input)
    identicals = start_parameters .== target_parameters
    differents = .!(identicals)
    # @show subs(P_all,reduce(vcat,variable_groups)=>startsols[1],parameters=>start_parameters)
    diff_parameters = parameters[differents]
    if isempty(diff_parameters)
        @warn("Identical start and target parameters.")
        P = P_all
        diff_start_parameters = start_parameters
        diff_target_parameters = target_parameters
        Psys = System(P;variable_groups,parameters=parameters)
    else
        P = subs(P_all,parameters[identicals]=>start_parameters[identicals])
        diff_start_parameters = start_parameters[differents]
        diff_target_parameters = target_parameters[differents]
        Psys = System(P;variable_groups,parameters=diff_parameters)
    end


    # forward_once(Psys,startsols,start_parameters,target_parameters)
    target_sol = forward_once(Psys,startsols,diff_start_parameters,diff_target_parameters)
    check_slackness(inv.(target_sol.s),target_parameters_input.u)
    target_sol
end

function forward_sequence(Psys::HomotopyContinuation.System,
                        startsols_input,
                        start_parameters_input,
                        target_parameters_input;n=1)

    start_parameters = reduce(vcat,start_parameters_input)
    target_parameters = reduce(vcat,target_parameters_input)
    diff_parameters = (target_parameters.-start_parameters)./n

    q0,s0,λ0 = startsols_input
    plens = [length(p) for p in start_parameters_input]
    start_point = merge((q=q0,s=s0,λ=λ0,nstep=0,res=zero(eltype(q0))),start_parameters_input)
    # solseq = initialize_sequence(start_point)
    solseq = StructArray([start_point])

    for k = 1:n
        pᵏ⁻¹ = start_parameters .+ (k-1).*diff_parameters
        pᵏ = start_parameters .+ k.*diff_parameters
        xᵏ⁻¹ = [[solseq[k].q; solseq[k].s; solseq[k].λ]]
        sol = forward_once(Psys,xᵏ⁻¹,pᵏ⁻¹,pᵏ)
        parameters = deepcopy(start_parameters_input)
        foreach((x,y)-> x[:] = y, parameters, split_by_lengths(pᵏ,plens))
        check_slackness(inv.(sol.s),parameters.u)
        StructArrays.append!!(solseq,[merge(sol,parameters)])
    end
    solseq
end

function forward_sequence(tg::TensegrityStructure,startsols,
                        start_parameters,
                        target_parameters,
                        mode=PrimalMode();F=reshape(build_G(tg),:,1))
    P,variable_groups,parameters = forward_system(tg,mode;F)
    Psys = System(P;variable_groups,parameters)
    forward_sequence(Psys,startsols,start_parameters,target_parameters)
end

function forward_multi_sequence(Psys::HomotopyContinuation.System,startsols_input,
                        parameter_points,mode=PrimalMode();n=1)

    parameter_point1 = parameter_points[1]
    startsols_inputs = [startsols_input]
    [begin
        @info "Forwarding the $(i)th sequence."
        seq = forward_sequence(Psys,
                            startsols_inputs[i],
                            parameter_points[i],
                            parameter_points[i+1];n)
        push!(startsols_inputs,(q=seq[end].q,s=seq[end].s,λ=seq[end].λ))
        seq
    end
    for (i,parameter_point) in enumerate(parameter_points[begin:end-1])
    ]
end

function forward_multi_sequence(tg::TensegrityStructure,startsols,
                        parameter_points,mode=PrimalMode();
                        F=reshape(build_G(tg),:,1),n=1)
    P,variable_groups,parameters = forward_system(tg,mode;F)
    Psys = System(P;variable_groups,parameters)
    forward_multi_sequence(Psys,startsols,parameter_points,mode;n)
end

function get_start_sol(bot)
    @unpack tg = bot
    q,_ = get_q(tg)
    λ,u = inverse_for_restlength(bot,bot)
    # λ = inverse_for_multipliers(bot,bot); u = get_cables_restlen(bot)
    ℓ = get_cables_len(bot)
    s = inv.(ℓ)
    (q=q,s=s,λ=λ),u
end

function get_start_system(bot,mode=PrimalMode();F=reshape(build_G(bot.tg),:,1))
    @unpack tg = bot
    start_sol,u = get_start_sol(bot)
    g = zeros(get_numbertype(bot),size(F,2))
    if typeof(mode)<:PrimalMode
        start_parameters = (u=u,g=g)
    elseif typeof(mode)<:StiffMode
        k = get_cables_stiffness(tg)
        start_parameters = (k=k,u=u,g=g)
    elseif typeof(mode)<:DeformMode
        d = get_d(tg)
        start_parameters = (d=d,u=u,g=g)
    elseif typeof(mode)<:AllMode
        d = get_d(tg)
        k = get_cables_stiffness(tg)
        start_parameters = (d=d,k=k,u=u,g=g)
    else
        error("Invalid mode")
    end
    start_sol, start_parameters
end
