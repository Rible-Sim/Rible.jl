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
    real_path_results = results(result; only_real=true,
                                only_nonsingular=true,
                                only_finite=true)
    @assert length(real_path_results) == 1
    path_result1 = real_path_results[1]
    @assert is_success(path_result1)
    return_code = path_result1.return_code
    nstep = steps(path_result1)
    sol = real(solution(path_result1))
    q,s,λ = split_by_lengths(sol,lens)
    (q=q,s=s,λ=λ,nstep=nstep)
end

function check_and_retrieve(result,Psys::HomotopyContinuation.System)
    check_and_retrieve(result,length.(Psys.variable_groups))
end

function forward_system(tg,mode=PrimalMode();F=reshape(build_G!(tg),:,1))
    @var q[1:tg.ncoords]
    @var s[1:tg.nstrings]
    @var λ[1:tg.nconstraint]
    @var d[1:tg.nconstraint]
    @var k[1:tg.nstrings]
    @var u[1:tg.nstrings]
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

    vg = [q,s,λ]

    if typeof(mode)<:PrimalMode
        P = [transpose(A(polyq))*polyλ - Q̃*U(polys,polyu)*polyq - F*g;
            S(polyq,polys);
            Φ(polyq)]
        Psys = System(P; variable_groups=vg, parameters = [u;g])
    elseif typeof(mode)<:StiffMode
        P = [transpose(A(polyq))*polyλ - Q̃*U(polys,polyu,polyk)*polyq - F*g;
            S(polyq,polys);
            Φ(polyq)]
        Psys = System(P; variable_groups=vg, parameters = [k;u;g])
    elseif typeof(mode)<:DeformMode
        P = [transpose(A(polyq))*polyλ - Q̃*U(polys,polyu)*polyq - F*g;
            S(polyq,polys);
            Φ(polyq,polyd)]
        Psys = System(P; variable_groups=vg, parameters = [d;u;g])
    elseif typeof(mode)<:AllMode
        P = [transpose(A(polyq))*polyλ - Q̃*U(polys,polyu,polyk)*polyq - F*g;
            S(polyq,polys);
            Φ(polyq,polyd)]
        Psys = System(P; variable_groups=vg, parameters = [d;k;u;g])
    else
        error("Invalid mode")
    end
    # PPP = [subs(f, u=>l,g=>0.0) for f in P]
    vars = (q=q,s=s,λ=λ,d=d,u=u,k=k,g=g)
    Psys,vars
end

function forward_once(Psys::HomotopyContinuation.System,
                        startsols,start_parameters,target_parameters)

    result = HomotopyContinuation.solve(Psys, startsols; start_parameters, target_parameters, threading = false)
    # result = HomotopyContinuation.solve(Fsys, startsols)
    check_and_retrieve(result,Psys)
end

function forward_once(tg::TensegrityStructure,startsols_input,
                        start_parameters_input,
                        target_parameters_input,
                        mode=PrimalMode();F=reshape(build_G!(tg),:,1))

    # @show maximum(abs.(to_number.(subs(P, q=>q0, s=>s0, λ=>λ0, u=>u0, g=>g0))))
    # reset_forces!(tg)
    # update_strings_apply_forces!(tg)
    # l = [s.state.length for s = tg.strings]
    # u1 = u0
    # g0 = 1.0
    # g1 = 0.0
    Psys,_ = forward_system(tg,mode;F)
    q0,s0,λ0 = startsols_input
    startsols = [[q0; s0; λ0]]
    start_parameters = reduce(vcat,start_parameters_input)
    target_parameters = reduce(vcat,target_parameters_input)

    forward_once(Psys,startsols,start_parameters,target_parameters)
end

function forward_sequence(Psys::HomotopyContinuation.System,
                        startsols_input,
                        start_parameters_input,
                        target_parameters_input;n=10)

    start_parameters = reduce(vcat,start_parameters_input)
    target_parameters = reduce(vcat,target_parameters_input)
    diff_parameters = (target_parameters.-start_parameters)./n

    q0,s0,λ0 = startsols_input
    plens = [length(p) for p in start_parameters_input]
    start_point = merge((q=q0,s=s0,λ=λ0,nstep=0),start_parameters_input)
    # solseq = initialize_sequence(start_point)
    solseq = StructArray([start_point])

    for k = 1:n
        pᵏ⁻¹ = start_parameters .+ (k-1).*diff_parameters
        pᵏ = start_parameters .+ k.*diff_parameters
        xᵏ⁻¹ = [[solseq[k].q; solseq[k].s; solseq[k].λ]]
        sol = forward_once(Psys,xᵏ⁻¹,pᵏ⁻¹,pᵏ)
        parameters = deepcopy(start_parameters_input)
        foreach((x,y)-> x[:] = y, parameters, split_by_lengths(pᵏ,plens))
        StructArrays.append!!(solseq,[merge(sol,parameters)])
    end
    solseq
end

function forward_sequence(tg::TensegrityStructure,startsols,
                        start_parameters,
                        target_parameters,
                        mode=PrimalMode();F=reshape(build_G!(tg),:,1))
    Psys,_ = forward_system(tg,mode;F)
    forward_sequence(Psys,startsols,start_parameters,target_parameters)
end

function forward_multi_sequence(Psys,startsols_input,
                        parameter_points,mode=PrimalMode();n=10)

    parameter_point1 = parameter_points[1]
    startsols_inputs = [startsols_input]
    [begin
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

function get_start_sol(bot,Y)
    @unpack tg = bot
    q,_ = get_q(tg)
    λ,u,_ = inverse(bot,deepcopy(bot),Y)
    ℓ = get_strings_len(bot)
    s = inv.(ℓ)
    (q=q,s=s,λ=λ),u
end

function get_start_system(bot,Y,mode=PrimalMode())
    @unpack tg = bot
    start_sol,u = get_start_sol(bot,Y)
    g = [zero(get_numbertype(bot))]
    if typeof(mode)<:PrimalMode
        start_parameters = (u=u,g=g)
    elseif typeof(mode)<:StiffMode
        k = get_strings_stiffness(tg)
        start_parameters = (k=k,u=u,g=g)
    elseif typeof(mode)<:DeformMode
        d = get_d(tg)
        start_parameters = (d=d,u=u,g=g)
    elseif typeof(mode)<:AllMode
        d = get_d(tg)
        k = get_strings_stiffness(tg)
        start_parameters = (d=d,k=k,u=u,g=g)
    else
        error("Invalid mode")
    end
    start_sol, start_parameters
end
