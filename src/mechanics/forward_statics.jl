abstract type ForwardMode end
struct PrimalMode <: ForwardMode end
struct StiffMode <: ForwardMode end
struct DeformMode <: ForwardMode end
struct AllMode <: ForwardMode end

function initialize_sequence(nt,len=1)
    StructArray{typeof(nt)}(undef,len)
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

function forward_system(st,mode=PrimalMode();F̌=reshape(build_Ǧ(st),:,1))
    (;num_of_cstr) = st
    (;num_of_free_coords) = st.connectivity.indexed
    (;nc) = st.connectivity.numbered
    ns = st.tensiles.cables |> length
    nλ = num_of_cstr
    @polyvar q̌[1:num_of_free_coords]
    @polyvar s[1:ns]
    @polyvar λ[1:nλ]
    @polyvar d[1:nλ]
    @polyvar c[1:nc]
    @polyvar k[1:ns]
    @polyvar u[1:ns]
    @polyvar g[1:size(F̌,2)]
    polyq̌ = 1.0q̌ .+ 0.0
    polys = 1.0s .+ 0.0
    polyλ = 1.0λ .+ 0.0
    polyd = 1.0d .+ 0.0
    polyc = 1.0c .+ 0.0
    polyk = 1.0k .+ 0.0
    polyu = 1.0u .+ 0.0
    polyg = 1.0g .+ 0.0
    q0 = get_coords(st)
    Φ = make_cstr_function(st,q0)
    A = make_cstr_jacobian(st,q0)
    Q̌ = make_Q̌(st,q0)
    S = make_S(st,q0)

    variable_groups = [q̌,s,λ]
    var_lens = length.(variable_groups)
    scaling = 1
    if mode isa PrimalMode
        P = [-transpose(A(polyq̌))*polyλ - Q̌(polyq̌,polys,polyu) - F̌*g;
            S(polyq̌,polys);
            Φ(polyq̌)]
        parameters = [u;g]
    elseif mode isa StiffMode
        P = [-transpose(A(polyq̌))*polyλ - Q̌(polyq̌,polys,polyu,polyk) - F̌*g;
            S(polyq̌,polys);
            Φ(polyq̌)]
        parameters = [k;u;g]
    elseif mode isa DeformMode
        P = [-transpose(A(polyq̌))*polyλ - Q̌(polyq̌,polys,polyu) - F̌*g;
            S(polyq̌,polys);
            Φ(polyq̌,polyd)]
        parameters = [d;u;g]
    elseif mode isa AllMode
        P = [-transpose(A(polyq̌,polyc))*polyλ -  Q̌(polyq̌,polys,polyu,polyk,polyc) - F̌*g;
            S(polyq̌,polys,polyc);
            Φ(polyq̌,polyd,polyc)]
        parameters= [d;c;k;u;g]
    else
        error("Invalid mode")
    end
    # PPP = [subs(f, u=>l,g=>0.0) for f in P]
    # vars = (q̌=q̌,s=s,λ=λ,d=d,u=u,k=k,g=g)
    P,var_lens,parameters
end

function forward_once(Psys::HomotopyContinuation.System,
                        var_lens,startsols,start_parameters,target_parameters)
    tracker_options = TrackerOptions(;extended_precision=true,parameters=:default)
    result = HomotopyContinuation.solve(Psys, startsols; start_parameters, target_parameters, tracker_options, threading = false)
    # result = HomotopyContinuation.solve(Fsys, startsols)
    check_and_retrieve(result,var_lens)
end

function find_diff_system(
            P,parameters,
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
    # Psys = System(P;variable_groups,parameters=diff_parameters)
    Psys = System(P;parameters=diff_parameters)
    Psys,ide_pindx
end

function forward_once(st::Structure,
                        startsols_input,
                        start_parameters_input,
                        target_parameters_input,
                        mode=PrimalMode();F=reshape(build_G(st),:,1))

    P_all,var_lens,parameters = forward_system(st,mode;F)
    Psys, diff_parameter_points  = find_diff_system(
                P_all,parameters,
                start_parameters_input,
                target_parameters_input
            )
    q0,s0,λ0 = startsols_input
    startsols = [[q0; s0; λ0]]


    # forward_once(Psys,startsols,start_parameters,target_parameters)
    target_sol = forward_once(Psys,var_lens,startsols,diff_start_parameters,diff_target_parameters)
    check_slackness(inv.(target_sol.s),target_parameters_input.u)
    target_sol
end

function forward_sequence(Psys::HomotopyContinuation.System,
                        var_lens,startsols_input,
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
        sol = forward_once(Psys,
                    var_lens,xᵏ⁻¹,
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

function forward_sequence(st::Structure,
        startsols,
        start_parameters,
        target_parameters,
        mode=PrimalMode();
        F̌=reshape(build_Ǧ(st),:,1),
        n=1
    )
    P,var_lens,parameters = forward_system(st,mode;F̌)
    parameter_points = [
        start_parameters,
        target_parameters,
    ]
    Psys, ide_pindx = find_diff_system(
                P,parameters,
                parameter_points
            )
    # Psys = System(P;parameters)
    forward_sequence(Psys,var_lens,startsols,start_parameters,target_parameters,ide_pindx;n)
end

function forward_multi_sequence(Psys::HomotopyContinuation.System,
                        var_lens,startsols_input,
                        parameter_points, ide_pindx, mode=PrimalMode();n=1)

    # parameter_point1 = parameter_points[1]
    startsols_inputs = [startsols_input]
    [
        begin
            @debug "Forwarding the $(i)th sequence."
            seq = forward_sequence(Psys,
                                var_lens,
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

function forward_multi_sequence(st::Structure,startsols,
                        parameter_points,mode=PrimalMode();
                        F̌=reshape(build_Ǧ(st),:,1),n=1)
    P,var_lens,parameters = forward_system(st,mode;F̌)
    # q̌0,s0,λ0 = startsols
    # q̌,s,λ = variable_groups
    # Pz = map(P) do z
    #     z(q̌=>q̌0,s=>s0,λ=>λ0,parameters=>reduce(vcat,parameter_points[1]))
    # end
    # Psys = System(P;variable_groups,parameters)

    Psys, ide_pindx = find_diff_system(
                P,parameters,
                parameter_points
            )
    seqs = forward_multi_sequence(Psys,var_lens,startsols,parameter_points, ide_pindx,mode;n)
    [recover.(seq,Ref(st)) for seq in seqs]
end

function recover(state::NamedTuple,st::AbstractStructure)
    q = get_coords(st)
    (;sys_free_idx) = st.connectivity.indexed
    q[sys_free_idx] = state.q̌
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
    (;st) = bot
    q̌ = get_free_coords(st)
    isequ, λ = check_static_equilibrium_output_multipliers(bot.structure)
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

function get_start_system(bot,mode=PrimalMode();F=reshape(build_Ǧ(bot.structure),:,1))
    (;st) = bot
    start_sol,u = get_start_sol(bot)
    g = zeros(get_numbertype(bot),size(F,2))
    if typeof(mode)<:PrimalMode
        start_parameters = @eponymtuple(u,g)
    elseif typeof(mode)<:StiffMode
        k = get_cables_stiffness(st)
        start_parameters = @eponymtuple(k,u,g)
    elseif typeof(mode)<:DeformMode
        d = get_d(st)
        start_parameters = @eponymtuple(d,u,g)
    elseif typeof(mode)<:AllMode
        d = get_d(st)
        c = get_local_coords(st)
        k = get_cables_stiffness(st)
        start_parameters = @eponymtuple(d,c,k,u,g)
    else
        error("Invalid mode")
    end
    start_sol, start_parameters
end


function check_stability!(bot::Robot,Ň;
    gravity=false,
    scaling=0.01,
    scalings=nothing
)
(;st,traj) = bot
static_equilibrium,λ = check_static_equilibrium_output_multipliers(st)
@assert static_equilibrium
q̌ = get_free_coords(st)
_, Ň0, er = check_stability(bot.structure,λ,Ň;verbose=true)
resize!(traj,1)
for i in 1:length(er.values)
    push!(traj,deepcopy(traj[end]))
    traj.t[end] = er.values[i]
    δq̌i = Ň0*er.vectors[:,i]
    # @show δq̌i, er.vectors[:,i]
    if scalings isa Nothing
        si = scaling
    else
        si = scalings[i]
    end
    ratio = norm(δq̌i) / norm(q̌) 
    traj.q̌[end] .= q̌ .+ si.*δq̌i/ratio
end
bot
end

function get_poly(bot_input;
    Ň
)
bot = deepcopy(bot_input)
(;st) = bot
# (;num_of_dof,num_of_cstr,connectivity) = bot.structure
# (;cables) = st.tensiles
# (;num_of_full_coords,num_of_free_coords) = connectivity.indexed
# ncables = length(cables)
# nλ = num_of_cstr
gue = get_initial(st)
Φ = make_cstr_function(st,gue.q)
A = make_cstr_jacobian(st,gue.q)
Q̌ = make_Q̌(st,gue.q)
S = make_S(st,gue.q)
Ǩm_Ǩg = make_Ǩm_Ǩg(st,gue.q)

pv = get_polyvar(st)

pnq̌ = 1.0pv.q̌ .+ 0.0
pns = 1.0pv.s .+ 0.0
pnλ = 1.0pv.λ .+ 0.0
pnd = 1.0pv.d .+ 0.0
pnc = 1.0pv.c .+ 0.0
pnk = 1.0pv.k .+ 0.0
pnμ = 1.0pv.μ .+ 0.0
polyΦ = Φ(pnq̌,pnd,pnc)
polyA = A(pnq̌,pnc)
polyQ̌ = Q̌(pnq̌,pns,pnμ,pnk,pnc)
polyS = S(pnq̌,pns,pnc)
polyQ̌a = transpose(polyA)*pnλ
polyǨa = reduce(hcat,differentiate.(-polyQ̌a,Ref(pv.q̌))) |> transpose
polyǨm, polyǨg = Ǩm_Ǩg(pnq̌,pns,pnμ,pnk,pnc)
polyǨ = polyǨm .+ polyǨg .+ polyǨa
polyŇ = Ň(pnq̌,pnc)
poly𝒦 = transpose(polyŇ)*polyǨ*polyŇ

polyP = [
    - polyQ̌ .- transpose(polyA)*pnλ ;
    polyS;
    polyΦ;
    # poly𝒦*pnξ.-pnζ.*pnξ;
    # transpose(pnξ)*pnξ-1;
]

# Ǩ0 = RB.build_Ǩ(bot.structure,gue.λ)
# Ǩx = map(polyǨ) do z
# 		z(
# 			pv.q̌=>gue.q̌,
# 			pv.s=>gue.s,
# 			pv.λ=>gue.λ,
# 			pv.μ=>gue.μ,
# 			pv.k=>gue.k,
# 			pv.d=>gue.d,
# 			pv.c=>gue.c
# 		)
# 	end
# # @show Ǩ0
# @show Ǩ0.- Ǩx |> norm

# P0 = map(polyP) do z
# 	z(
# 		pvq̌=>q̌0,
# 		pvs=>s0,
# 		pvλ=>λ0,
# 		# pvξ=>ξ0,
# 		pvμ=>μ0,
# 		pvk=>k0,
# 	    pvd=>d0,
# 		pvc=>c0,
# 		# pv.ζ=>ζ0
# 	)
# end
# @show P0[                 1:num_of_free_coords] |> norm
# @show P0[           num_of_free_coords+1:num_of_free_coords+ncables] |> norm
# @show P0[   num_of_free_coords+ncables+1:num_of_free_coords+ncables+nλ] |> norm
# @show P0[num_of_free_coords+ncables+nλ+1:num_of_free_coords+ncables+nλ+num_of_dof]
# @show P0[end]
polyP,poly𝒦,gue,pv
end

function pinpoint(bot_input;
    Ň
)
polyP, poly𝒦, gue, pv = get_poly(bot_input;Ň)
ň = length(pv.q̌)
ns = length(pv.s)
nλ = length(pv.λ)
function make_bf()
    function inner_pp!(f,x)
        q̌x = @view x[        1:ň]
        sx = @view x[      ň+1:ň+ns]
        λx = @view x[   ň+ns+1:ň+ns+nλ]
        Px = map(polyP) do z
            z(
                pv.q̌=>q̌x,
                pv.s=>sx,
                pv.λ=>λx,
                pv.μ=>gue.μ,
                pv.k=>gue.k,
                pv.d=>gue.d,
                pv.c=>gue.c,
            )
        end

        f .= Px
    end
end
f_holder = zeros(ň+ns+nλ)
x_initial = vcat(gue.q̌,gue.s,gue.λ)
pp! = make_bf()

pp = nlsolve(pp!,x_initial,ftol=1e-10,iterations=100,method=:newton)
# @show
pp!(f_holder,pp.zero)
# @show f_holder |> norm
# @show f_holder[                 1:ň+ns+nλ] |> norm
# @show f_holder[ň+ns+nλ+1:ň+ns+nλ+num_of_dof] |> norm
# @show f_holder[end]
# @show  pp.zero[ň+ns+nλ+1:ň+ns+nλ+num_of_dof]
# @show  pp.zero[end]
q̌ = pp.zero[        1:ň]
s = pp.zero[      ň+1:ň+ns]
λ = pp.zero[   ň+ns+1:ň+ns+nλ]
ini = @eponymtuple(
        q̌,s,λ,
        isconverged=converged(pp),
        d=gue.d, c=gue.c, μ=gue.μ, k=gue.k
)
polyP, poly𝒦, ini, pv
end

function path_follow(bot_input;Ň)
polyP, poly𝒦, ini, pv = pinpoint(bot_input;Ň)
variable_groups = [pv.q̌,pv.s,pv.λ]
parameters = [pv.d;pv.c;pv.k;pv.μ]
startsols = [[ini.q̌;ini.s;ini.λ]]
start_parameters = [ini.d;ini.c;ini.k;ini.μ]
target_parameters = [ini.d;ini.c;ini.k;ini.μ.+1.0]
Psys = System(polyP;parameters)
result = HomotopyContinuation.solve(
        Psys,
        startsols;
        start_parameters,
        target_parameters,
        threading = false
)
path_results = results(result)
if length(path_results) != 1
    @show failed(result)
    error("Tracking failed.")
end
path_result1 = path_results[1]
sol = real(solution(path_result1))
q̌,s,λ = split_by_lengths(sol,length.(variable_groups))
@eponymtuple(q̌,s,λ)
end

function path_follow_critical(bot_input)
polyP, ini, pv = pinpoint_critical(bot_input)
variable_groups = [pv.q̌,pv.s,pv.λ,pv.ξ,[pv.ζ]]
parameters = [pv.d;pv.c;pv.k;pv.μ]
startsols = [[ini.q̌;ini.s;ini.λ;ini.ξ;ini.ζ]]
start_parameters = [ini.d;ini.c;ini.k;ini.μ]
target_parameters = [ini.d;ini.c;ini.k;ini.μ.+1.0]
Psys = System(polyP;parameters)
result = HomotopyContinuation.solve(
        Psys,
        startsols;
        start_parameters,
        target_parameters,
        threading = false
)
path_results = results(result)
if length(path_results) != 1
    @show failed(result)
    error("Tracking failed.")
end
path_result1 = path_results[1]
sol = real(solution(path_result1))
q̌,s,λ,ξ,ζ = RB.split_by_lengths(sol,length.(variable_groups))
@eponymtuple(q̌,s,λ,ξ,ζ)
end
