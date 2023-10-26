using LinearAlgebra
using SparseArrays
using StaticArrays
using StructArrays
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
using LaTeXStrings
using RecursiveArrayTools
using Cthulhu
using XLSX, DataFrames
using Unitful
using GeometryBasics
using TypeSortedCollections
using EponymTuples
using Match
using HomotopyContinuation
using Random
using Evolutionary
using Printf
using Revise
import Rible as RB
cd("examples/manipulator")
includet("man_define.jl")
includet("../analysis.jl")
includet("../vis.jl")
includet("man_plotting.jl")
k = 453.09; c = 100.0;
restlen = 0.174-3.9061/k
man_inv = man_ndof(8,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=false); man_inv_plot = deepcopy(man_inv)

function sup!(ax)
    cir = GeometryBasics.Circle(Point2(0.061260,-0.21),0.1)
    mesh!(ax, cir)
end

plot_traj!(man_inv;sup!)

function dist!(st,state)
    target = GeometryBasics.Circle(Point2(0.061260,-0.21),0.1)
    T = RB.get_numbertype(st)
    RB.update_points!(st,state.c)
    RB.update_rigids!(st,state.q)
    d = Ref(zero(T))
    foreach(st.rigidbodies) do rb
        rbid = rb.prop.id
        if rbid in collect(1:2:11)
            d[] += abs(norm(rb.state.rps[7] .- origin(target)) - radius(target))
        end
    end
    d[]
end

function Evolutionary.trace!(record::Dict{String,Any}, objfun, state, population, method::ES, options)
    s = strategy(state)
    record["σ"] = s.σ
    # record["population"] = population[1]
end

function get_actuate_seqs(bot,g0=[0.0])
    (;st) = bot
    (;nrigids) = st
    (;mem2num,num2sys) = st.connectivity.numbered
    tg_end = deepcopy(st)
    nsegs = nrigids-1
    startsols,parameters0 = RB.get_start_system(bot,RB.AllMode())

    c0 = parameters0.c
    u0 = parameters0.u
    d0 = parameters0.d
    k0 = parameters0.k
    g1 = g0

    indx_d = reduce(
        vcat,
        [
            [3(i-1)+1,3(i-1)+2] for i = 1:nrigids
        ]
    )
    # indx_d = collect(1:length(d0))
    indx_c = reduce(
        vcat,
        [
            reduce(vcat,num2sys[mem2num[i][4:7]])
            for i = 1:nrigids
        ]
    )
    # indx_c = collect(1:length(c0))
    indx_u = collect(1:length(u0))
    d1 = copy(d0)
    d1[indx_d].+=1
    c1 = copy(c0)
    c1[indx_c].+=1
    u1 = copy(u0)
    u1[indx_u].+=1

    dummy_parameters1 = (
        d=d1,
        c=c1,
        k=k0,
        u=u1,
        g=g0
    )
    dummy_parameter_points = [parameters0,dummy_parameters1]

    P,variable_groups,parameters = RB.forward_system(st,RB.AllMode();)
    Psys, ide_pindx = RB.find_diff_system(P,variable_groups,parameters,dummy_parameter_points)

    indr1p7 = num2sys[mem2num[1][7]]
    indr2p5 = num2sys[mem2num[2][5]]
    @show indr1p7
    @show c0[indr1p7]

    d1 .= d0
    c1 .= c0
    u1 .= u0
    a1 = 0.01
    u1[collect(1:2:2nsegs)] .+= a1
    u1[collect(2:2:2nsegs)] .-= a1


    xlens = length.([indx_d,indx_c])
    function f_core(x)
        d1[indx_d], c1[indx_c] = RB.split_by_lengths(x,xlens)
        parameters1 = (d=d1, c=c1, k=k0, u=u1, g=g1)
        seq_raw = RB.forward_sequence(Psys,
                            startsols,
                            parameters0,
                            parameters1,
                            ide_pindx;n=1)
        seq_rc = RB.recover.(seq_raw,Ref(st))
    end
    function f(x)
        seq_rc = f_core(x)
        sumdis = dist!(tg_end,seq_rc[end])
        # @show
        # seqs_rc
    end
    x0 = vcat(d0[indx_d],c0[indx_c])
    @show d0[16]
    csts = BoxConstraints(
        x0.-0.05,
        x0.+0.05
    )
    # @show csts.bounds.bx[[31,32]]
    seq_ref0 = f_core(x0)
    @show dist!(tg_end,seq_ref0[end])
    Random.seed!(2);
    N = length(x0)
    es = ES(
        initStrategy=IsotropicStrategy(0.0002, 1.0/sqrt(2N), 1.0/sqrt(2*sqrt(N))),
        recombination=average,
        srecombination=average,
        mutation=gaussian,
        smutation=gaussian,
        μ=1,
        ρ=1,
        λ=5,
        selection=:plus
    )
    opts = Evolutionary.Options(
        show_trace = true,
        successive_f_tol = 10,
        # parallelization = :thread,
    )
    rst = Evolutionary.optimize(f, csts, x0, es, opts)
    f_core(rst.minimizer),seq_ref0
end
seq1, seq1_ref = get_actuate_seqs(man_inv)
@descend_code_warntype get_actuate_seqs(man_inv)
seq1[end].d - seq1[begin].d .|> abs |> findmax
seq1[end].d[16]
seq1[begin].d[16]
seq1[end].c - seq1[begin].c .|> abs |> maximum
seqs_diff = StructArray([seq1_ref[end],seq1[end]])
RB.apply!(man_inv_plot,seqs_diff)
dist!(man_inv_plot.st,seqs_diff[end])
plot_traj!(man_inv_plot;sup!)



fig = Figure()
ax = Axis(fig[1,1])
ax.limits



rosenbrock(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
Random.seed!(2);
N = 60
es = ES(
    initStrategy=IsotropicStrategy(N),
    recombination=average,
    srecombination=average,
    mutation=gaussian,
    smutation=gaussian,
    μ=1,
    ρ=1,
    λ=100,
    selection=:plus
)
opts = Evolutionary.Options(
    show_trace = true,
    parallelization = :thread,
)
Evolutionary.optimize(rosenbrock, randn(N), es, opts)

csts = BoxConstraints(
    -2ones(N),
     ones(N)
)
csts.bounds.ineqx
csts.bounds.σx

randn(N)

function plot_angles(bot,seq,seq_g)
    (;traj, st) = bot

    df = DataFrame(XLSX.readtable("弹簧2验证水平0924.xlsx", "Sheet2";infer_eltypes=true))
    Δu_read = 100df.a
    θ_read = df.theta
    θ_g_read = df.theta_G

    # seq 1
    u = [u[1] for u in seq.u]
    RB.set_new_initial!(bot,seq.q[1])
    append!(traj.q, seq.q[2:end])
    append!(traj.t,     u[2:end])
    θ = [
        begin
            RB.update_rigids!(st,q)
            get_angles(bot)[1]
        end
        for q in seq.q
    ]

    # seq 2
    u_g = [u[1] for u in seq_g.u]
    RB.set_new_initial!(bot,seq_g.q[1])
    append!(traj.q, seq_g.q[2:end])
    append!(traj.t,     u_g[2:end])

    θ_g = [
        begin
            RB.update_rigids!(st,q)
            get_angles(bot)[1]
        end
        for q in seq_g.q
    ]


    fig_width = 245.71811u"pt"
    fig_height = 110.89u"pt"
    font_size = 10u"pt"
    dpi=300/(1u"inch")
    font_pixel = to_resolution(dpi,font_size)
    resolution = to_resolution.(dpi,(fig_width,fig_height))
    # @show resolution
    fig = Figure(;resolution)
    ax1 = Axis(fig[1,1])
    ax2 = Axis(fig[1,2])
    scatter!(ax1,Δu_read,θ,markersize=20)
    scatter!(ax1,Δu_read,θ_read,markersize=20,marker='▲')
    scatter!(ax2,Δu_read,θ_g,markersize=20,label="Sim.")
    scatter!(ax2,Δu_read,θ_g_read,markersize=20,marker='▲',label="Exp.")
    axislegend(ax2,position = :rt)
    xlims!(ax1,-1.2,1.2)
    xlims!(ax2,-1.2,1.2)
    ylims!(ax1,-36,36)
    ylims!(ax2,-36,36)
    ax2.yticksvisible = false
    ax2.yticksvisible = false
    ax2.yticklabelsvisible = false
    ax1.xlabel = "Δu (cm)"
    ax1.ylabel = "θ (°)"
    ax2.xlabel = "Δu (cm)"
    fig
end

GM.activate!(); plot_angles(man_inv,seq,seq_g)
GM.activate!(); plot_angles(man_inv,seq,seq_g,"fig_angles")

plot_traj!(man_inv)


k = 1203.09; c = 100.0;
restlen = 0.174-3.9061/k
man_inv = man_ndof(4,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=true)
plotstructure(man_inv)

function get_actuate_seqs(bot,g=0.0)
    # plotstructure(man_inv)
    Y = RB.build_Y(bot)
    # q,_ = RB.get_q(bot.st)
    # λ,u,a = RB.inverse(bot,deepcopy(bot),Y)
    start_sol,parameters0 = RB.get_start_system(bot,Y)
    # start_sol
    # parameters0.u
    # parameters0.g
    # RB.actuate!(bot,a)
    # RB.check_static_equilibrium(bot.st,q,λ)
    # 0.182-5.2024/292.14
    u1 = copy(parameters0.u)
    u2 = copy(parameters0.u)
    a2 = 0.01
    u2[[1,4,5,8]] .+= a2
    u2[[2,3,6,7]] .-= a2
    u3 = copy(parameters0.u)
    u4 = copy(parameters0.u)
    a4 = -0.01
    u4[[1,4,5,8]] .+= a4
    u4[[2,3,6,7]] .-= a4
    # u1,u2

    g1 = 1.0g
    g2 = 1.0g
    g3 = 1.0g
    g4 = 1.0g

    # result = RB.forward_once(bot.st,(q,s,λ),(u,0.0),(u1,g1))

    parameters1 = (u=u1, g=[g1])
    parameters2 = (u=u2, g=[g2])
    parameters3 = (u=u3, g=[g3])
    parameters4 = (u=u4, g=[g4])
    parameters_seq = [parameters0,parameters1,parameters2,parameters3,parameters4]

    man_sys,_ = RB.forward_system(bot.st)

    # seq = RB.forward_sequence(man_sys,start_sol,start_parameters,target_parameters)
    seqs = RB.forward_multi_sequence(man_sys,start_sol,parameters_seq;n=10)

    # RB.set_new_initial!(bot,result.q)
    # plotstructure(bot)
    # get_angles(bot)
    # plotstructure(bot,seqs[1])
    # plotstructure(bot,seqs[2])
    seqs
end

k = 403.09; c = 100.0;
restlen = 0.174-3.9061/k
man_inv = man_ndof(4,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=true)
seqs1 = get_actuate_seqs(man_inv,0.0)

man_inv = man_ndof(4,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=true)
seqs2 = get_actuate_seqs(man_inv,1.0)

k = 1203.09; c = 100.0;
man_inv = man_ndof(4,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=true)
seqs3 = get_actuate_seqs(man_inv,1.0)

plotstructure(man_inv,seqs[4])

function plot_sweeping!(ax,bot,seq1,seq2)
    @unpack st, traj = bot

    for i in 1:length(seq1.q)
        RB.distribute_q_to_rbs!(st,seq1.q[i])
        lower_bars = get_lower_bars(st)
        for bar in lower_bars
            linesegments!(ax, bar, color = :grey, linewidth = 3)
        end
    end

    RB.set_new_initial!(bot,seq1.q[begin])
    append!(traj.qs,seq1.q[begin+1:end])
    endpoint_line = get_trajectory(bot,5,2)
    # @show endpoint_line
    lines!(ax, endpoint_line, color = :grey, linewidth = 3)

    for i in 1:length(seq2.q)
        RB.distribute_q_to_rbs!(st,seq2.q[i])
        upper_bars = get_upper_bars(st)
        for bar in upper_bars
            linesegments!(ax, bar, color = :grey, linewidth = 3)
        end
    end

    RB.set_new_initial!(bot,seq2.q[begin])
    append!(traj.qs,seq2.q[begin+1:end])
    endpoint_line = get_trajectory(bot,5,2)
    # @show endpoint_line
    lines!(ax, endpoint_line, color = :grey, linewidth = 3)

    bars,cables = plotstructure!(ax,st)
    update_scene!(st,bars,cables,seq1.q[begin])
    ax.aspect = DataAspect()
end

function plot_sweeping_compare(bot,seqs1,seqs2,seqs3)
    fig_width = 1.5*245.71811u"pt"
    fig_height = 150.89u"pt"
    font_size = 10u"pt"
    dpi=300/(1u"inch")
    font_pixel = to_resolution(dpi,font_size)
    resolution = to_resolution.(dpi,(fig_width,fig_height))
    fig = Figure(;resolution)
    ax1 = fig[1, 1] = Axis(fig)
    ax2 = fig[1, 2] = Axis(fig)
    ax3 = fig[1, 3] = Axis(fig)
    plot_sweeping!(ax1,bot,seqs1[2],seqs1[4])
    plot_sweeping!(ax2,bot,seqs2[2],seqs2[4])
    plot_sweeping!(ax3,bot,seqs3[2],seqs3[4])
    for ax in [ax1,ax2,ax3]
        xlims!(ax,-0.1,0.7)
        ylims!(ax,-0.5,0.4)
        ax.xticks = [0.0,0.2,0.4,0.6]
        ax.yticks = [-0.4,-0.2,0.0,0.2]

    end
    ax2.yticksvisible = false
    ax2.yticksvisible = false
    ax2.yticklabelsvisible = false
    ax3.yticksvisible = false
    ax3.yticksvisible = false
    ax3.yticklabelsvisible = false
    ax1.ylabel = "y (m)"
    ax1.xlabel = "x (m)"
    ax2.xlabel = "x (m)"
    ax3.xlabel = "x (m)"
    for (ilabel,label) in enumerate(["(a)", "(b)", "(c)"])
        Label(fig.layout[1, ilabel, TopLeft()], label,
            textsize = font_pixel,
            # font = noto_sans_bold,
            padding = (0, 0, 0, 0),
            halign = :right)
    end
    fig
end

fig_sweeping_compare = plot_sweeping_compare(man_inv,seqs1,seqs2,seqs3)
GLMakie.activate!()
CairoMakie.activate!(type="svg")
CairoMakie.save(raw"D:\OneDrive - 中山大学\Papers\ForwardStatics\ieeeconf\images\fig_sweeping_compare.svg", fig_sweeping_compare)

plot_sweeping(man_inv,seqs[2],seqs[4])


k = 453.09; c = 100.0;
restlen = 0.174-3.9061/k
man_inv = man_ndof(8,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=true); man_inv_plot = deepcopy(man_inv)

function get_actuate_video(bot_input;g0=[0.0])
    bot = deepcopy(bot_input)
    (;st) = bot
    (;nrigids) = st
    (;mem2num,num2sys) = st.connectivity.numbered
    tg_end = deepcopy(st)
    nsegs = nrigids-1
    startsols,parameters0 = RB.get_start_system(bot,RB.DeformMode())

    u0 = parameters0.u
    d0 = parameters0.d
    g1 = g0

    indx_d = reduce(
        vcat,
        [
            [3(i-1)+1,3(i-1)+2] for i = 1:nrigids
        ]
    )
    indx_u = collect(1:length(u0))
    d1 = copy(d0)
    d1[indx_d].+=0.0001
    u1 = copy(u0)
    u1[indx_u].+=0.0001

    dummy_parameters1 = (
        d=d1,
        u=u1,
        g=g0
    )
    dummy_parameter_points = [parameters0,dummy_parameters1]

    P,variable_groups,parameters = RB.forward_system(st,RB.DeformMode();)
    Psys, ide_pindx = RB.find_diff_system(P,variable_groups,parameters,dummy_parameter_points)

    # precompile
    RB.forward_once(Psys,
                        [[startsols.q̌; startsols.s; startsols.λ]],
                        reduce(vcat,(d=d0[indx_d], u=u0, )),
                        reduce(vcat,(d=d1[indx_d], u=u1, ))
        )

    set_theme!(theme_try;fontsize=6 |> pt2px)
    fig = Figure(resolution=(1920,1080))
    ax = Axis(fig[1,1])
    ax.aspect = DataAspect()

    bars,cables = plotstructure!(ax,st)

    sg_d = SliderGrid(fig[1,2],
        [
            (
                label = latexstring("d_{$i}"),
                range = LinRange(d0[i]-0.01,d0[i]+0.01,11),
                startvalue = d0[i]
            )
            for i in indx_d
        ]...
    )

    sg_a = SliderGrid(fig[1,3],
        [
            (
                label = latexstring("\\Delta u_{$i}"),
                range = LinRange(-0.01,0.01,11),
                startvalue = 0.0
            )
            for i in collect(1:2:2nsegs)
        ]...
    )
    colsize!(fig.layout,1,Relative(0.6))

    sg_values = vcat(
        [s.value for s in sg_d.sliders],
        [s.value for s in sg_a.sliders]
    )
    sg_lens = length.([sg_d.sliders,sg_a.sliders])
    onany(sg_values...) do vals_tuple...
        vals = [v for v in vals_tuple]
        dv, av = RB.split_by_lengths(vals,sg_lens)
        d1 .= d0
        d1[indx_d] = dv
        u1 .= u0
        # a1 = 0.01
        u1[collect(1:2:2nsegs)] .+= av
        u1[collect(2:2:2nsegs)] .-= av

        rst = RB.forward_once(Psys,
                            [[startsols.q̌; startsols.s; startsols.λ]],
                            reduce(vcat,(d=d0[indx_d], u=u0, )),
                            reduce(vcat,(d=d1[indx_d], u=u1, ))
            )
        rst_rc = RB.recover(rst,st)
        update_scene!(st,bars,cables,rst_rc.q)
    end
    xlims!(ax,-0.1,1.2)
    ylims!(ax,-1.1,0.2)
    fig
end
get_actuate_video(man_inv)
seq1 = get_actuate_video(man_inv)
RB.apply!(man_inv_plot,StructArray([seq1]))
plot_traj!(man_inv_plot;)

plotstructure(man_inv_plot)
