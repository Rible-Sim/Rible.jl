using LinearAlgebra
using SparseArrays
using StaticArrays
using Makie
import GLMakie as GM
import CairoMakie as CM
GM.activate!()
# import PyPlot as plt
# plt.pygui(true)
using LaTeXStrings
using RecursiveArrayTools
using Cthulhu
using XLSX, DataFrames
using Unitful
using BenchmarkTools
# using NLsolve
# using DifferentialEquations
using HomotopyContinuation
using Printf
using Statistics
using GeometryBasics
# using FunctionWrappers
using EponymTuples
using Match
using TypeSortedCollections
using Revise
import TensegrityRobots as TR
cd("examples/DualTriangles")
includet("man_define.jl")
# includet("man_plotting.jl")
includet("../analysis.jl")
includet("../vis.jl")

figdir = raw"C:\Users\luo22\OneDrive\Papers\Ph.D.Thesis\ff"

k = 1275; c = 100.0;
restlen = 0.163-2.7246/318.16
man_inv = man_ndof(2,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=false); man_inv_plot = deepcopy(man_inv)
man_inv |> typeof |> isconcretetype
plot_traj!(man_inv)
man_inv.tg.rigidbodies.data[1]
@code_warntype man_ndof(1,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=false)
@descend_code_warntype man_ndof(1,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=false)
rb1 = man_inv.tg.rigidbodies.data[1][1]
TR.update_points!(man_inv.tg)


(;mem2num,num2ID,num2sys) = man_inv.tg.connectivity.numbered
TR.get_c(man_inv)[indr1p2]

# !!!! Check inverse statics first !!!!
function get_testall_seq(bot;g=0.0,amax=0.01,amin=-0.01,n=20)
    (;mem2num,num2ID,num2sys) = bot.tg.connectivity.numbered
    indr1p2 = num2sys[mem2num[1][2]]
    indr2p5 = num2sys[mem2num[2][5]]
    start_sol,parameters0 = TR.get_start_system(bot,TR.AllMode())
    c0 = parameters0.c
    u0 = parameters0.u
    d0 = parameters0.d
    k0 = parameters0.k
    @show indr1p2,indr2p5
    @show c0[indr1p2]
    @show c0[indr2p5]
    u1 = copy(u0)
    a1 = amax
    u1[[1,3]] .+= a1
    u1[[2,4]] .-= a1
    u2 = copy(u0)
    a2 = amin
    u2[[1,3]] .+= a2
    u2[[2,4]] .-= a2

    g1 = g
    g2 = g

    c1 = copy(c0)
    c2 = copy(c0)
    c2[indr1p2] = [0.11,0.02]
    c2[indr2p5] = [0.0 ,0.09]
    parameters1 = (d=d0, c=c1, k=k0, u=u0, g=[g1], )
    parameters2 = (d=d0, c=c2, k=k0, u=u0, g=[g2], )
    parameters_seq = [parameters0,parameters1,parameters2]
    seqs = TR.forward_multi_sequence(bot.tg,start_sol,parameters_seq,TR.AllMode();n)

    seqs[2]
end
seq = get_testall_seq(man_inv)
VectorOfArray(seq.c)[[25,26],:]
TR.apply!(man_inv_plot,seq)
plot_traj!(man_inv_plot)
# !!!! Check inverse statics first !!!!
function get_actuate_seq(bot;g=0.0,amax=0.01,amin=-0.01,n=20)
    # plot_traj!(man_inv)
    # q,_ = TR.get_q(bot.tg)
    # λ,u,a = TR.inverse(bot,deepcopy(bot),Y)
    # @show Y
    start_sol,parameters0 = TR.get_start_system(bot)
    @show parameters0.u
    # @show TR.get_strings_len(bot)
    # @show start_sol
    # @show parameters0
    # start_sol
    # parameters0.u
    # parameters0.g
    # TR.actuate!(bot,a)
    # TR.check_static_equilibrium(bot.tg,q,λ)
    # 0.182-5.2024/292.14
    u1 = copy(parameters0.u)
    a1 = amax
    u1[[1,3]] .+= a1
    u1[[2,4]] .-= a1
    u2 = copy(parameters0.u)
    a2 = amin
    u2[[1,3]] .+= a2
    u2[[2,4]] .-= a2
    # u1,u2

    g1 = g
    g2 = g

    # result = TR.forward_once(bot.tg,(q,s,λ),(u,0.0),(u1,g1))

    parameters1 = (u=u1, g=[g1])
    parameters2 = (u=u2, g=[g2])
    parameters_seq = [parameters0,parameters1,parameters2]

    # seq = TR.forward_sequence(man_sys,start_sol,start_parameters,target_parameters)
    seqs = TR.forward_multi_sequence(bot.tg,start_sol,parameters_seq;n)

    # u = [u[1] for u in seqs[2].u]
    # TR.set_new_initial!(bot,seqs[2].q[1])
    # append!(bot.traj.q, seqs[2].q[2:end])
    # append!(bot.traj.t,         u[2:end])
    # plot_traj!(bot)
    # get_angles(bot)
    # plot_traj!(bot,seqs[1])
    # plot_traj!(bot,seqs[2])
    seqs[2]
end
seq = get_actuate_seq(man_inv)

seq_g = get_actuate_seq(man_inv;g=1.0,amax=0.00,amin=-0.015,n=15)
plot_traj!(man_inv,seq_g)

function cable_angle(bot_input,seq)
    bot = deepcopy(bot_input)
    (;tg) = bot
    θs = Vector{Float64}()
    for q in seq.q
        TR.update_rigids!(tg,q)
        TR.update_tensiles!(tg)
        sdir = tg.tensiles.cables[2].state.direction
        θ = atan(sdir[2],sdir[1])
        push!(θs,θ)
    end
    # θs
    (maximum(θs)-minimum(θs))*4.75
end
θs = cable_angle(man_inv,seq)



man_inv = man_ndof(1,[1.0,0.0];θ=deg2rad(-2.2571854917431766),k,c,restlen,isvirtual=true)
q, = TR.get_q(man_inv.tg)

df = DataFrame(XLSX.readtable("三次实验平均1-9.xlsx", "Sheet4";infer_eltypes=true)...)
theta1_raw = sort(df.theta1_raw)
theta2_raw = sort(df.theta2_raw)
theta1_raw_diff = theta1_raw[begin+1:end] .- theta1_raw[begin:end-1]
theta2_raw_diff = theta2_raw[begin+1:end] .- theta2_raw[begin:end-1]
theta_raw = sort(vcat(theta1_raw_diff,theta2_raw_diff))
scatter(theta_raw)

function plot_angles(bot,seq,seq_g,figname=nothing)
    (;traj, tg) = bot

    df   = DataFrame(XLSX.readtable("三次实验平均1-9.xlsx", "Sheet2";infer_eltypes=true))
    df_g = DataFrame(XLSX.readtable("三次实验平均1-9.xlsx", "Sheet3";infer_eltypes=true))
    Δu_read = 1000df.a
    θ1_read = df.theta1
    θ2_read = df.theta2
    Δu_g_read = 1000df_g.a_g
    θ1_g_read = df_g.theta1_g
    θ2_g_read = df_g.theta2_g

    # seq 1
    u = [u[1] for u in seq.u]
    TR.set_new_initial!(bot,seq.q[1])
    append!(traj.q, seq.q[2:end])
    append!(traj.t,     u[2:end])

    θ1 = [
        begin
            TR.update_rigids!(tg,q)
            get_angles(bot)[1]
        end
        for q in seq.q
    ]

    θ2 = [
        begin
            TR.update_rigids!(tg,q)
            get_angles(bot)[2]
        end
        for q in seq.q
    ]

    # seq 2
    u_g = [u[1] for u in seq_g.u]
    TR.set_new_initial!(bot,seq_g.q[1])
    append!(traj.q, seq_g.q[2:end])
    append!(traj.t,     u_g[2:end])

    θ1_g = [
        begin
            TR.update_rigids!(tg,q)
            get_angles(bot)[1]
        end
        for q in seq_g.q
    ]

    θ2_g = [
        begin
            TR.update_rigids!(tg,q)
            get_angles(bot)[2]
        end
        for q in seq_g.q
    ]


    # fig_width = 245.71811u"pt"
    # fig_height = 155.89u"pt"
    # font_size = 10u"pt"
    # dpi=300/(1u"inch")
    # font_pixel = to_resolution(dpi,font_size)
    # resolution = to_resolution.(dpi,(fig_width,fig_height))
    # @show resolution
    with_theme(theme_pub;
            resolution = (0.7tw,0.4tw),
            Scatter = (
                markersize = fontsize,
            )
        ) do
        fig = Figure(;)
        ax1 = Axis(fig[1,1], ylabel = L"\theta_1(\degree)")
        ax2 = Axis(fig[2,1], ylabel = L"\theta_2(\degree)")
        ax3 = Axis(fig[1,2])
        ax4 = Axis(fig[2,2])
        @show maximum(abs.(θ1.-θ1_read))
        @show maximum(abs.(θ2.-θ2_read))
        @show maximum(abs.(θ1_g.-θ1_g_read))
        @show maximum(abs.(θ2_g.-θ2_g_read))
        scatter!(ax1,Δu_read,θ1,marker=:circle)
        scatter!(ax1,Δu_read,θ1_read,marker=:diamond)
        scatter!(ax2,Δu_read,θ2,marker=:circle)
        scatter!(ax2,Δu_read,θ2_read,marker=:diamond)
        scatter!(ax3,Δu_g_read,θ1_g,marker=:circle)
        scatter!(ax3,Δu_g_read,θ1_g_read,marker=:diamond)
        scatter!(ax4,Δu_g_read,θ2_g,label="Sim.",marker=:circle)
        scatter!(ax4,Δu_g_read,θ2_g_read,marker=:diamond,label="Exp.")
        axislegend(ax4,position = :lb)
        xlims!(ax1,-11,11); xlims!(ax2,-11,11)
        xlims!(ax3,-17,02); xlims!(ax4,-17,02)
        ylims!(ax1,-38,36); ylims!(ax2,-38,36)
        ylims!(ax3,-38,36); ylims!(ax4,-38,36)
        ax1.xticklabelsvisible = false
        ax3.xticklabelsvisible = false
        ax3.yticklabelsvisible = false
        ax4.yticklabelsvisible = false
        # ax3.yticksvisible = false
        # ax3.yticksvisible = false
        # ax2.xlabel = "Δμ (cm)"
        for (ilabel,label) in enumerate(["(a)", "(b)"])
            Label(fig.layout[ilabel, 1, TopLeft()], label,)
        end
        for (ilabel,label) in enumerate(["(c)", "(d)"])
            Label(fig.layout[ilabel, 2, TopLeft()], label,)
        end
        # Label(fig[1:2, 0], "θ (°)", rotation = pi/2)
        Label(fig[3, 1:2], L"\Delta\mu (\mathrm{mm})", halign = :center)

        savefig(fig,figname)
    end
end

GM.activate!(); plot_angles(man_inv,seq,seq_g)
CM.activate!(); plot_angles(man_inv,seq,seq_g,"fig_angles")

k = 558.16; c = 100.0;
restlen = 0.163-2.7246/318.16
man_inv = man_ndof(3,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=true)
GM.activate!(); plot_traj!(man_inv)

function get_actuate_seqs(bot,g=0.0)

    start_sol,parameters0 = TR.get_start_system(bot)
    # start_sol.λ ./= 1e12
    display(start_sol.λ)
    u1 = copy(parameters0.u)
    u2 = copy(parameters0.u)
    a2 = 0.01
    u2[[1,3,5]] .+= a2
    u2[[2,4,6]] .-= a2
    u3 = copy(parameters0.u)
    u4 = copy(parameters0.u)
    a4 = -0.01
    u4[[1,3,5]] .+= a4
    u4[[2,4,6]] .-= a4

    g1 = 1.0g
    g2 = 1.0g
    g3 = 1.0g
    g4 = 1.0g

    parameters1 = (u=u1, g=[g1])
    parameters2 = (u=u2, g=[g2])
    parameters3 = (u=u3, g=[g3])
    parameters4 = (u=u4, g=[g4])
    parameters_seq = [parameters0,parameters1,parameters2,parameters3,parameters4]

    seqs = TR.forward_multi_sequence(bot.tg,start_sol,parameters_seq;n=10)

    seqs
end

k = 558.16; c = 100.0;
restlen = 0.163-2.7246/318.16
man_inv = man_ndof(3,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=true); man_video = deepcopy(man_inv)
seqs1 = get_actuate_seqs(man_inv,0.0)
seqs1[2].nstep[2:end]
seqs1[4].nstep[2:end]
plot_traj!(man_inv,seqs1[2])
# TR.get_strings_stiffness(man_inv)
man_inv = man_ndof(3,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=true)
seqs2 = get_actuate_seqs(man_inv,1.0)
# plot_traj!(man_inv,seqs2[1])

k = 1200.0; c = 100.0;
man_inv = man_ndof(3,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=true)
seqs3 = get_actuate_seqs(man_inv,1.0)
seqs3[2].nstep[2:end]
seqs3[4].nstep[2:end]
nsteps = reduce(vcat,[seqs1[2].nstep[2:end],seqs2[2].nstep[2:end],seqs3[2].nstep[2:end],
                      seqs1[4].nstep[2:end],seqs2[4].nstep[2:end],seqs3[4].nstep[2:end]])
stats_nstep = nsteps |> (x) -> [minimum(x), maximum(x), median(x), mean(x)]
stats_nstep
seqs3[4].q
TR.apply!(man_inv,seqs3[4])
plot_traj!(man_inv)

function plot_sweeping!(ax,bot,seq1,seq2,colormap)
    (;tg, traj) = bot
    imid = cld(length(colormap),2)
    for i in 1:length(seq1.q)
        TR.update_rigids!(tg,seq1.q[i])
        lower_bars = get_lower_bars(tg)
        for bar in lower_bars
            linesegments!(ax, bar, color = colormap.colors[imid-1+i], linewidth = 3)
        end
    end

    TR.set_new_initial!(bot,seq1.q[begin])
    append!(traj.qs,seq1.q[begin+1:end])
    endpoint_line = get_trajectory(bot,4,2)
    endpoint_linesegments = VectorOfArray(reduce(vcat,[[ep,ep] for ep in endpoint_line.u])[begin+1:end-1])
    endpoint_colors =  reduce(vcat,[[colori, colori] for colori in colormap.colors[imid:end]])[begin+1:end-1]
    linesegments!(ax, endpoint_linesegments, color = endpoint_colors, linewidth = 3)

    for i in 1:length(seq2.q)
        TR.update_rigids!(tg,seq2.q[i])
        upper_bars = get_upper_bars(tg)
        for bar in upper_bars
            linesegments!(ax, bar, color = colormap.colors[imid+1-i], linewidth = 3)
        end
    end

    TR.set_new_initial!(bot,seq2.q[begin])
    append!(traj.qs,seq2.q[begin+1:end])
    endpoint_line = get_trajectory(bot,4,2)
    endpoint_linesegments = VectorOfArray(reduce(vcat,[[ep,ep] for ep in endpoint_line.u])[begin+1:end-1])
    endpoint_colors =  reduce(vcat,[[colori, colori] for colori in colormap.colors[imid:-1:1]])[begin+1:end-1]
    linesegments!(ax, endpoint_linesegments, color = endpoint_colors, linewidth = 3)

    # bars,strings = plot_traj!(ax,tg)
    # update_scene!(tg,bars,strings,seq1.q[begin])
    ax.aspect = DataAspect()
end

function plot_sweeping_compare(bot,seqs1,seqs2,seqs3)
    fig_width = 1.55*245.71811u"pt"
    fig_height = 150.89u"pt"
    font_size = 10u"pt"
    dpi=300/(1u"inch")
    font_pixel = to_resolution(dpi,font_size)
    resolution = to_resolution.(dpi,(fig_width,fig_height))
    fig = Figure(;resolution)
    ax1 = fig[1, 1] = Axis(fig)
    ax2 = fig[1, 2] = Axis(fig)
    ax3 = fig[1, 3] = Axis(fig)
    colormap = cgrad(:RdYlGn, 21, categorical = true)
    cb = Colorbar(fig[1, 4]; limits = (-10, 10),
                colormap,
                size = 30,
                label = "Δμ(mm)")
    edges = range(-10, 10, length = 22)
    centers = (edges[1:5:21] .+ edges[2:5:22]) .* 0.5
    cb.ticks = (centers, string.(-10:5:10))
    plot_sweeping!(ax1,bot,seqs1[2],seqs1[4],colormap)
    plot_sweeping!(ax2,bot,seqs2[2],seqs2[4],colormap)
    plot_sweeping!(ax3,bot,seqs3[2],seqs3[4],colormap)
    for ax in [ax1,ax2,ax3]
        xlims!(ax,-0.1,0.7)
        ylims!(ax,-0.5,0.5)
        ax.xticks = [0.0,0.2,0.4,0.6]
        ax.yticks = [-0.4,-0.2,0.0,0.2,0.4]

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
CairoMakie.activate!(type="pdf")
CairoMakie.save(texroot*raw"\OneDrive - 中山大学\Papers\ForwardStatics\ieeeconf\images\fig_sweeping_compare.pdf", fig_sweeping_compare)

plot_sweeping(man_inv,seqs[2],seqs[4])

function get_actuate_video(bot_input)
    bot = deepcopy(bot_input)
    fig = Figure(resolution = (1280, 720))
    ax1 = fig[1, 1] = Axis(fig)
    ax1.aspect = DataAspect()
    ax1.tellheight = false
    bars,strings = plot_traj!!(ax1,bot.tg)
    xlims!(ax1,-0.1,0.7)
    ylims!(ax1,-0.5,0.4)
    ls_g = labelslider!(fig, "Gravity:", 0:0.1:1;
                format = x -> "$(x)×9.8m/s²", valuekw = Dict(:font=>"CMU Serif"))
    set_startindex(ls_g.slider,0.0)

    Y = TR.build_Y(bot)
    start_sol,parameters0 = TR.get_start_system(bot,Y,TR.StiffMode())
    s = start_sol.s
    ℓ = 1 ./s
    k0,u0,g0 = parameters0
    # na = size(Y,2)
    # alabels = ["a$i:" for i in 1:na]
    # amin, amax = -0.01, 0.01
    # aranges = Ref(amin:((amax-amin)/10):amax)
    # aformats = Ref(x -> begin @sprintf "%+.5fm" x end)
    # ls_a = labelslidergrid!(fig,alabels,aranges;aformats)
    # set_startindex.(ls_a.sliders,zeros(na))
    # fig[:,2] = ls_a.layout

    Δumin, Δumax = -0.01, 0.01
    kmin, kmax = 300.0, 1200.0
    Δustart = 0.0
    kstart = kmin
    sliderlabels = [L"Δu:",L"k:"]
    sliderranges = [Δumin:0.001:Δumax,kmin:100.0:kmax]
    sliderformats = [x -> begin
                        if x >= 0
                            @sprintf "+%.3fm" x
                        else
                            @sprintf "–%.3fm" -x
                        end
                    end
                    , x -> begin @sprintf "%6.1fN/m" x end]
    ls_slider = labelslidergrid!(fig,sliderlabels,sliderranges;formats=sliderformats,valuekw = Dict(:font=>"CMU Serif"))
    set_startindex.(ls_slider.sliders,[Δustart,kstart])


    # nk = length(k0)
    # klabels = ["k$i:" for i in 1:nk]
    # kmin, kmax = k0[1], k0[1]*10
    # kranges = Ref(kmin:((kmax-kmin)/10):kmax)
    # kformats = Ref(x -> begin @sprintf "%.2fN/m" x end)
    # @show k0
    # ls_k = labelslidergrid!(fig,klabels,kranges;kformats)
    # set_startindex.(ls_k.sliders,k0)
    fig[1,2] = ls_slider.layout
    fig[2,1] = ls_g.layout

    targets_ob = vcat(
        [sl.value for sl in ls_slider.sliders],
        [ls_g.slider.value]
    )
    onany(targets_ob...) do targets...
        Δu = targets[1]
        a1 = [Δu,-Δu,Δu]
        u1 = copy(u0)
        u1 += Y*a1
        k  = targets[2]
        k1 = copy(k0)
        k1 .= k
        g1 = targets[end]
        # @show a1
        parameters1 = (k=k1, u=u1, g=[g1])
        result =  TR.forward_once(bot.tg,start_sol,parameters0,parameters1,TR.StiffMode())
        update_scene!(bot.tg,bars,strings,result.q)
    end
    Label(fig[0, :], text = "Forward statics interactive simulation \n of a 3-segment tensegrity manipulator", textsize = 45)
    rowsize!(fig.layout, 2, Relative(4/6))
    colsize!(fig.layout, 1, Relative(1/2))
    fig
end

fig_video = get_actuate_video(man_video)

includet("../geodesicDR.jl")
function compare(nseg)
    k = 558.16; c = 100.0;
    # k = 1e-30; c = 100.0;
    restlen = 0.163-2.7246/318.16
    man_compare = man_ndof(nseg,[1.0,0.0];θ=0.0,k,c,restlen,isvirtual=true)
    # @show TR.get_strings_restlen(man_compare)
    @show TR.get_strings_stiffness(man_compare)
    # F = reshape(TR.build_G(man_compare.tg),:,1)
    start_sol,parameters0 = TR.get_start_system(man_compare)
    # start_sol.λ ./=1e2
    u1 = copy(parameters0.u)
    g1 = copy(parameters0.g)
    a1 = -0.01
    u1[[1,3,5]] .+= a1
    u1[[2,4,6]] .-= a1
    parameters1 = (u=u1, g=g1)

    # @info "Testing DR"
    # TR.set_restlen!(man_compare.tg,u1)
    # @info "Testing for ϵ=1e-13"
    # GDR(man_compare;β=1e-1,maxiters=10000,ϵ=1e-13,N=50,ξ=1e-13)
    # rs13,_,_,spine_DR13 = @btime $GDR($man_compare;β=5e-1,maxiters=10000,ϵ=1e-13,N=50,ξ=1e-13)
    # @info "Step = $(length(rs13)); Residual = $(rs13[end])"
    # # @show TR.get_strings_stiffness(spine_DR13)
    # # spine_DR13
    @info "Testing HC"
    @info "Precompiling"
    _ = TR.forward_once(man_compare.tg,start_sol,parameters0,parameters1,TR.PrimalMode())
    @info "Benchmarking"
    sol2 = @btime $TR.forward_once($man_compare.tg,$start_sol,$parameters0,$parameters1,$TR.PrimalMode())
    @info "Step = $(sol2.nstep); Residual = $(sol2.res)"
    TR.set_new_initial!(man_compare,sol2.q)
    man_compare
    # sol2
    # qdiff = norm(spine_DR13.traj.qs[end] .- sol2.q) |> maximum
    # @info "Recheck: $qdiff"
end

spine_DR = compare(4)
sol_HC = compare(4)
TR.set_new_initial!(spine_DR,sol_HC.q)

plot_traj!(spine_DR)

function compare_loop(irange)
    for i = irange
        @info "Testing $i-segment."
        compare(i)
    end
end

compare_loop(3:3:12)
compare_loop([12])

""
