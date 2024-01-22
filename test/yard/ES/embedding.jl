using Revise #jl
import Rible as RB
include(joinpath(pathof(RB),"../../test/yard/tensegrity.jl"))
using AbbreviatedStackTraces #jl
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true #jl
ENV["JULIA_STACKTRACE_MINIMAL"] = true #jl

figdir::String = joinpath(pathof(RB),"../../tmp")
if Sys.iswindows() #src
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\DynamicTensegrity\ES" #src
elseif Sys.isapple() #src
    figdir::String = raw"/Users/jacob/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/Papers/DynamicTensegrity/ES" #src
end #src

include(joinpath(pathof(RB),"../../test/vis.jl"))
includet(joinpath(pathof(RB),"../../test/vis.jl")) #jl

include(joinpath(pathof(RB),"../../examples/bodies/make_3d_tri.jl"))
includet(joinpath(pathof(RB),"../../examples/bodies/make_3d_tri.jl")) #jl
include(joinpath(pathof(RB),"../../examples/bodies/make_3d_bar.jl"))
includet(joinpath(pathof(RB),"../../examples/bodies/make_3d_bar.jl")) #jl
include(joinpath(pathof(RB),"../../examples/bodies/make_3d_plate.jl"))
includet(joinpath(pathof(RB),"../../examples/bodies/make_3d_plate.jl")) #jl
include(joinpath(pathof(RB),"../../examples/robots/twotre3d.jl"))
includet(joinpath(pathof(RB),"../../examples/robots/twotre3d.jl")) #jl
include(joinpath(pathof(RB),"../../examples/robots/embed3d.jl"))
includet(joinpath(pathof(RB),"../../examples/robots/embed3d.jl")) #jl


fontsize = 8 
markersize = fontsize
linewidth = 0.5 
cablewidth = 0.75 
barwidth = 1.5 
cw = 469 
tw = 469 

#-- embedding
#without outer and gravity, to obtain rest length
gravity = false

m = 3
α = 2π/m
θ = 1.25α
n = 4
b = 0.14
r = 0.04*sqrt(2)

twotre1 = twotre3d(;
    r1= 0.03*sqrt(2),
    r,b,m,α,θ,n = 1
) 

prism1 = embed3d(;
    r1= 0.03*sqrt(2),
    r,b,m,α,θ,n = 1,
    isprism = true,
)

newembed11 = embed3d(;
    r1= 0.03*sqrt(2),
    r,b,m,α,θ,n = 1,
    outer = true,
    # isprism = true,
)

newembed1o = embed3d(;
    r1= 0.03*sqrt(2),
    r,b,m,α,θ,n,
    outer = true,
)

prism0 = embed3d(;
    r1= 0.06*sqrt(2),
    r,b,m,α,θ,n = 1,
    isprism = true,
)

newembed0 = embed3d(;
    r1= 0.06*sqrt(2),
    r,b,m,α,θ,n = 1,
    outer = true,
    # isprism = true,
)

plot_traj!(newembed11)
GM.activate!(); with_theme(theme_pub;
        size = (0.6cw,0.2cw),
        figure_padding = (0,0,-fontsize/2,0),
        Axis3 = (
            azimuth = 3.907532586451984,
            elevation = 0.18379626575974045
        )
    ) do
    fig = Figure()
    xlims = (-1.8e-1,1.8e-1)
    ylims = (-1.8e-1,1.8e-1)
    zlims = (-1e-3,2.4e-1)
    doslide = false
    showpoints = false
    showlabels = false
    g0 = fig[1,1] = GridLayout(;tellheight=false,)
    g1 = fig[1,2] = GridLayout(;tellheight=false,)
    # g11 = g1[1,1] = GridLayout()
    # g12 = g1[2,1] = GridLayout()
    g2 = fig[1,3] = GridLayout(;tellheight=false,)
    # g21 = g2[1,1] = GridLayout()
    # g22 = g2[2,1] = GridLayout()
    g3 = fig[1,4] = GridLayout()
    plot_traj!(
        twotre1;
        AxisType = Axis3,
        fig = g0,
        xlims,
        ylims,
        zlims,
        doslide,
        showpoints,
        showlabels,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("(a) ", font=:bold),
            )
        end,
        sup! = (ax,tgob,sgi) -> begin
            hidexyz(ax)
        end,
    )
    plot_traj!(
        prism1;
        AxisType = Axis3,
        fig = g1,
        xlims,
        ylims,
        zlims,
        doslide,
        showpoints,
        showlabels,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("(b) ", font=:bold),
            )
        end,
        sup! = (ax,tgob,sgi) -> begin
            hidexyz(ax)
        end,
    )
        # plot_traj!(
        #     prism0;
        #     AxisType = Axis3,
        #     fig = g12,
        #     xlims,
        #     ylims,
        #     zlims = (-1e-3,1e-1),
        #     doslide,
        #     showpoints,
        #     showlabels,
        #     titleformatfunc = (sgi,tt)-> begin
        #         rich(
        #             rich("(d) ", font=:bold),
        #         )
        #     end,
        #     sup! = (ax,tgob,sgi) -> begin
        #         hidexyz(ax)
        #     end,
        # )
    plot_traj!(newembed11;
        AxisType = Axis3,
        fig = g2,
        xlims,
        ylims,
        zlims,
        doslide,
        showpoints,
        showlabels,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("(c) ", font=:bold),
            )
        end,
        sup! = (ax,tgob,sgi) -> begin
            hidexyz(ax)
        end,
    )
        # plot_traj!(newembed0;
        #     AxisType = Axis3,
        #     fig = g22,
        #     xlims,
        #     ylims,
        #     zlims = (-1e-3,1e-1),
        #     doslide,
        #     showpoints,
        #     showlabels,
        #     titleformatfunc = (sgi,tt)-> begin
        #         rich(
        #             rich("(e) ", font=:bold),
        #         )
        #     end,
        #     sup! = (ax,tgob,sgi) -> begin
        #         hidexyz(ax)
        #     end,
        # )
    plot_traj!(newembed1o;
        AxisType = Axis3,
        fig = g3,
        xlims,
        ylims,
        zlims = (-1e-3,8.5e-1),
        doslide,
        showpoints,
        showlabels,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("(d) ", font=:bold),
            )
        end,
        sup! = (ax,tgob,sgi) -> begin
            hidexyz(ax)
        end,
    )
    colgap!(fig.layout,0)
    rowsize!(g0,1,Fixed(0.10cw))
    rowsize!(g1,1,Fixed(0.10cw))
    rowsize!(g2,1,Fixed(0.10cw))
    # rowgap!(g1,-fontsize)
    # rowsize!(g1,2,Relative(0.35))
    # rowgap!(g2,-fontsize)
    # rowsize!(g2,2,Relative(0.35))
    savefig(fig,"embedding")
    fig
end
newembed1 = embed3d(;
    r1= 0.03*sqrt(2),
    r,b,m,α,θ,n,
    outer = false,
)
plot_traj!(newembed1)
μ1 = RB.inverse_for_restlength(
    newembed1,newembed1;gravity,fmin=100.0
)
barplot(μ1)
# (optional) check μ1 
RB.set_restlen!(newembed1.st,μ1)
RB.update!(newembed1.st)
f = RB.get_cables_tension(newembed1.st) 
f |> extrema
RB.check_static_equilibrium_output_multipliers(newembed1.st;gravity)
RB.undamped_eigen(newembed1.st;gravity)
RB.undamped_eigen!(newembed1;gravity)
RB.reset!(newembed1)
#-- without outer, see see gravity
RB.GDR!(newembed1;gravity=true)

plot_traj!(newembed1)

#-- folded without outer and gravity, to obtain rest length
newembed0 = embed3d(;
    r1= 0.06*sqrt(2),
    r,b,m,α,θ,n
)

plot_traj!(newembed0)

μ0 = RB.inverse_for_restlength(
    newembed0,newembed0;gravity,fmin=100.0
)
barplot(μ0)
#-- with outer, using the folded and unfolded 
outer = true
gravity = true
newembed0_outer = embed3d(;
    r1= 0.06*sqrt(2),
    r,b,m,α,θ,n,outer = true
)
plot_traj!(newembed0_outer)
μ0o = RB.get_cables_restlen(newembed0_outer)
μ0o[begin:3n*m] .= μ0
μ1o = deepcopy(μ0o)
μ1o[begin:3n*m] .= μ1
RB.set_restlen!(newembed0_outer.st,μ0o)
RB.GDR!(newembed0_outer;gravity)
GM.activate!();plot_traj!(newembed0_outer;)
RB.set_new_initial!(newembed0_outer,newembed0_outer.traj.q[end])
RB.update!(newembed0_outer.st;gravity)
RB.check_static_equilibrium_output_multipliers(newembed0_outer.st;gravity)
RB.undamped_eigen(newembed0_outer.st;gravity)
RB.undamped_eigen!(newembed0_outer;gravity)
friction_coefficients = [μ0o,μ1o]

barplot(μ0o)
barplot!(μ1o)
pretty_table(
    Dict(
        [
            ("Folded, DistanceSpringDampers 1~6", μ0o[1]),
            ("Deployed, DistanceSpringDampers 1~6", μ0o[7]),
            ("Folded, DistanceSpringDampers 1~6", μ1o[1]),
            ("Deployed, DistanceSpringDampers 1~6", μ1o[7]),
            ("Outer DistanceSpringDampers 1~6", μ1o[end-2:end]),
        ]
    )
)

Td = 8.0
tend = 10.0

friction_coefficients = [
    begin
        μ = deepcopy(μ0o)
        μ[begin:3m*j] .= μ1[begin:3m*j]
        μ
    end
    for j = 0:n
]

function make_multi_stages_pres_actor(friction_coefficients;start = 0.0, stop = 10.0, len = 2, )
    nμ = length(friction_coefficients[begin])

    function itp(t)
        scaled_itps = extrapolate(
            Interpolations.scale(
                interpolate(
                    reduce(hcat,friction_coefficients),
                    (NoInterp(),BSpline(Linear()))
                    # (NoInterp(),BSpline(Quadratic(Flat(OnGrid()))))
                ),
                1:nμ, range(;start, stop, length = len,)
            ),
            (Throw(),Flat())
        )
        [scaled_itps(j,t) for j in 1:nμ]
    end

    RB.PrescribedActuator(
        1,
        RB.ManualActuator(1,collect(1:nμ),zeros(nμ),RB.Uncoupled()),
        itp
    )
end
make_multi_stages_pres_actor(friction_coefficients;start=0.0,stop=Td,len=n+1)

newembed0_outer_deploy = RB.Robot(
    deepcopy(newembed0_outer.st),
    (actuators=[make_multi_stages_pres_actor(friction_coefficients;start=0.0,stop=Td,len=length(friction_coefficients))],)
)

RB.solve!(
    RB.DynamicsProblem(
        newembed0_outer_deploy,
        (x)->dynfuncs(x;actuate=true,gravity=true)
    ),
    RB.Zhong06();
    dt=1e-3,
    tspan=(0.0,tend),
    ftol=1e-7,
    exception=true
)

plot_traj!(newembed0_outer_deploy)
GM.activate!();with_theme(theme_pub;
        figure_padding = (0,fontsize,-fontsize/2,fontsize/2),
        # size = (0.9cw,0.25cw),
        Axis3 = (
            azimuth = 4.995530633326985,
            elevation = 0.18269908169872415
        )
    ) do
    plot_traj!(
        newembed0_outer_deploy;
        AxisType=Axis3,
        figsize = (0.7cw,0.25cw),
        gridsize = (1,5),
        attimes = range(0,Td,5),
        xlims = (-2e-1,2e-1),
        ylims = (-4e-1,2e-1),
        zlims = (-1e-3,8e-1),
        showinfo = false,
        doslide = false,
        # dorecord = true,
        slack_linestyle = :solid,
        actuate = true,
        showpoints = false,
        showlabels = false,
        sup! = (ax,tgob,sgi) -> begin
            hidex(ax)
            hidey(ax)
            if sgi != 5
                hidez(ax)
            end
        end,
        # rowgap=2fontsize,
        colgap=0,
        figname = "newembed0_outer_deploy"
    )
end
nb = newembed0_outer_deploy.st.nbodies
ps = [RB.get_trajectory!(newembed0_outer_deploy,nb,j) for j = m+1:m+m]
GM.activate!();with_theme(theme_pub;
        size = (0.7cw,0.2cw),
        figure_padding = (0,fontsize,0,0)
    ) do 
    fig = Figure()
    (;t) = newembed0_outer_deploy.traj
    for i = 1:m
        ax = Axis(fig[1,i], 
            xlabel = tlabel, 
            ylabel = latexstring("""$(["x","y","z"][i])~(\\mathrm{m})""")
        )
        lines!(ax,t,ps[2][i,:])
        # lines!(ax,t,ps[2][i,:])
        # lines!(ax,t,ps[3][i,:])
        xlims!(t[begin],t[end])
        Label(
            fig[1,i,TopLeft()],
            "($(alphabet[i]))",
            font = "CMU Serif Bold",
            padding = (0, 0, 0, 0),
            halign = :right,
            valign = :bottom
        )
    end
    savefig(fig,"newembed0_outer_deploy_traj")
    fig
end
