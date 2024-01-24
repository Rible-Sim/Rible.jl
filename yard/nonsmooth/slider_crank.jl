using Revise #jl
import Rible as RB
include(joinpath(pathof(RB),"../../yard/nonsmooth.jl"))
using AbbreviatedStackTraces #jl
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true #jl
ENV["JULIA_STACKTRACE_MINIMAL"] = true #jl
import Rible as RB
figdir::String = joinpath(pathof(RB),"../../tmp")
if Sys.iswindows() #src
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\IMSD 2025\LaTex_Abstract" #src
elseif Sys.isapple() #src
    figdir::String = raw"/Users/jacob/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/Papers/IMSD 2024/LaTex_Abstract" #src
end #src
include(joinpath(pathof(RB),"../../test/vis.jl"))
includet(joinpath(pathof(RB),"../../test/vis.jl")) #jl

tw = 455.8843 #pt |> pt2px
scalefactor = 4

#--  slider crank 
include(joinpath(pathof(RB),"../../examples/robots/slider_crank.jl"))
includet(joinpath(pathof(RB),"../../examples/robots/slider_crank.jl"))#jl

# natural coordinates
sc = slider_crank(;coordsType=RB.NCF.NC)

# Quaternion coordinates
sc = slider_crank(;coordsType=RB.QCF.QC)

# Contact Surfaces
planes = RB.StaticContactSurfaces(
    [
        RB.Plane([0,0, 1.0],[0,0,-0.026]),
        RB.Plane([0,0,-1.0],[0,0, 0.026])
    ]
)

GM.activate!(;scalefactor); plot_traj!(sc;showground=false)

RB.has_constant_mass_matrix(sc)

dt = 5e-4
tspan = (0.0,2*683dt)
tspan = (0.0,1.0)

# No Contact Dynamics
sc = slider_crank(;coordsType=RB.NCF.NC)

prob = RB.DynamicsProblem(sc,)
RB.solve!(
    prob,
    RB.DynamicsSolver(
        RB.Zhong06()
        # RB.Moreau(0.5)
    );
    dt,tspan,ftol=1e-12,maxiters=50,verbose=true,exception=true,progress=false,
)
GM.activate!(;scalefactor=1); plot_traj!(sc;showground=false)

RB.solve!(
    prob,
    RB.DynamicsSolver(RB.Moreau(0.5));
    dt,tspan,ftol=1e-12,maxiters=3,verbose=true,exception=true,progress=false,
)


dts = [
    1e-3,5e-4,2e-4,1e-4,1e-5
]

sc_dts = [
    RB.solve!(
        RB.DynamicsProblem(
            slider_crank(;coordsType=RB.NCF.NC),
        ),
        RB.DynamicsSolver(solver);
        dt,
        tspan=(0.0,0.1),
        ftol=1e-14,maxiters=50,verbose=false,exception=true,progress=true,
    ).prob.bot
    for dt in dts, solver in (RB.Zhong06(), RB.Moreau(1.0)) 
]

fig = Figure() #src
ax = Axis(fig[1,1]) #src
for bot in sc_dts[begin+1:end,1] #src
    b1r2 = RB.get_mid_velocity!(bot,4,2) #src
    lines!(ax,RB.get_mid_times(bot),b1r2[2,:]) #src
end #src
for bot in sc_dts[begin+1:end,2] #src
    b1r2 = RB.get_mid_velocity!(bot,4,2) #src
    lines!(ax,RB.get_mid_times(bot),b1r2[2,:]) #src
end #src
fig #src


GM.activate!(;scalefactor); with_theme(theme_pub;
        size = (0.9tw,0.2tw)
    ) do
    fig = Figure()
    ax1 = Axis(fig[1,1],ylabel = "Abs. Err. (rad)")
    ax2 = Axis(fig[1,2],ylabel = "Abs. Err. (rad/s)")
    _,err_avg_moreau_traj = RB.get_err_avg(vcat(sc_dts[begin:end-1,2],sc_dts[end,1]);bid=4,pid=2,di=2,field=:traj)
    plot_convergence_order!(ax1,dts[begin:end-1],err_avg_moreau_traj;show_orders=false,marker=:utriangle,color=:blue,label="Moreau")
    _,err_avg_nmsi_traj = RB.get_err_avg(sc_dts[:,1];bid=4,pid=2,di=2,field=:traj)
    plot_convergence_order!(ax1,dts[begin:end-1],err_avg_nmsi_traj,label="NMSI")
    _,err_avg_moreau_vel = RB.get_err_avg(vcat(sc_dts[begin:end-1,2],sc_dts[end,1]);bid=4,pid=2,di=3,field=:midvel)
    plot_convergence_order!(ax2,dts[begin:end-1],err_avg_moreau_vel;show_orders=false,marker=:utriangle,color=:blue,label="Moreau")
    _,err_avg_nmsi_vel = RB.get_err_avg(sc_dts[:,1];bid=4,pid=2,di=3,field=:midvel)
    plot_convergence_order!(ax2,dts[begin:end-1],err_avg_nmsi_vel,label="NMSI")
    # axislegend(ax2,position=:rb)
    Label(fig[1,1,TopLeft()],"($(alphabet[1]))",font=:bold)
    Label(fig[1,2,TopLeft()],"($(alphabet[2]))",font=:bold)
    # savefig(fig,"sc_convergence";backend=CM)
    fig
end

# No Contact Dynamics
prob = RB.DynamicsProblem(sc,)

RB.solve!(
    prob,
    RB.DynamicsSolver(RB.Zhong06());
    dt,tspan,ftol=1e-12,maxiters=50,verbose=true,exception=true,progress=false,
)

GM.activate!(;scalefactor = 1);plot_traj!(sc;showground=false)

me = RB.mechanical_energy!(sc)
lines(me.E)
# Frictionless Contact Dynamics

dt = 1e-4
tspan = (0.0,0.1)

prob = RB.DynamicsProblem(
    sc,
    planes,
    RB.RestitutionFrictionCombined(
        RB.NewtonRestitution(),
        RB.Frictionless(),
    )
)

# two-layer
simulator = RB.solve!(
    prob,
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        )
    );
    dt,tspan,ftol=1e-12,maxiters=50,verbose=true,exception=true,progress=false,
)

# Mono
simulator = RB.solve!(
    prob,
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.MonolithicContactSolver(
            RB.InteriorPointMethod()
        )
    );
    dt,
    tspan = (0.0,0.01),
    ftol=1e-12,maxiters=50,verbose_contact=true,exception=true,progress=false,
)

# comparison 
#at different timestep
# number of outer iterations: average, maximum, (contact)/(nocontact) 
# number of inner iterations: average maximum,
# total computation time/per timestep
dts = [
    1e-3,5e-4,2e-4,1e-4
]

sim_sc_iters = [
    @time RB.solve!(
        RB.DynamicsProblem(
            slider_crank(;coordsType=RB.QCF.QC),
            planes,
            RB.RestitutionFrictionCombined(
                RB.NewtonRestitution(),
                RB.Frictionless(),
            )
        ),
        RB.DynamicsSolver(
            RB.Zhong06(),
            contact_solver,
        );
        dt,
        tspan=(0.0,0.01),
        ftol=1e-12,maxiters=50,verbose=false,exception=true,progress=true,
    )
    for dt in dts, contact_solver in [
        # RB.InnerLayerContactSolver(
        #     RB.InteriorPointMethod()
        # ),
        RB.MonolithicContactSolver(
            RB.InteriorPointMethod()
        ),
    ]
]

function get_mean_and_max(iters)
    mean(iters), maximum(iters; init=-Inf)
end

function get_mean_and_min(iters)
    mean(iters), minimum(iters; init=Inf)
end

function get_two_itertaions_stats(simulator)
    (;data) = simulator.solver_history
    (;num_of_contacts,outer_iteration,inner_iteration,outer_condition_number,inner_condition_number) = data
    indices = [
        num_of_contacts .== i
        for i = 0:2
    ]
    outer_mean_max = map(
        get_mean_and_max,
        [
            outer_iteration[idx]
            for idx in indices
        ]
    ) 
    inner_mean_max = map(
        get_mean_and_max,
        [
            inner_iteration[idx]
            for idx in indices
        ]
    )
    outer_cd_mean_max = map(
        get_mean_and_max,
        [
            outer_condition_number[idx]
            for idx in indices
        ]
    )
    inner_cd_mean_max = map(
        get_mean_and_max,
        [
            inner_condition_number[idx]
            for idx in indices
        ]
    )
    @eponymtuple(outer_mean_max, inner_mean_max, outer_cd_mean_max, inner_cd_mean_max)
end

function get_mono_itertaions_stats(simulator)
    (;data) = simulator.solver_history
    (;num_of_contacts,iteration,stepsizes,condition_number) = data
    indices = [
        num_of_contacts .== i
        for i = 0:2
    ]
    mean_max_iter = map(
        get_mean_and_max,
        [
            iteration[idx]
            for idx in indices
        ]
    )
    mean_min_step1 = map(
        get_mean_and_min,
        [
            VectorOfArray(stepsizes[idx])[1,:]
            for idx in indices
        ]
    )
    # mean_min_step2 = map(
    #     get_mean_and_min,
    #     [
    #         VectorOfArray(stepsizes[idx])[2,:]
    #         for idx in indices
    #     ]
    # )
    cd_mean_max = map(
        get_mean_and_max,
        [
            condition_number[idx]
            for idx in indices
        ]
    )
    mean_max_iter, mean_min_step1, cd_mean_max
end

get_mono_itertaions_stats.(sim_sc_iters[:,1])

dts_strings = [
    "\$h=\\qty{$dt}{\\second}\$"
    for dt in [
        "1e-3", "5e-4", "2e-4", "1e-4"
    ]
]

stats = StructArray(get_two_itertaions_stats.(sim_sc_iters[:,1]))

t1 = TableCol("", dts_strings, string.(VectorOfArray(stats.outer_mean_max)[1,:]))
t2 = TableCol("outer iteration", dts_strings, string.(VectorOfArray(stats.outer_mean_max)[3,:]))
t3 = TableCol("inner iteration", dts_strings, string.(VectorOfArray(stats.inner_mean_max)[3,:]))

stats_table = join_table( "no contact"  =>hcat(t1),
            "contacting"  =>hcat(t2, t3))
to_tex(stats_table) |> print

stepsizes = VectorOfArray(simulator.solver_history.data.stepsizes)

# convergence

dts = [
    1e-3,5e-4,2e-4,1e-4,1e-5
]
# Quaternion coordinates
sc_mono_dt = [
    RB.solve!(
        RB.DynamicsProblem(
            slider_crank(;coordsType=RB.QCF.QC),
            planes,
            RB.RestitutionFrictionCombined(
                RB.NewtonRestitution(),
                RB.Frictionless(),
            )
        ),
        RB.DynamicsSolver(
            RB.Zhong06(),
            RB.InnerLayerContactSolver(
                RB.InteriorPointMethod()
            )
        );
        dt,
        tspan=(0.0,0.1),
        ftol=1e-12,maxiters=50,verbose=false,exception=true,progress=true,
    ).prob.bot
    for dt in dts
]
_, err_avg = RB.get_err_avg(sc_mono_dt;bid=4,pid=1,di=3)

GM.activate!(;px_per_unit=4,scalefactor = 4);with_theme(theme_pub;
        fonts = (; 
            regular = "CMU Serif", 
            bold = "CMU Serif Bold",
            italic = "CMU Serif Italic",
            math = "NewComputerModern 10 Italic"
        ),
        size = (tw,0.2tw),
        figure_padding = (fontsize,fontsize,0,-fontsize*2),
        Axis3 = (
            azimuth = 6.338616570826984,
            elevation = 0.1470350191987241,
        )
    ) do 
    fig = Figure()
    gl1 = fig[1,1] = GridLayout()
    gl2 = fig[1,2] = GridLayout()
    gl3 = fig[1,3] = GridLayout()
    r1 = RB.get_trajectory!(sc,4,1)
    r2 = RB.get_trajectory!(sc,4,2)
    r3 = RB.get_trajectory!(sc,4,3)
    r4 = RB.get_trajectory!(sc,4,4)
    thestep = RB.time2step(0.013,sc.traj.t)
    plot_traj!(
        sc,
        AxisType = Axis3,
        fig = gl1,
        atsteps = [thestep],
        showinfo = true,
        doslide=false,
        showground = false,
        showpoints = false,
        showtitle = false,
        xlims = (-0.05,0.05),
        ylims = (-0.05,0.6),
        zlims = (-0.05,0.15),
        sup! = (ax,tgob,sgi) -> begin
            hidexdecorations!(ax)
            scatter!(ax,
                Point3(r1[thestep]),
                markersize = fontsize/2,
            )
        end
    )
    # traj
    ax = Axis(gl2[1,1];
        xlabel=tlabel, ylabel = L"z~(\mathrm{m})"
    )
    lines!(ax,sc.traj.t,r1[3,:])
    xlims!(ax,extrema(sc.traj.t)...)

    ax3 = Axis(gl3[1,1];
        ylabel="Abs. Err. (m)"
    )
    plot_convergence_order!(ax3,dts[begin:end-1],err_avg;show_orders=true)
    # axislegend(ax3;position=:rb)

    Label(gl1[1,1,TopLeft()],"($(alphabet[1]))",font=:bold)
    Label(gl2[1,1,TopLeft()],"($(alphabet[2]))",font=:bold)
    Label(gl3[1,1,TopLeft()],"($(alphabet[3]))",font=:bold)
    colsize!(fig.layout,2,Relative(0.30))
    colsize!(fig.layout,3,Relative(0.23))
    colgap!(fig.layout,fontsize/2)
    colgap!(fig.layout,1,-fontsize*3)
   # ylims!(ax,-0.)
    savefig(fig,"slider_crank")
    fig
end

r1_dt1 = RB.get_trajectory!(sc_mono_dt[5],4,1)
r1_dt1[3,:] |> lines

# Frictional Contact Dynamics

dt = 1e-4
tspan = (0.0,1.0)

prob = RB.DynamicsProblem(
    sc,
    planes,
    RB.RestitutionFrictionCombined(
        RB.NewtonRestitution(),
        RB.CoulombFriction(),
    )
)

# two-layer
RB.solve!(
    prob,
    RB.DynamicsSolver(
        RB.Zhong06(),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        )
    );
    dt,tspan,ftol=1e-14,maxiters=50,verbose_contact=true,exception=true,progress=false,
)


# two-layer
RB.solve!(
    prob,
    RB.DynamicsSolver(
        RB.Moreau(1.0),
        RB.InnerLayerContactSolver(
            RB.InteriorPointMethod()
        )
    );
    dt,tspan,ftol=1e-14,maxiters=50,verbose=true,exception=true,progress=false,
)

GM.activate!(;px_per_unit=2,scalefactor = 2);plot_traj!(sc;showground=false)

me = RB.mechanical_energy!(sc)
Makie.lines(me.E)