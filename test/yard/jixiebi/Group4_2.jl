botc1 = deepcopy(bot1)
botc2 = deepcopy(bot2)
botc3 = deepcopy(bot3)
botc4 = deepcopy(bot4)
nt = 601
q_bot1 = bot1.traj.q[nt]
s_bot1 = bot1.traj.s[nt]
q_bot2 = bot2.traj.q[nt]
s_bot2 = bot2.traj.s[nt]
q_bot3 = bot3.traj.q[nt]
s_bot3 = bot3.traj.s[nt]
q_bot4 = bot4.traj.q[nt]
s_bot4 = bot4.traj.s[nt]
botc1.tg.state.system.q̇ .= 0.0
botc2.tg.state.system.q̇ .= 0.0
botc3.tg.state.system.q̇ .= 0.0
botc4.tg.state.system.q̇ .= 0.0
botc1.traj.q[1] .= q_bot1
botc1.traj.s[1] .= s_bot1
botc2.traj.q[1] .= q_bot2
botc2.traj.s[1] .= s_bot2
botc3.traj.q[1] .= q_bot3
botc3.traj.s[1] .= s_bot3
botc4.traj.q[1] .= q_bot4
botc4.traj.s[1] .= s_bot4
TR.update!(botc1.tg)
TR.update_orientations!(botc1.tg)
TR.update!(botc2.tg)
TR.update_orientations!(botc2.tg)
TR.update!(botc3.tg)
TR.update_orientations!(botc3.tg)
TR.update!(botc4.tg)
TR.update_orientations!(botc4.tg)

probc1 = TR.SimProblem(botc1,dynfuncs1)
probc2 = TR.SimProblem(botc2,dynfuncs1)
probc3 = TR.SimProblem(botc3,dynfuncs1)
probc4 = TR.SimProblem(botc4,dynfuncs1)

dt = 1e-2; T = 6.0

TR.solve!(probc1, TR.FBZhong06(),
        (
            actuate! = actuate_3!,
            prescribe! = nothing
            );
        dt=dt,tspan=(0.0,T),ftol=1e-7,verbose=true)
TR.solve!(probc2, TR.FBZhong06(),
        (
            actuate! = actuate_2!,
            prescribe! = nothing
            );
        dt=dt,tspan=(0.0,T),ftol=1e-7,verbose=true)
TR.solve!(probc3, TR.FBZhong06(),
        (
            actuate! = actuate_2!,
            prescribe! = nothing
            );
        dt=dt,tspan=(0.0,T),ftol=1e-7,verbose=true)
TR.solve!(probc4, TR.FBZhong06(),
        (
            actuate! = actuate_3!,
            prescribe! = nothing
            );
        dt=dt,tspan=(0.0,T),ftol=1e-7,verbose=true)



with_theme(theme_pub;
    Axis3 = (
        # azimuth = -2π/3,
        # elevation = -2π/3,
        azimuth = -2.1043951023931977,
        elevation = 0.7902036732051045
        # azimuth = -pi/2,
        # elevation =  - pi/4
    ),
    Mesh = (
        # color = :black,
        transparency = false,
    ),
    ) do
plot_traj_with_jixiepintai!(botc1,botc2,botc3,botc4;
    # AxisType=LScene,
    # doslide=false,
    doslide=true,
    showinfo=true,
    axtitle=false,
    # hidezdecorations=true,
    hideydecorations=true,
    fontsize=10,
    # gridsize=(2,2),
    # attimes=[0.0,2.25,3.75,6.0],
    # atsteps=nothing,
    showground=false,
    showlabels=false,
    rigidcolor=:grey,
    xlims=(-1600,2700),
    ylims=(-800,800),
    zlims=(-1400,1400),
    cylinder_x = 1300,
    cylinder_r = 700,
)
end
