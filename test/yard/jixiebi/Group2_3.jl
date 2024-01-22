bot2 = deepcopy(bot1)
q_bot1 = bot1.traj.q[601]
s_bot1 = bot1.traj.s[601]
q̇_bot1 = bot1.traj.q̇[601]
bot2.tg.state.system.q̇ .= q̇_bot1
bot2.traj.q̇[1] .= q̇_bot1
bot2.traj.q[1] .= q_bot1
bot2.traj.s[1] .= s_bot1
TR.update!(bot2.tg)
TR.update_orientations!(bot2.tg)

# plot_tg!(bot2.tg)

prob3 = TR.SimProblem(bot2, dynfuncs1)
dt = 1e-2; T = 9.0

TR.solve!(prob3,TR.FBZhong06(),
        (
            actuate! = actuate2_1!,
            prescribe! = nothing
            );
        dt=dt,tspan=(0.0,T),ftol=1e-7,verbose=true)

with_theme(theme_pub;
    Axis3 = (
        azimuth = -pi/2,
        elevation = 0.0,
        # azimuth = 5.812388980384689,
# elevation = 0.5100000000000002
    ),
    Mesh = (
        # color = :black,
        transparency = false,
    ),
    ) do
plot_traj_with_cylinder!(bot2;
# AxisType=LScene,
# doslide=false,
doslide=true,
showinfo=false,
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
xlims=(-1800,3200),
ylims=(-350,350),
zlims=(-1800,600),
cylinder_z = -420,
cylinder_x = 1550,
cylinder_r = 220,
)
end