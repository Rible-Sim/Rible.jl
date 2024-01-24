bot1 = deepcopy(bot)
q_bot = bot.traj.q[901]
s_bot = bot.traj.s[901]
bot1.tg.state.system.q̇ .= 0.0
bot1.traj.q[1] .= q_bot
bot1.traj.s[1] .= s_bot
TR.update!(bot1.tg)
TR.update_orientations!(bot1.tg)

plot_tg!(bot1.tg)

prob2 = TR.SimProblem(bot1, dynfuncs1)
dt = 1e-2; T = 6.0

TR.solve!(prob2,TR.FBZhong06(),
        (
            actuate! = actuate2!,
            prescribe! = nothing
            );
        dt=dt,tspan=(0.0,T),ftol=1e-7,verbose=true)

with_theme(theme_pub;
    Axis3 = (
        azimuth = 3π/2,
        elevation = 0.0,
    ),
    Mesh = (
        # color = :black,
        transparency = false,
    ),
    ) do
plot_traj!(bot1;
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
    xlims=(-100,3500),
    ylims=(-50,50),
    zlims=(-1800,1800),
)
end