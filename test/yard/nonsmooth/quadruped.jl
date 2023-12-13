
#----------- quadruped ---------------
quadbot = quad(10.0)
quadbot.st.connectivity.numbered

plot_traj!(quadbot;
    zlims = (0,1),
    showground = false
)

function quad_dynfuncs(bot)
    contact_dynfuncs(bot;
        flatplane = RB.Plane([0,0,1.0],[0,0,0.0])
    )
end

tspan = (0.0,1.91)
tspan = (0.0,1.0)
h = 1e-3

prob = RB.DynamicsProblem(quadbot,quad_dynfuncs)

RB.solve!(prob,RB.ZhongCCP();tspan,dt=h,ftol=1e-14,maxiters=50,exception=true)

RB.solve!(prob,RB.AlphaCCP(0.8);tspan,dt=h,ftol=1e-10,maxiters=50,exception=true)

plot_traj!(quadbot;
    zlims = (0,1),
    showground = false
)

me = RB.mechanical_energy!(quadbot)

me.E |> lines
contacts_traj_voa = VectorOfArray(quadbot.contacts_traj)
