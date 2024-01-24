
bar = make_bar(;)

plot_traj!(bar;)

function bar_contact_dynfuncs(bot)
end


tspan = (0.0,0.4)
h = 1e-3


prob = RB.DynamicsProblem(bar,bar_contact_dynfuncs)

RB.solve!(prob,RB.ZhongCCP();tspan,dt=h,ftol=1e-12,maxiters=50,exception=false)

contacts_traj_voa = VectorOfArray(bar.contacts_traj)
for (i,c) in enumerate(contacts_traj_voa[1,1:end])
 check_Coulomb(i,c)
end

plot_traj!(bar;showinfo=false)

me = RB.mechanical_energy!(bar)
lines(me.E)

[c[1].state.active for c in bar.contacts_traj] |> lines
[c[1].state.persistent for c in bar.contacts_traj] |> lines

[c[2].state.active for c in bar.contacts_traj] |> lines
[c[2].state.persistent for c in bar.contacts_traj] |> lines