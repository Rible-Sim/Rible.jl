with_theme(theme_pub;
    Axis3 = (
        azimuth = -π/2,
        elevation = 0.0,
        # azimuth = 5.812388980384689,
# elevation = 0.5100000000000002
    ),
    Mesh = (
        # color = :black,
        transparency = false,
    ),
    ) do
    record_traj!(bot, 1:450;
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
        xlims=(-100,3400),
        ylims=(-800,800),
        zlims=(-800,800),
        # cylinder_z = -420,
        # cylinder_x = 1550,
        # cylinder_r = 220,
    )
end 

with_theme(theme_pub;
    Axis3 = (
        # azimuth = 3π/2,
        # azimuth = 0.0,
        # elevation = π/2,
        # elevation = 0.0,
        # azimuth = π/4,
        # elevation = π/6,
    #     azimuth = -2.1043951023931977,
    # elevation = 0.7902036732051045,  
    azimuth = -1.244395102393198,
elevation = 0.5502036732051043
    ),
    Mesh = (
        # color = :black,
        transparency = false,
    ),
    ) do
    record_traj!(botc1,botc2,botc3,botc4, 336:475;
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
        xlims=(-1600,2700),
        ylims=(-800,800),
        zlims=(-1400,1400),
        cylinder_x = 1300,
        cylinder_r = 700,
    )
end

with_theme(theme_pub;
    Axis3 = (
        # azimuth = 3π/2,
        # azimuth = 0.0,
        # elevation = π/2,
        # elevation = 0.0,
        # azimuth = π/4,
        # elevation = π/6,
    #     azimuth = -2.1043951023931977,
    # elevation = 0.7902036732051045,  
    azimuth = -0.9943951023931974,
elevation = 0.7202036732051046
    ),
    Mesh = (
        # color = :black,
        transparency = false,
    ),
    ) do
    record_traj_jixiepingtai!(bot, 1:701;
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
        xlims=(-1600,2700),
        ylims=(-800,800),
        zlims=(-1400,1400),
        cylinder_x = 1300,
        cylinder_r = 700,
    )
end