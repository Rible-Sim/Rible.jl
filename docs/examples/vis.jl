function plot_local_points(rb::TR.AbstractRigidBody)
    (;r̄g,r̄ps,nr̄ps) = rb.prop
    fig = Figure(resolution=(1280,720))
    ndim = TR.get_ndim(rb)
    if ndim == 2
        ax = Axis(fig[1,1])
        ax.aspect = DataAspect()
    else
        ax = LScene(fig[1,1], scenekw = (clear=true,))
    end
    scatter!(ax,r̄ps,markersize=5)
    scatter!(ax,r̄g,markersize=5)
    text!(ax,
        ["$i $(string(r̄p))" for (i,r̄p) in enumerate(r̄ps)] ,
        position = r̄ps,
        color = :darkblue,
        align = (:left, :baseline),
        offset = (-5, 10)
    )
    text!(ax,
        "r̄g $(string(r̄g))",
        position = r̄g,
        color = :darkred,
        align = (:left, :baseline),
        offset = (-5, 10)
    )
    fig
end

function plot_tg!(tg::TR.TensegrityStructure)
    ndim = TR.get_ndim(tg)
    fig = Figure(resolution=(1280,720))
    if ndim == 2
        ax = Axis(fig[1,1], aspect=:data)
    else
        ax = Axis3(fig[1,1], aspect=:data)
    end
    tgob = Observable(deepcopy(tg))
    init_plot!(ax,tgob)
    fig
end

function plot_traj!(bot::TR.TensegrityRobot;showlabels=true,textsize=10,actuate=false)
    (;tg,traj) = bot
    ndim = TR.get_ndim(tg)
    fig = Figure(resolution=(1280,720))
    axtitle = Observable("")
    if ndim == 2
        ax = Axis(fig[1,1],title=axtitle)
        ax.aspect = DataAspect()
    else
        # ax = Axis3(fig[1,1],title=axtitle, aspect=:data)
        # ax = LScene(fig[1,1], show_axis=false, scenekw = (clear=true,))
        ax = LScene(fig[1,1], ) #scenekw = (clear=true,)
        # cam = Makie.camera(ax.scene)
        # cam = cam3d!(ax.scene, projectiontype=Makie.Orthographic)
        # update_cam!(ax.scene, cam)
    end
    grid1 = fig[2,1] = GridLayout()
    tg.state.system.q .= traj.q[begin]
    TR.update!(tg)
    tgob = Observable(deepcopy(tg))
    init_plot!(ax,deepcopy(tgob);isref=true)
    init_plot!(ax,tgob;showlabels,textsize)
    ls_step = labelslider!(fig,"Step",1:length(traj))
    grid1[1,1] = ls_step.layout
    on(ls_step.slider.value) do this_step
        tg.state.system.q .= traj.q[this_step]
        # tg.state.system.q̇ .= traj.q̇[this_step]
        if actuate
            TR.actuate!(bot,[traj.t[this_step]])
        end
        TR.update!(tg)
        axtitle[] = string(traj.t[this_step])
        # @show TR.mechanical_energy(tg)
        #
        # camera = cameracontrols(ax.scene)
        # lookat = camera.lookat[]
        # eyepos = camera.eyeposition[]
        # @show lookat, eyepos
        TR.analyse_slack(tg,true)
        tgob[] = tg
    end
    fig
end

function init_plot!(ax,tgob;isref=false,showlabels=true,textsize=10)
    (;ncables,nrigids) = tgob[]
    if isref
        showlabels = false
        cablecolor =
        cablelabelcolor =
        rigidcolor =
        rigidlabelcolor = :grey
    else
        cablecolor=:deepskyblue
        # cablecolor=:red
        cablelabelcolor=:darkgreen
        rigidcolor=:black
        rigidlabelcolor=:darkblue
    end
    # cables
    if ncables !== 0
    # if false
        linesegs_cables = @lift begin
            (;connected) = $tgob.connectivity
            ndim = TR.get_ndim($tgob)
            T = TR.get_numbertype($tgob)
            ret = Vector{Pair{Point{ndim,T},Point{ndim,T}}}()
            mapreduce(
                (scnt)->
                Point(scnt.end1.rbsig.state.rps[scnt.end1.pid]) =>
                Point(scnt.end2.rbsig.state.rps[scnt.end2.pid]),
                vcat,
                connected
                ;init=ret
            )
        end
        linesegments!(ax, linesegs_cables, color = cablecolor, linewidth = 2)
        rcs_by_cables = @lift begin
            (;connected) = $tgob.connectivity
            ndim = TR.get_ndim($tgob)
            T = TR.get_numbertype($tgob)
            ret = Vector{MVector{ndim,T}}()
            mapreduce(
                (scnt)->
                [(
                    scnt.end1.rbsig.state.rps[scnt.end1.pid].+
                    scnt.end2.rbsig.state.rps[scnt.end2.pid]
                )./2],
                vcat,
                connected
                ;init=ret
            )
        end
        # @show rcs_by_cables
        if showlabels
            text!(ax,
                ["c$(i)" for (i,rc) in enumerate(rcs_by_cables[])] ,
                position = rcs_by_cables,
                textsize = textsize,
                color = cablelabelcolor,
                align = (:left, :top),
                offset = (-5, -10)
            )
        end
    end
    # rigid bars
    rigidbars = TR.get_rigidbars(tgob[])
    if !isempty(rigidbars)
        linesegs_bars = @lift begin
            bars = TR.get_rigidbars($tgob)
            ndim = TR.get_ndim($tgob)
            T = TR.get_numbertype($tgob)
            ret = Vector{Pair{Point{ndim,T},Point{ndim,T}}}()
            mapreduce(
                (bar)->
                Point(bar.state.rps[1]) =>
                Point(bar.state.rps[2]),
                vcat,
                bars
                ;init=ret
            )
        end
        linesegments!(ax, linesegs_bars, color = rigidcolor, linewidth = 4)
    end
    if showlabels
        # mass centers
        rg_by_rbs = @lift begin
            rbs = TR.get_rigidbodies($tgob)
            [rb.state.rg for rb in rbs]
        end
        scatter!(ax,rg_by_rbs, color = rigidlabelcolor, markersize = 4)
        # points on rigidbodies
        rps_by_rbs = [
            @lift begin
                rbs = TR.get_rigidbodies($tgob)
                rbs[rbid].state.rps
            end
            for rbid = 1:nrigids
        ]
        text!(ax,
            ["r$(i)g" for (i,rg) in enumerate(rg_by_rbs[])] ,
            position = rg_by_rbs,
            textsize = textsize,
            color = rigidlabelcolor,
            align = (:left, :top),
            offset = (-5, -10)
        )
        for (rbid,rps) in enumerate(rps_by_rbs)
            # if rbid in [8,9]
                scatter!(ax,rps)
                text!(ax,
                    # ["r$(rbid)p$pid $(string(rp))" for (pid,rp) in enumerate(rps[])],
                    ["r$(rbid)p$pid" for (pid,rp) in enumerate(rps[])],
                    position = rps,
                    textsize = textsize,
                    color = :darkred,
                    align = (:left, :top),
                    offset = (20(rand()-0.5), 20(rand()-0.5))
                )
            # end
        end
    end
end

function get_trajectory!(bot::TR.TensegrityRobot,rbid::Int,pid::Int,step_range=:)
    (; tg, traj)= bot
    rp = VectorOfArray(Vector{Vector{Float64}}())
    rbs = TR.get_rigidbodies(tg)
    rb = rbs[rbid]
    for q in traj.q
        TR.update_rigids!(tg,q)
        if pid == 0
            push!(rp,rb.state.rg)
        else
            push!(rp,rb.state.rps[pid])
        end
    end
    rp
end

function get_velocity!(bot::TR.TensegrityRobot,rbid::Int,pid::Int,step_range=:)
    (; tg, traj)= bot
    ṙp = VectorOfArray(Vector{Vector{Float64}}())
    rbs = TR.get_rigidbodies(tg)
    rb = rbs[rbid]
    for (q,q̇) in zip(traj.q, traj.q̇)
        TR.update_rigids!(tg,q,q̇)
        push!(ṙp,rb.state.ṙps[pid])
    end
    ṙp
end

function get_angular!(bot::TR.TensegrityRobot,rbid::Int,step_range=:)
    (; tg, traj)= bot
    ω = VectorOfArray(Vector{Vector{Float64}}())
    rbs = TR.get_rigidbodies(tg)
    rb = rbs[rbid]
    for (q,q̇) in zip(traj.q, traj.q̇)
        TR.update_rigids!(tg,q,q̇)
        TR.update_orientations!(tg)
        push!(ω,rb.state.ω)
    end
    ω
end

function get_time_mids(bot::TR.TensegrityRobot)
    (;t) = bot.traj
    (t[begin+1:end] .+ t[begin:end-1])./2
end

function get_tension!(bot::TR.TensegrityRobot,cid::Int,step_range=:)
    (; tg, traj)= bot
    (; cables) = tg
    T = TR.get_numbertype(tg)
    f = Vector{T}()
    h = traj.t[2] - traj.t[1]
    q_mids = [(traj.q[k] .+ traj.q[k+1])./2 for k = 1:length(traj)-1]
    q̇_mids = [(traj.q[k] .- traj.q[k+1])./h for k = 1:length(traj)-1]
    for (q,q̇) in zip(q_mids, q̇_mids)
        TR.update!(tg,q)
        push!(f,cables[cid].state.tension)
    end
    f
end

function parse_Adams_dynamic(url)
    adams_xml = readxml(url)
    ns = namespaces(adams_xml.root)
    ns[1] = "adams" => ns[1][2]
    entities = findall("//adams:Entity[@entType='ADAMS_Variable']", adams_xml.root, ns)
    names = Symbol.(vcat("time",[replace(ele["entity"],"."=>"_") for ele in entities]))
    step_initialConditions = findfirst("//adams:Step[@type='initialConditions']",adams_xml.root, ns)
    initialConditions = (;zip(names,parse.(Float64,split(step_initialConditions.content,"\n")[2:end-1]))...)
    ret = StructArray([initialConditions])
    step_dynamic = findall("//adams:Step[@type='dynamic']",adams_xml.root, ns)
    for step in step_dynamic
        step_values = parse.(Float64,split(step.content,"\n")[2:end-1])
        # @assert length(variable_names) == length(step_data)
        push!(ret,(;zip(names,step_values)...))
    end
    ret
end

function parse_Adams_static(url)
    adams_xml = readxml(url)
    ns = namespaces(adams_xml.root)
    ns[1] = "adams" => ns[1][2]
    entities = findall("//adams:Entity[@entType='ADAMS_Variable']", adams_xml.root, ns)
    names = Symbol.(vcat("time",[replace(ele["entity"],"."=>"_") for ele in entities]))
    step_initialConditions = findfirst("//adams:Step[@type='initialConditions']",adams_xml.root, ns)
    step_static = findfirst("//adams:Step[@type='static']",adams_xml.root, ns)
    step_values = parse.(Float64,split(step_static.content,"\n")[2:end-1])
    (;zip(names,step_values)...)
end

function parse_Adams_frequency(url)
    adams_xml = readxml(url)
    ns = namespaces(adams_xml.root)
    ns[1] = "adams" => ns[1][2]
    frequency_node = findfirst("//adams:ModalEnergyHeader[@type='undamped_natural_freqs']",adams_xml.root, ns)
    frequency_string = split(frequency_node.content,"\n")[2]
    frequency = parse.(Float64,split(frequency_string," "))
end
