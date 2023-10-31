
function set_ax!(ax,zup=0.20;elev=16,azim=-50)
    # ax.set_xticklabels([])
    # ax.set_yticklabels([])
    # ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    xmarg = ymarg = 0.0000
    zmarg = 0.000000
    x_min = -0.12+xmarg; x_max = 0.12-xmarg
    y_min = -0.12+ymarg; y_max = 0.12-ymarg
    z_min =  0.0+zmarg; z_max = zup-zmarg
    xspan = x_max - x_min
    yspan = y_max - y_min
    zspan = z_max - z_min
    ax.set_xlabel("x (m)",labelpad=-6)
    ax.set_ylabel("y (m)",labelpad=-6)
    ax.set_zlabel("z (m)")
    ax.set_xticks([-0.10,0,0.10])
    ax.set_yticks([-0.10,0,0.10])
    ax.set_zticks([ 0,0.10,0.20,0.30])
    ax.tick_params(axis = "x", pad=-5)
    ax.tick_params(axis = "y", pad=-5)
    ax.set_xlim3d(x_min,x_max)
    ax.set_ylim3d(y_min,y_max)
    ax.set_zlim3d(z_min,z_max)
    ax.set_box_aspect((xspan, yspan, zspan))
    ax.view_init(;elev, azim)
end

function position_plot(θs;row=2,nstage=2,zup=0.2,elev=16,azim=-50,y=0.0)
    n = length(θs)
    h = 0.066
    fig = plt.figure(figsize=(8,5))
    axs = [fig.add_subplot(row,ceil(Int,n/row),i,projection="3d") for i = 1:n]
    θstrings = [latexstring("\\theta=0")]
    for i = 1:n-1
        raco = 1//12*i/(n-1)
        raco_n = numerator(raco)
        raco_d = denominator(raco)
        if raco_n == 1
            push!(θstrings, latexstring("\\theta=\\pi/$raco_d"))
        else
            push!(θstrings, latexstring("\\theta=$raco_n\\pi/$raco_d"))
        end
    end
    alphabet = join('a':'z')
    for i = 1:n
        iθ = θs[i]
        R = RotY(iθ)
        reflinkn = links(nstage,h,R)
        ax = axs[i]
        bars,strings = bars_and_strings_segs_3D(reflinkn)
        ax.add_collection3d(strings)
        ax.add_collection3d(bars)
        ax.set_title("($(alphabet[i])) "*θstrings[i];y)
        set_ax!(ax,zup;elev,azim)
    end
    fig.tight_layout()
    fig
end

function inverse_analysis(θs)
    n = length(θs)
    nstage = 2
    h = 0.066
    linkn = links(nstage,h,RotY(0.0);k=3e1)
    # return linkn
    Y = Array(build_Y(linkn))

    # ωs = Vector{Vector{Float64}}()
    # Zs = Vector{Matrix{Float64}}()
    rls = Vector{Vector{Float64}}()
    # epes = Vector{Float64}()


    # show(θstrings)
    for i = 1:n
        iθ = θs[i]
        @show i,iθ
        R = RotY(iθ)
        reflinkn = links(nstage,h,R)
        refqi,_,_ = RB.get_initial(reflinkn)
        # actlinkn,lhs,rhs = RB.statics_equation_for_actuation(linkn,reflinkn,Y,gravity=false,scale=false)
        # λpt = zeros(actlinkn.num_of_cstr)
        # lengths = [s.state.length for s in actlinkn.strings]
        # u0 = [s.original_restlen for s in actlinkn.strings]
        # apt = lengths - u0
        # @show lhs*vcat(λpt,apt) ≈ rhs

        refλi,_,ai= RB.inverse_for_actuation(linkn,reflinkn,Y,gravity=false,scale=false)
        actlinkn = deepcopy(linkn)
        RB.actuate!(actlinkn,ai)
        rli = [s.state.restlen for s in actlinkn.strings]
        push!(rls,rli)
        # RB.reset_forces!(actlinkn)
        # RB.distribute_q_to_rbs!(actlinkn,refqi)
        # RB.update_strings_apply_forces!(actlinkn)
        # tensions_i = [s.state.tension for s in actlinkn.strings]
        # @show tensions_i
        #
        # epei = RB.elastic_potential_energy(actlinkn,refqi)
        # push!(epes,epei)


        # ωi,Zi = RB.undamped_eigen(actlinkn,refqi,refλi)
        # push!(ωs,ωi)
        # push!(Zs,Zi)
    end
    # fig.tight_layout()
    # linkn, ωs, Zs, rls, epes, fig
    linkn,rls
end


function build_EPE_from_actuation(tr,Y,index=:)
    @unpack st = tr
    ℓ = RB.build_ℓ(st)
    u0 = RB.get_original_restlen(tr)
    k = [s.k for s in st.strings]
    function inner_V(a)
        Δℓ = ℓ[index]-u0[index]-(Y*a)[index]
        ret = 1/2*transpose(Δℓ)*Diagonal(k[index])*Δℓ
    end
end

function EPE_from_actuation(tgstruct,Y,a)
    Vfunc = build_EPE_from_actuation(tgstruct,Y)
    Vfunc(a)
end


function inverse2energy(tgstruct,refstruct,Y)
    actstruct,lhs,rhs = RB.statics_equation_for_actuation(tgstruct,refstruct,Y,gravity=false,scale=false)
    λpt = zeros(actstruct.num_of_cstr)
    lengths = [s.state.length for s in actstruct.strings]
    u0 = [s.original_restlen for s in actstruct.strings]
    apt = lengths - u0
    x0 = vcat(λpt,apt)
    @info @show lhs*x0 ≈ rhs
    xp,nb = RB.get_solution_set(lhs,rhs)
    x(ξ) = x0 + nb*ξ
    # ξ0 = nb\(x0-xp)
    # apt = lengths - u0
    # 0 < u0 < lengths
    # -lengths < -u0 < 0
    # 0 < apt < lengths
    ξ_limit = -u0./nb[actstruct.num_of_cstr+1:end]
    ξ_min_index = argmin(abs.(ξ_limit))
    ξ_max = ξ_limit[ξ_min_index]
    @show ξ_max
    nξ = 100
    # ξ_range =
    epeξ = Vector{Float64}()
    rlξ = Vector{Vector{Float64}}()
    for (i,ξ) in enumerate(LinRange(0,ξ_max,nξ))
        # @show i,ξ
        xi = x(ξ)
        ai = xi[actstruct.num_of_cstr+1:end]
        RB.actuate!(actstruct,ai)
        epei = RB.elastic_potential_energy(actstruct)
        epei_a = EPE_from_actuation(actstruct,Y,ai)
        @assert epei_a ≈ epei
        push!(epeξ, epei)
        rli = u0 + Y*ai
        push!(rlξ,rli)
    end
    rlξ, epeξ
end

function inverse2energy_plot(θs,rl_indexes)
    n = length(θs)
    nstage = 2
    h = 0.066
    linkn = links(nstage,h,RotY(0.0);k=3e1)
    Y = Array(build_Y(linkn))
    for i = 1:n
        iθ = θs[i]
        @show i,iθ
        R = RotY(iθ)
        reflinkn = links(nstage,h,R)
        rlξ_raw, epeξ = inverse2energy(linkn,reflinkn,Y)
        rlξ = [[rlξ_raw[i][j] for i in 1:length(rlξ_raw)] for j in 1:length(rlξ_raw[1])]
        fig,ax = plt.subplots(1,1,figsize=(4.5,3))
        for group in rl_indexes[i]
            label = latexstring(join(["s_$i" for i in group],","))
            ax.plot(epeξ,rlξ[group[1]];label)
        end
        ax.set_xlim(-0.01,0.20)
        ax.set_ylim(0.00,0.10)
        ax.set_xlabel("Elastic Potential Energy (J)")
        ax.set_ylabel("Rest Length (m)")
        ax.legend(loc="upper right")
        # ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0) )
        ax.grid("on")
        fig.tight_layout()
        fig.savefig("link2_inv_pose$(i)_energy.png",dpi=300,bbox_inches="tight")
    end
    linkn
end


function inverse_with_energy(tgstruct,refstruct,Y,epe)
    actstruct,lhs,rhs = RB.statics_equation_for_actuation(tgstruct,refstruct,Y,gravity=false,scale=false)
    xp,nb = RB.get_solution_set(lhs,rhs)
    Vfunc = build_EPE_from_actuation(refstruct,Y)
    function f!(ret,ξ)
        x = xp + nb*ξ
        a = x[tgstruct.num_of_cstr+1:end]
        ret[1] = Vfunc(a) - epe
    end
    nl_result = nlsolve(f!,[0.0], autodiff = :forward)
    if !converged(nl_result)
        @error "Not converged!"
    end
    x_result = (xp + nb*nl_result.zero)
    λ_result = x_result[1:actstruct.num_of_cstr]
    a_result = x_result[actstruct.num_of_cstr+1:end]
    RB.actuation_check(actstruct,Y,a_result)
    u0 = [s.original_restlen for s in actstruct.strings]
    rl_result = u0 + Y*a_result
    λ_result,rl_result,a_result
end

function inverse_with_energy_plot(θs,rl_index)
    n = length(θs)
    nstage = 2
    h = 0.066
    linkn = links(nstage,h,RotY(0.0);k=3e1)
    Y = Array(build_Y(linkn))
    rls = [zeros(n) for i = 1:linkn.nstrings]
    for i = 1:n
        iθ = θs[i]
        @show i,iθ
        R = RotY(iθ)
        reflinkn = links(nstage,h,R)
        epe = 0.05
        _,rli,_ = inverse_with_energy(linkn,reflinkn,Y,epe)
        for j = 1:linkn.nstrings
            rls[j][i] = rli[j]
        end
    end
    fig, ax = plt.subplots(1,1,figsize=(4.5,3))
    for group in rl_index
        label = latexstring(join(["s_$i" for i in group],","))
        ax.plot(θs,rls[group[1]];label)
    end
    ax.set_ylim(0.0,0.08)
    ax.set_xlabel(L"\theta\ (\mathrm{rad})")
    ax.set_ylabel("Rest Length (m)")
    # ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.0))
    ax.legend(loc="lower left")
    ax.grid("on")
    fig.tight_layout()
    fig.savefig("link2_inv_by_energy.png",dpi=300,bbox_inches="tight")
end


function inverse_with_energy_3(tr,reftr,Y,epes)
    acttg,lhs,rhs,c = RB.statics_equation_for_actuation(tr.st,reftr.st,Y,gravity=false,scale=false)
    acttr = RB.Robot(acttg,deepcopy(tr.hub))
    xp,nb = RB.get_solution_set(lhs,rhs)
    Vfuncs = [build_EPE_from_actuation(acttr,Y,6(i-1)+1:6i) for i = 1:3]
    function f!(ret,ξ)
        x = xp + nb*ξ
        a = x[acttg.num_of_cstr+1:end]
        for i = 1:3
            ret[i] = Vfuncs[i](a) - epes[i]
        end
    end
    nl_result = nlsolve(f!,zeros(3), autodiff = :forward)
    if !converged(nl_result)
        @error "Not converged!"
    end
    x_result = (xp + nb*nl_result.zero)
    λ_result = x_result[1:acttg.num_of_cstr].*c
    a_result = x_result[acttg.num_of_cstr+1:end]
    RB.actuation_check(acttr,Y,a_result)
    refq,_ = RB.get_coords(reftr.st)
    RB.actuate!(acttr,a_result)
    @assert RB.check_static_equilibrium(acttg,refq,λ_result)
    u0 = RB.get_original_restlen(acttr)
    rl_result = u0 + Y*a_result
    λ_result,rl_result,a_result
end

function inverse_with_energy_3_plot(θs,rl_indexes)
    n = length(θs)
    nstage = 4
    h = 0.066
    linkn = links(nstage,h,RotY(0.0);k=3e1)
    Y = Array(build_Y(linkn))
    rls = [zeros(n) for i = 1:linkn.nstrings]
    for i = 1:n
        iθ = θs[i]
        @show i,iθ
        R = RotY(iθ)
        reflinkn = links(nstage,h,R)
        epes = [0.07,0.06,0.05]
        # epes = fill(0.05,3)
        _,rli,_ = inverse_with_energy_3(linkn,reflinkn,Y,epes)
        for j = 1:linkn.nstrings
            rls[j][i] = rli[j]
        end
    end
    fig, axs = plt.subplots(1,3,figsize=(11,4))
    alphabet = join('a':'z')
    for i = 1:3
        ax = axs[i]
        for group in rl_indexes[i]
            label = latexstring(join(["s_{$i}" for i in group],","))
            ax.plot(θs,rls[group[1]];label)
        end
        ax.set_ylim(0.02,0.07)
        ax.set_xlabel(L"\theta\ (\mathrm{rad})")
        ax.set_ylabel("Rest Length (m)")
        ax.legend(loc="lower left")
        ax.set_title("($(alphabet[i]))",y=-0.1)
        ax.grid("on")
    end
    fig.tight_layout()
    fig.savefig("link2_inv_by_energy_3.png",dpi=300,bbox_inches="tight")
end
