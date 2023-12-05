using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using Parameters
using Makie
AbstractPlotting.__init__()
plot(rand(10))


function eigen_freq(θs = [0.0], epes = [0.05]; mode_index = nothing)
    R = RotY(0.0)
    nstage = 4
    h = 0.066
    linkn = links(nstage,h,R;k=3e1)
    Y = build_Y(linkn)

    ncase = length(epes)
    default_cycler = plt.matplotlib.rcParams["axes.prop_cycle"]
    default_colors = default_cycler.by_key()["color"]
    freqs = VectorOfArray(Vector{Vector{Float64}}())
    for (i,(θ,epe)) in enumerate(zip(θs,epes))
        reflinkn = links(nstage,h,RotY(θ);k=3e1)
        q,_ = RB.get_coords(reflinkn)
        λ,_,a = inverse_with_energy_3(linkn,reflinkn,Y,fill(epe,3))
        # @show a[1:3:end]
        RB.actuate!(linkn,a)
        ω,Z = RB.undamped_eigen!(linkn,q,λ)
        tension = [s.state.tension for s in linkn.strings]
        @show tension[1:3:end]
        push!(freqs,ω)
    end
    nfirstmodes = linkn.num_of_dof

    # mode
    # fig,ax = plt.subplots(1,1,figsize=(5,4))
    # ω = freqs[:,1]
    # ax.plot(1:nfirstmodes,ω,marker="o",linestyle="none")
    # ax.set_ylim(0,90)
    # ax.set_xticks(1:nfirstmodes)
    # ax.set_xlabel("Mode")
    # ax.set_ylabel("Frequency (Hz)")

    # V_E
    # fig,ax = plt.subplots(1,1,figsize=(6.5,5))
    # ax.set_prop_cycle(color=vcat(default_colors,default_colors),
    #                  marker=vcat(fill("o",10),fill("v",10))
    #                  )
    # for group = mode_index
    #     label = "Mode "*join(["$m" for m in group],",")
    #     ax.plot(epes,freqs[group[1],:],fillstyle="none";label)
    #     ax.set_ylabel("Frequency (Hz)")
    #     ax.set_xlabel(L"V_E\ (\mathrm{J})")
    #     ax.set_xticks([0.05,0.1,0.15,0.2])
    #     # ax.set_xlim(1,20)
    #     # ax.set_ylim(0,100)
    # end
    # ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0))

    # θ
    fig,ax = plt.subplots(1,1,figsize=(6.5,5))
    ax.set_prop_cycle(color=vcat(default_colors,default_colors),
                     marker=vcat(fill("o",10),fill("v",10))
                     )
    for i = 1:nfirstmodes
        ax.plot(θs,freqs[i,:],fillstyle="none",label="Mode $i")
        ax.set_ylabel("Frequency (Hz)")
        ax.set_xlabel(L"\theta\ (\mathrm{rad})")
        ax.set_xticks(LinRange(0.0,0.25,6))
        # ax.set_xlim(1,20)
        # ax.set_ylim(0,100)
    end
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.0))


    ax.grid("on")
    fig.tight_layout()
    linkn,fig
end
linkn,fig = eigen_freq()
fig.savefig("link4_freq_zero.png",dpi=300,bbox_inches="tight")

mode_index = [[1,2],[3],[4],[5,6],[7],[8],[9],[10],[11],[12],[13],[14,15],[16],[17,18]]
mode_index = [[i] for i = 1:18]
linkn,fig = eigen_freq(LinRange(0.0,0.0,9),LinRange(0.05,0.20,9);mode_index)
fig.savefig("link4_freq_energy.png",dpi=300,bbox_inches="tight")

linkn,fig = eigen_freq(LinRange(0.0,π/12,9),LinRange(0.05,0.05,9))
fig.savefig("link4_freq_theta.png",dpi=300,bbox_inches="tight")


function plot_mode(epe=0.05)
    R = RotY(0.0)
    nstage = 4
    h = 0.066
    linkn = links(nstage,h,R;k=3e1)
    Y = build_Y(linkn)

    reflinkn = deepcopy(linkn)
    q,_ = RB.get_coords(reflinkn)
    λ,_,a = inverse_with_energy_3(linkn,reflinkn,Y,fill(epe,3))
    RB.actuate!(linkn,a)
    ω,Z = RB.undamped_eigen!(linkn,q,λ)

    fig = plt.figure(figsize=(12,9))
    modeindex = 1:6
    q0,_ = RB.get_coords(linkn)
    mode_tgstruct = deepcopy(linkn)
    alphabet = join('a':'z')
    for (i,figlabel) = zip(modeindex,alphabet)
        ax = fig.add_subplot(2,Int(ceil(length(modeindex)/2)),i,projection="3d")
        qmode = q0 + 0.01Z[:,i]
        RB.distribute_q_to_rbs!(mode_tgstruct,qmode)
        bars,strings = bars_and_strings_segs_3D(mode_tgstruct)
        refbars,refstrings = bars_and_strings_segs_3D(linkn;ref=true)
        ax.add_collection3d(refstrings)
        ax.add_collection3d(refbars)
        ax.add_collection3d(strings)
        ax.add_collection3d(bars)
        ax.set_title("($figlabel) Mode $i",y=-0.1)
        set_ax!(ax,0.35)
    end
    fig.tight_layout()
    fig
end
fig = plot_mode()
fig.savefig("link4_mode.png",dpi=300,bbox_inches="tight")
