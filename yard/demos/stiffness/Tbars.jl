figdir::String = ""
if Sys.iswindows()
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\FrictionalContact\CMAME"
elseif Sys.isapple()
    figdir::String = raw"/Users/jacob/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/Papers/FrictionalContact/CMAME"
end
#-- end

include("deps.jl")
using AbbreviatedStackTraces #jl
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true #jl
ENV["JULIA_STACKTRACE_MINIMAL"] = true #jl
import Rible as RB
include("../../vis.jl")
includet("../../vis.jl") #jl

#-- T bars
include("../../../examples/robots/Tbars.jl")
includet("../../../examples/robots/Tbars.jl")#jl
tbbot = Tbars(;θ=π/4)
bot = tbbot
@myshow bot.structure.num_of_dof
bodies = RB.get_bodies(bot)
body1 = bodies[1]
dt = 1e-3
tspan = (0.0,5.0)
prob = RB.DynamicsProblem(bot,)
solver = RB.Zhong06()
RB.solve!(
    prob,
    RB.DynamicsSolver(
        solver
    );
    dt,tspan,ftol=1e-10,maxiters=50,verbose=false,exception=true,progress=false,
)

plot_traj!(bot;showarrows = false, showground=false)

#todo slider with no hook
GM.activate!();with_theme(theme_pub;
        resolution = (1tw,0.15tw),
        figure_padding = (0,0,fontsize/2,0),
        Axis3 = (
            # azimuth = -1.8701322643948965,
            # elevation = 0.6128666392948965,
            azimuth = -1.9234135143948998,
            elevation = 0.22103070179489467,
            perspectiveness = 0.3,
        )
    ) do 
    fig = Figure()
    gd0 = fig[1,1] = GridLayout(;tellheight=false)
    gd1 = fig[1,2] = GridLayout(;tellheight=false)
    # rowsize!(gd0,1,Fixed(0.15tw))
    # rowsize!(gd1,1,Fixed(0.15tw))
    bot0 = Tbars(;θ=0)
    plot_traj!(
        bot0,
        fig = gd0,
        AxisType = Axis3,
        showpoints = false,
        showlabels = false,
        showground = false,
        doslide = false,
        # showinfo = true,
        xlims = (-2.2,1.2),
        ylims = (-1.2,1.2),
        zlims = (0,0.06),
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[sgi])) ", font=:bold),
                "Initial"
            )
        end,
        sup! = (ax,tgob,sgi)-> begin
            # cables
            hidez(ax)
        end
    )
    bot1 = Tbars(;θ=π/4)
    plot_traj!(
        bot1,
        fig = gd1,
        AxisType = Axis3,
        showpoints = false,
        showlabels = false,
        showground = false,
        doslide = false,
        xlims = (-2.2,1.2),
        ylims = (-1.2,1.2),
        zlims = (0,0.06),
        titleformatfunc = (sgi,tt)-> begin
            rich(
                rich("($(alphabet[sgi+1])) ", font=:bold),
                "Movement"
            )
        end,
        sup! = (ax,tgob,sgi)-> begin
            # cables
            hidez(ax)
        end
    )
    savefig(fig,"Tbars")
    DataInspector(fig)
    fig
end

tbbot = Tbars(;θ=0.0)
bot = tbbot

RB.check_static_equilibrium_output_multipliers(bot.structure)

function nullspace_on_free(st)    
    (;sys_free_idx,bodyid2sys_full_coords,bodyid2sys_dof_idx) = st.connectivity.indexed
    q = RB.get_coords(bot.structure)
    Nin = RB.intrinsic_nullspace(st,q)[
        sys_free_idx,
        reduce(vcat,bodyid2sys_dof_idx[2:end])
    ]
    I3 = I(3)
    O3 = zero(I3)
    o3 = O3[:,1]
    qbar = q[bodyid2sys_full_coords[4]]
    ri = @view qbar[1:3]
    u = @view qbar[4:6]
    v,w = RB.NCF.HouseholderOrthogonalization(u)
    @myshow v,w
    x,_,_ = ri
    Nslider1 = [
        o3;
        o3;;
    ] |> sparse
    Nslider2 = [
         0;
        -x;
         0;
        o3;;
    ] |> sparse
    Nbar = [
        o3;
        0;
        1;;
    ] |> sparse
    Nex = vcat(Nslider1,Nslider2,Nbar)
    Nin*Nex
end

q = RB.get_coords(bot.structure)
q̌ = RB.get_free_coords(bot.structure)
Ǎ = RB.make_cstr_jacobian(bot.structure)(q)
Ň_ = RB.nullspace(Ǎ)
Ň = RB.modified_gram_schmidt(Ň_)

Ň = nullspace_on_free(bot.structure)

# done construct null space 
#note only work in θ = 0
rank(Ň)

Ǎ*Ň |> norm

Q̃ = RB.build_Q̃(bot.structure)
L̂ = RB.build_L̂(bot.structure)

# Left hand side
Q̃L̂ = Q̃*L̂

Bᵀ = -Q̃L̂
ℬᵀ = transpose(Ň)*Bᵀ
ℬ = transpose(ℬᵀ)
S,D = RB.static_kinematic_determine(ℬᵀ;atol=1e-14)
S 
D
nk = size(D,2)
ns = size(S,2)

k = RB.get_cables_stiffness(bot.structure)

l = RB.get_cables_len(bot.structure)

struct𝒦 = [
    begin
        s = S[:,i]        
        Ǩm = RB.build_material_stiffness_matrix!(bot.structure,q,s)
        𝒦m = transpose(Ň)*Ǩm*Ň 
        # s = S\f
        # @show s
        λ = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*s
        @show λ
        Ǩa = - RB.cstr_forces_jacobian(bot.structure,q,λ)

        Ǩg = RB.build_geometric_stiffness_matrix!(bot.structure,q,s)

        𝒦g = transpose(Ň)*(Ǩg)*Ň

        𝒦a = transpose(Ň)*(Ǩa)*Ň

        𝒦p = 𝒦g .+ 𝒦a
        @eponymtuple(𝒦m, 𝒦g, 𝒦a, 𝒦p,)
    end
    for i = 1:ns
] |> StructArray

mat𝒦ms = reduce(hcat,struct𝒦.𝒦m)
mat𝒦gs = reduce(hcat,struct𝒦.𝒦g)
mat𝒦as = reduce(hcat,struct𝒦.𝒦a)
mat𝒦ps = reduce(hcat,struct𝒦.𝒦p)

GM.activate!();with_theme(RB.theme_pub;
        resolution = (0.95tw,0.18tw),
        figure_padding = (0,0,0,0),
        Axis3 = (
            azimuth = -π/2-1e-10,
            elevation = π/2,
        )
    ) do 
    maxS = maximum(abs.(S))
    rtol = 1e-10
    Sbool = S.> maxS*rtol
    S[.!Sbool] .= 0.0
    fig = Figure()
    gd = fig[1,1:4] = GridLayout(;tellheight=false)
    bot0 = Tbars(;θ=0)
    plot_traj!(
        bot0,
        fig = gd,
        AxisType=Axis3,
        gridsize = (1,4),
        showpoints = false,
        showlabels = false,
        showground = false,
        doslide = false,
        showcables = false,
        xlims = (-2.2,1.2),
        ylims = (-1.2,1.2),
        rowgap=0,
        titleformatfunc = (sgi,tt)-> begin
            rich(
                    rich("($(alphabet[sgi])) ", font=:bold),
                    [
                        "Self-stress State 1",
                        "Self-stress State 2",
                        "Self-stress State 3",
                        "Self-stress State 4"
                    ][sgi]
                )
        end,
        sup! = (ax,tgob,sgi)-> begin
            # cables
            RB.hidexyz(ax)
            @myshow Sbool[:,sgi]
            linesegs_cables = @lift begin
                get_linesegs_cables($tgob;)[Sbool[:,sgi]]
            end
            linesegments!(ax, 
                linesegs_cables, 
                color = :red, 
                # linewidth = cablewidth
                )
            rcs_by_cables = @lift begin
                (;tensioned) = $tgob.connectivity
                num_of_dim = RB.get_num_of_dims($tgob)
                T = RB.get_numbertype($tgob)
                ret = Vector{MVector{num_of_dim,T}}()
                mapreduce(
                    (scnt)->
                    [(
                        scnt.hen.rbsig.state.loci_states[scnt.hen.pid].position.+
                        scnt.egg.rbsig.state.loci_states[scnt.egg.pid].position
                    )./2],
                    vcat,
                    tensioned.connected
                    ;init=ret
                )
            end
            # @show rcs_by_cables
            Stext = [
                    @sprintf "%4.2f"  S[i,sgi] 
                    for i in axes(S,1)
                    if Sbool[i,sgi]
                ]
            @myshow Stext
            # scatter!(
            #     ax,
            #     rcs_by_cables[][Sbool[:,sgi]],
            #     marker = :rect, 
            #     markersize = 12 |> pt2px, 
            #     color = :white
            # )
            text!(
                ax,
                Stext,
                position = rcs_by_cables[][Sbool[:,sgi]],
                fontsize = 5 |> pt2px,
                color = :red,
                align = (:left, :top),
                offset = (fontsize/4, 0)
            )
        end
    )
    savefig(fig,"Tbars_states")
    DataInspector(fig)
    fig
end

ᾱs = [
   [1,0,0,0],
   [0,1,0,0],
   [0,0,1,0],
   [0,0,0,1],
]

@myshow mat𝒦ps
σs = 0:0.1:10
Vals =  [
    begin
        [   
            begin
                𝒦 = mat𝒦ms*ᾱ + σ*mat𝒦ps*ᾱ
                𝒦[1]
            end
            for ᾱ in ᾱs
        ]
    end
    for σ in σs
] |> VectorOfArray

GM.activate!();with_theme(RB.theme_pub;
        resolution = (0.45tw,0.2tw),
        figure_padding = (0,fontsize,0,fontsize),
    ) do 
    fig = Figure()
    ax1 = Axis(fig[1,1],
        xlabel = L"\sigma",
        ylabel = L"\rho_{(\mathrm{1})}"
    )
    lines!(ax1,σs,Vals[1,:],label=("Self-stress State 1 and 2"))
    lines!(ax1,σs,Vals[3,:],label=("Self-stress State 3"))
    lines!(ax1,σs,Vals[4,:],label=("Self-stress State 4"))
    xlims!(ax1,0,10)

    # ylims!(ax1,-0,6)
    # scatter!(
    #     ax1,
    #     [σ_max,σ_zero],
    #     [ρ_max,ρ_zero]
    # )
    
    # ax2 = Axis(fig[1,2],
    #     xlabel = L"\sigma",
    #     ylabel = L"\rho"
    # )
    # Legend(
    #     fig[1,3],
    #     ax2
    # )
    # xlims!(ax2,0,1700)
    # ylims!(ax2,-20,400)
    Legend(fig[1,2],
        ax1
    )
    savefig(fig,"Tbars_curve")
    fig
end

#-- end Tbars 