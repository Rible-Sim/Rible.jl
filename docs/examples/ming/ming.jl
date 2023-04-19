using LinearAlgebra
using SparseArrays
using StaticArrays
using Rotations
using Parameters
using Makie
import PyPlot as plt
plt.pygui(true)
using LaTeXStrings
#using ForwardDiff
# using DynamicPolynomials
# using Printf
using RecursiveArrayTools
using XLSX
using DataFrames
using Loess
using AxisArrays
using EasyFit
using CoordinateTransformations
# using HomotopyContinuation
using TypeSortedCollections
using Interpolations
using GeometryBasics
using FileIO
using NLsolve
using Revise
import TensegrityRobots as TR
cd("examples/ming")
includet("links_define.jl")
includet("links_plotting.jl")
includet("../analysis.jl")
includet("../plot_helpers.jl")
set_pyplot()
set_theme!(font = "Times New Roman")
function fit_data()
    sma_data = XLSX.readxlsx("SMA_data.xlsx")["Sheet1"]
    sma_df = DataFrame(
        XLSX.readtable(
        "SMA_data.xlsx", "Sheet1",infer_eltypes=true)...)
    sorted_deforms = sort(unique(sma_df.形变量))
    exclude_deform_indices = findall((x)->x∈[14,15,16,17,25,26,35,39],sma_df.形变量)
    exclude_temp_indices = findall((x)->x>100,sma_df.实时温度)
    exclude_indices = union(exclude_deform_indices,exclude_temp_indices)
    forces = sma_df[Not(exclude_indices),:回复力]
    temperatures = sma_df[Not(exclude_indices),:实时温度]
    deformations = sma_df[Not(exclude_indices),:形变量]
    include_defrom_indices = findall((x)->!(x∈[14,15,16,17,25,26,35,39]),sorted_deforms)
    @show sorted_deforms[include_defrom_indices]
    # fig = Figure()
    # ax = fig[1,1] = LScene(fig, scenekw = (camera = cam3d!, raw = false))
    # scatter!(ax,temperatures,deformations,forces,axis=false, markersize = 300)

    fig = plt.figure(figsize=(Full_width,8))
    ax = fig.add_subplot(1,1,1,projection="3d")
    ax.scatter3D(temperatures,float.(deformations),forces)

    np = 100
    ts = collect(range(25,90; length = np))
    sets_deform = sort(unique(deformations))
    axis_forces = AxisArray(zeros(np, length(sets_deform)); temp = ts, deform = sets_deform)
    for (i,d) in enumerate(sets_deform)
        d_indices = findall((x)->x==d,deformations)
        d_temps = temperatures[d_indices]
        d_forces = forces[d_indices]
        loess_model = loess(d_temps,d_forces,span=0.3,degree=3)
        fs = predict(loess_model, ts)
        axis_forces[:,i] = fs
        # lines!(ax,ts,fill(float(d),np),fs)
        ax.plot3D(ts,fill(float(d),np),fs,color="black")
    end
    ds = range(extrema(sets_deform)...; length = np)
    ks = Vector{Float64}()
    F0s = Vector{Float64}()
    for (i,t) in enumerate(ts)
        t_forces = axis_forces[atvalue(t)]
        # loess_model = loess(sets_deform,t_forces,span=1.0,degree=2)
        # fs = predict(loess_model, ds)
        fit = fitlinear(sets_deform,t_forces)
        push!(ks,fit.a)
        push!(F0s,fit.b)
        # lines!(ax,fill(float(t),length(fit.x)),fit.x,fit.y)
        ax.plot3D(fill(float(t),length(fit.x)),fit.x,fit.y,color="orange")
    end
    # xlims!(ax.scene,100,20)
    # ylims!(ax.scene,50,0)
    # scale!(ax.scene, -1,-1,1)
    # fig

    k_linear_indices = findall((x)->x>30.5&&x<61,ts)
    k_fit = fitlinear(ts[k_linear_indices],ks[k_linear_indices])
    fig_k,ax_k = plt.subplots(1,1,figsize=(Full_width,5))
    ax_k.plot(ts,ks,ls="",marker="o",label="Data")
    ax_k.plot(k_fit.x,k_fit.y,label="Linear fit")
    ax_k.set_xlim(20,100)
    ax_k.set_xlabel("Temperature (T)")
    ax_k.set_ylabel("Stiffness (N/mm)")
    ax_k.grid("on")
    ax_k.legend()
    fig_k.tight_layout()
    fig_k.savefig("Temperature-Stiffness.png",dpi=300,bbox_inches="tight")

    F0_loess_model = loess(ts,F0s,span=0.2,degree=3)
    F0s_loess = predict(F0_loess_model, ts)
    fig_F0,ax_F0 = plt.subplots(1,1,figsize=(Full_width,5))
    ax_F0.plot(ts,F0s,ls="",marker="o",label="Data")
    ax_F0.plot(ts,F0s_loess,label="LOESS")
    ax_F0.set_xlim(20,100)
    ax_F0.set_xlabel("Temperature (T)")
    ax_F0.set_ylabel("Force (N)")
    ax_F0.grid("on")
    ax_F0.legend()
    fig_F0.tight_layout()
    fig_F0.savefig("Temperature-Force.png",dpi=300,bbox_inches="tight")

    F0_loess_model,k_fit,fig,ax
end
ENV["KMP_DUPLICATE_LIB_OK"]="TRUE"
F0_loess_model,k_fit,fig,ax = fit_data()
ax.set_xlim(20,100)
ax.set_ylim(0,50)
ax.set_zlim(0,30)
ax.set_box_aspect((80, 50, 30))
ax.view_init(;elev=27, azim=-119)
ax.set_xlabel("Temperature (T)")
ax.set_ylabel("Deformation (m)")
ax.set_zlabel("Force (N)")
fig.tight_layout()
fig.savefig("fitting.png",dpi=300,bbox_inches="tight")

function make_heating_law(k_fit,F0_model)
    function heating_law(T)
        F0 = predict(F0_model,float(T))
        k = (k_fit.a*T+k_fit.b)*1000
        if k < 0
            k = zero(k)
        end
        F0,k
    end
end

heating_law = make_heating_law(k_fit,F0_loess_model)
heating_law(50)

pyramid = read_rb_mesh()

spine3 = spine(5,0.1,RotZ(0.0),heating_law;c=10000.0)

fig,ax = plotstructure(spine3.tg,pyramid)
fig
function simulate_linkn(linkn;dt=0.01,tend=10.0,target_temps=[40,40,40],verbose=false)

    function dynfuncs(tr)
        @unpack tg = tr
        M = TR.build_massmatrix(tg)
        Φ = TR.build_Φ(tg)
        A = TR.build_A(tg)
        Q̃ = TR.build_Q̃(tg)
        td = 2.0
        T0 = 40.0
        alg = BSpline(Quadratic(Flat(OnGrid())))
        itp1 = interpolate([T0,target_temps[1]],alg)
        itp2 = interpolate([T0,target_temps[2]],alg)
        itp3 = interpolate([T0,target_temps[3]],alg)
        function F!(F,q,q̇,t)
            if t < td
                TR.heat!(tr.hub.heaters[1],tg,itp1(1+t/td))
                TR.heat!(tr.hub.heaters[2],tg,itp2(1+t/td))
                TR.heat!(tr.hub.heaters[3],tg,itp3(1+t/td))
            end
            TR.reset_forces!(tg)
            TR.distribute_q_to_rbs!(tg,q,q̇)
            TR.update_strings_apply_forces!(tg)
            TR.update_SMA_strings_apply_forces!(tg)
            TR.apply_gravity!(tg)
            # F .= Q̃*TR.fvector(tg)
            # rb3 = tg.rigidbodies[end]
            # f = [5.0,0.0,0.0]
            # rb3.state.Faps[13] .+= f
            TR.assemble_forces!(F,tg)
        end
        M,Φ,A,F!,nothing
    end

    reflinkn = deepcopy(linkn)
    q0,q̇0,λ0 = TR.get_initial(reflinkn.tg)

    prob = TR.DyProblem(dynfuncs(linkn),q0,q̇0,λ0,(0.0,tend))
    TR.solve!(linkn,prob,TR.Zhong06();dt,ftol=1e-14,verbose)

    linkn,reflinkn
end

_,_ = simulate_linkn(spine3,target_temps=[50,40,40])
_,_ = simulate_linkn(spine3,target_temps=[50,35,35])
_,_ = simulate_linkn(spine3,target_temp=45)
_,_ = simulate_linkn(spine3,target_temp=50)
_,_ = simulate_linkn(spine3,target_temp=55)
_,_ = simulate_linkn(spine3,target_temp=60)
plotstructure(spine3,pyramid)
plotstructure(spine3,pyramid;record_name="relax.mp4")
q_static = copy(spine3.traj.qs[end])
TR.set_new_initial!(spine3,q_static)
spine3_copy = deepcopy(spine3)
fig = plotstructure([spine3,spine3_copy],1,2,pyramid)
save("spines.png", fig)

TR.heat!(spine3.hub.heaters[1],spine3.tg,60.0)
_,_ = simulate_linkn(spine3)
plotstructure(spine3,pyramid)
plotstructure(spine3,pyramid;record_name="90.mp4")

function get_forces(tr,rbid)
    @unpack tg = tr
    rb = tg.rigidbodies[rbid]
    norm.(rb.state.Faps)
end

get_forces(spine3,2)

rp13 = get_trajectory(spine3,3,13)
x13 = rp13[1,:]
ts = spine3.traj.ts
fig,axs = plt.subplots(2,1,figsize=(Full_width,6))
axs[1].plot(ts,x13,ls=linestyles.xs[1],marker = markers.xs[1],label="a")
axs[1].legend(loc="upper left", bbox_to_anchor=(1.05, 1.04))
axs[2].plot(rand(10))
fig.tight_layout()
fig
