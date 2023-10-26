
function actuation_interpolation(target_t,initial_a, target_a,
                                Alg=Quadratic(Flat(OnGrid())))
    itps = [interpolate([a0,af],
            BSpline(Alg)) for (a0,af) in zip(initial_a,target_a)]
    function itp(t)
        [itp(t) for itp in itps]
    end
    function inner_aitp(t)
        if t < target_t
            itp(1+t/target_t)
        else
            target_a
        end
    end
end

function simulate_linearactuating(ndof=6;k,c,target_t=10.0,tend=20.0,Alg=Linear())
    manipulator = man_ndof(ndof,k=k,c=c)
    Y = build_Y(manipulator)

    refman = man_ndof(ndof,θ=-π/5,k=k) # reference
    _,_,initial_a = RB.inverse(manipulator,deepcopy(manipulator),Y)
    _,_,target_a = RB.inverse(manipulator,deepcopy(refman),Y)

    aitp = actuation_interpolation(target_t,initial_a,target_a,Alg)

    function linearactuate(tr,aitp)
        @unpack st = tr
        M = RB.build_massmatrix(st)
        Φ = RB.build_Φ(st)
        A = RB.build_A(st)
        # Q̃ = RB.build_Q̃(st)
        function F!(F,q,q̇,t)
            RB.actuate!(tr,aitp(t))
            RB.reset_forces!(st)
            RB.distribute_q_to_rbs!(st,q,q̇)
            RB.update_cables_apply_forces!(st)
            RB.assemble_forces!(F,st)
        end

        M,Φ,A,F!,nothing
    end

    dt = 0.01 # Same dt used for PID AND Dynamics solver
    prob = RB.DyProblem(linearactuate(manipulator,aitp),manipulator,(0.0,tend))

    RB.solve!(manipulator,prob,RB.Zhong06(),dt=dt,ftol=1e-14)
    actuation_trajs = aitp.(manipulator.traj.ts)
    for actuation_traj in actuation_trajs
        foreach((x,a)->push!(x.traj.us,a),manipulator.hub.actuators,actuation_traj)
    end
    manipulator,aitp
end

function plot_restlen(man,aitp_linear,aitp_quadra)
    ts = man.traj.ts
    u0 = RB.get_original_restlen(man)
    Y = build_Y(man)
    fig, ax = plt.subplots(1,1,figsize=(Single_width,2.5))
    rl_linear = VectorOfArray([Y*a + u0 for a in aitp_linear.(ts)])
    rl_quadra = VectorOfArray([Y*a + u0 for a in aitp_quadra.(ts)])
    ax.plot(ts,rl_linear[1,:], ls=linestyles.xs[1], label="Linear")
    # ax.plot(ts,rl_quadra[1,:], ls=linestyles.xs[2], label="Quadratic")
    ax.set_ylim(0.28,0.34)
    ax.set_xlim(0,25)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Rest Length (m)")
    ax.grid("on")
    # ax.legend()
    fig.tight_layout()
    fig
end


function pyplot2(man_linear,tstops)
    @unpack st, traj = man_linear
    steps = [findfirst((x)->x>i,traj.ts) for i = tstops]

    fig,axs_raw = plt.subplots(2,3,figsize=(Full_width,4.5))
    axs = permutedims(axs_raw,[2,1])
    for (i,step) in enumerate(steps)
        RB.distribute_q_to_rbs!(st,traj.qs[step])
        bars_segs,cables_segs = bars_and_cables_segs(st)
        ax = axs[i]
        ax.add_collection(bars_segs)
        ax.add_collection(cables_segs)
        ax.set_ylim(-1.00,0.20)
        ax.set_xlim(-0.25,1.35)
        ax.set_aspect("equal")
        ax.set_title("($(alphabet[i])) "*
                    latexstring("t=$(tstops[i])\\mathrm{s}"),
                        y = -0.35)
        ax.grid(true)
        # ax.set_xticks([0,0.5,1])
        # ax.set_yticks([-1.5,-1.0,-0.5,0])
        if i ∈[0]
            ax.set_xticklabels([])
            ax.xaxis.label.set_visible(false)
        end
        if i ∈[2,3,5,6]
            ax.set_yticklabels([])
            ax.yaxis.label.set_visible(false)
        end
    end
    fig.tight_layout()
    fig
end
