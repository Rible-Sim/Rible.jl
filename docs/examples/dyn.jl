function dynfuncs(bot;actuate=false,gravity=false)
    (;tg) = bot
    function F!(F,q,q̇,t)
        if actuate
            TR.actuate!(bot,[t])
        end
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_tensiles!(tg)
        if gravity
            TR.apply_gravity!(tg)
        end
        F .= TR.generate_forces!(tg)
    end
    function Jac_F!(∂F∂q̌,∂F∂q̌̇,q,q̇,t)
        ∂F∂q̌ .= 0
        ∂F∂q̌̇ .= 0
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.update_tensiles!(tg)
        TR.build_∂Q̌∂q̌!(∂F∂q̌,tg)
        TR.build_∂Q̌∂q̌̇!(∂F∂q̌̇,tg)
    end
    @eponymtuple(F!,Jac_F!)
end

function make_pres_actor(μ0,μ1,start,stop)
    nμ = length(μ0)

    function itp(t)
        scaled_itps = extrapolate(
            Interpolations.scale(
                interpolate(
                    hcat(μ0,μ1),
                    (NoInterp(),BSpline(Linear()))
                    # (NoInterp(),BSpline(Quadratic(Flat(OnGrid()))))
                ),
                1:nμ, start:stop-start:stop
            ),
            (Throw(),Flat())
        )
        [scaled_itps(j,t) for j in 1:nμ]
    end

    TR.PrescribedActuator(
        1,
        TR.ManualActuator(1,collect(1:nμ),zeros(nμ),TR.Uncoupled()),
        itp
    )
end
