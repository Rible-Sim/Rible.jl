material_properties = Table(
    [
        (name = "ABS",              density = 1045u"kg/m^3",     modulus_elas = 3.1u"GPa"),
        # (name = "Acetate", density = 2.0u"GPa"),
        (name = "Acrylic",          density = 1185u"kg/m^3",     modulus_elas = 3.5u"GPa"),
        # (name = "Cartilage", density = 2.0u"GPa"),
        (name = "Diamond",          density = 3300u"kg/m^3",     modulus_elas = 1200u"GPa"),
        (name = "Graphite fibre",   density = 2300u"kg/m^3",     modulus_elas = 700u"GPa"),
        (name = "Ceramics",         density = 3600u"kg/m^3",     modulus_elas = 276u"GPa"),
        # (name = "Carbides", density = 2.0u"GPa"),
        # (name = "CFR epoxy", density = 2.0u"GPa"),
        (name = "Concrete",         density = 2400u"kg/m^3",     modulus_elas = 14u"GPa"),
        # (name = "Douglas fir", density = 2.0u"GPa"),
        (name = "Epoxy",            density = 1150u"kg/m^3",     modulus_elas = 8u"GPa"),
        # (name = "GFR polyester", density = 2.0u"GPa"),
        (name = "Glasses",          density = 4100u"kg/m^3",     modulus_elas = 70u"GPa"),
        # (name = "HDPE", density = 2.0u"GPa"),
        # (name = "LDPE", density = 2.0u"GPa"),
        (name = "Nylon 66",         density = 1150u"kg/m^3",     modulus_elas = 3.3u"GPa"),
        (name = "Oak",              density =  650u"kg/m^3",     modulus_elas =  12u"GPa"),
        (name = "Teak",             density =  630u"kg/m^3",     modulus_elas =  12u"GPa"),
        (name = "Polyester",        density = 1350u"kg/m^3",     modulus_elas = 2.4u"GPa"),
        (name = "PC",               density = 1200u"kg/m^3",     modulus_elas = 2.4u"GPa"),
        (name = "PP",               density = 900u"kg/m^3",      modulus_elas = 1.6u"GPa"),
        (name = "PVC",              density = 1700u"kg/m^3",     modulus_elas = 4.1u"GPa"),
        (name = "Aluminium",        density = 2710u"kg/m^3",     modulus_elas = 71u"GPa"),
        (name = "Al-Cu alloy",      density = 2800u"kg/m^3",     modulus_elas = 75u"GPa"),
        (name = "Al-Mg alloy",      density = 2725u"kg/m^3",     modulus_elas = 71u"GPa"),
        (name = "Alloy steel",      density = 7900u"kg/m^3",     modulus_elas = 210u"GPa"),
        (name = "Brass",            density = 8500u"kg/m^3",     modulus_elas = 104u"GPa"),
        (name = "Bronze",           density = 8800u"kg/m^3",     modulus_elas = 117u"GPa"),
        (name = "Cast iron",        density = 7150u"kg/m^3",     modulus_elas = 97u"GPa"),
        (name = "Copper",           density = 8950u"kg/m^3",     modulus_elas = 117u"GPa"),
        (name = "Gold",             density = 19300u"kg/m^3",    modulus_elas = 79u"GPa"),
        (name = "Iron",             density = 7850u"kg/m^3",     modulus_elas = 206u"GPa"),
        (name = "Lead",             density = 11370u"kg/m^3",    modulus_elas = 18u"GPa"),
        (name = "Ni-Alloy",         density = 9000u"kg/m^3",     modulus_elas = 200u"GPa"),
        (name = "Platinum",         density = 21040u"kg/m^3",    modulus_elas = 164u"GPa"),
        (name = "Silver",           density = 10530u"kg/m^3",    modulus_elas = 78u"GPa"),
        (name = "Stainless steel",  density = 7930u"kg/m^3",     modulus_elas = 200u"GPa"),
        (name = "Titanium",         density = 4540u"kg/m^3",     modulus_elas = 118u"GPa"),
    ]
)

function dynfuncs(bot;actuate=false,gravity=false,(Fˣ!)=(F,t)->nothing)
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
        Fˣ!(F,t)
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
