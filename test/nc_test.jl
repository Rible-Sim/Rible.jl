rb1 = tail.rigidbodies[1]
rb2 = tail.rigidbodies[2]
L = 0.02
function make_cache(prop,L)
    @unpack mass,inertia,r̄g = prop
    lncs = TR.NaturalCoordinates.LocalNaturalCoordinates2P(L)
    cf = TR.NaturalCoordinates.CoordinateFunctions(lncs)
    @unpack mass,inertia,r̄g = prop
    M = TR.NaturalCoordinates.make_M(cf,mass,inertia,r̄g)
    cf,M
end

@code_warntype TR.NaturalCoordinates.LocalNaturalCoordinates2P(L)
cf,M = make_cache(rb1.prop,L)
@code_warntype make_cache(rb1.prop,L)

cfM == rb1.state.cache.M
