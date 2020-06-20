rb1 = tail.rigidbodies[1]
rb2 = tail.rigidbodies[2]
L = 0.02
function make_cache(prop,L)
    @unpack mass,inertia,CoM = prop
    bps = TR.NaturalCoordinates.BasicPoints2P(L)
    cf = TR.NaturalCoordinates.CoordinateFunctions(bps)
    @unpack mass,inertia,CoM = prop
    M = TR.NaturalCoordinates.make_M(cf,mass,inertia,CoM)
    cf,M
end

@code_warntype TR.NaturalCoordinates.BasicPoints2P(L)
cf,M = make_cache(rb1.prop,L)
@code_warntype make_cache(rb1.prop,L)

cfM == rb1.state.cache.M
