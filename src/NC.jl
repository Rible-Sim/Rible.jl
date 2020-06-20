using .NaturalCoordinates
struct NaturalCoordinatesCache{ArrayT,MT,fT}
    CG::ArrayT
    Cp::Vector{ArrayT}
    M::MT
    funcs::fT
end

function NaturalCoordinatesCache(prop,L::Real)
    bps = NaturalCoordinates.BasicPoints2P(L)
    cf = NaturalCoordinates.CoordinateFunctions(bps)
    @unpack mass,inertia,CoM,naps,aps = prop
    M = NaturalCoordinates.make_M(cf,mass,inertia,CoM)
    @unpack C,c = cf
    CG = C(c(CoM))
    Cp = [SMatrix{2,4}(C(c(aps[i])))
            for i in 1:naps]
    NaturalCoordinatesCache(CG,Cp,M,cf)
end
