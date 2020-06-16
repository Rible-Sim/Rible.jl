
function inertia2Z(inertia,Lij)
    z = 1/Lij^2*inertia
end

function form_mass_matrix(m,CoM,Lij,z)
    M = zeros(4,4)
    a = CoM./Lij
    M[1,1] = m - 2m*a[1] + z
    M[1,2] = 0.0
    M[1,3] = m*a[1] - z
    M[1,4] = -m*a[2]
    M[2,2] = m - 2m*a[1] + z
    M[2,3] = m*a[2]
    M[2,4] = m*a[1] - z
    M[3,3] = z
    M[3,4] = 0.0
    M[4,4] = z
    SymM = Matrix(Symmetric(M))
end

function C(c)
    ret = Matrix{eltype(c)}(undef,2,4)
    ret[1,1] = 1-c[1]
    ret[1,2] = c[2]
    ret[1,3] = c[1]
    ret[1,4] = -c[2]
    ret[2,1] = -c[2]
    ret[2,2] = 1-c[1]
    ret[2,3] = c[2]
    ret[2,4] = c[1]
    ret
end

function NCaux(prop,ri,rj)
    @unpack mass,CoM,inertia,anchorpoints = prop
    Lij = norm(ri-rj)
    z = inertia2Z(inertia,Lij)
    M = form_mass_matrix(mass,CoM,Lij,z)
    X̄⁻¹ = 1/Lij*Matrix(1.0I,2,2)
    function c(r̄)
        ret = X̄⁻¹*r̄
    end
    CG = C(c(CoM))
    Q = zeros(4)
    function Φ(q)
        xi,yi,xj,yj = q
        (xj-xi)^2 + (yj-yi)^2 - Lij^2
    end
    function Φq(q)
        xi,yi,xj,yj = q
        ret = similar(q)
        ret[1] =  2(xi-xj)
        ret[2] =  2(yi-yj)
        ret[3] =  2(xj-xi)
        ret[4] =  2(yj-yi)
        ret
    end
    nap = prop.number_aps
    Cp = [SMatrix{2,4}(C(c(anchorpoints[i])))
            for i in 1:nap]
    aux = NaturalCoordinatesAuxiliaries2D(
    SMatrix{4,4}(M),
    SMatrix{2,4}(CG),
    Cp,
    MVector{4}(Q),
    Lij,c,Φ,Φq
    )
end
