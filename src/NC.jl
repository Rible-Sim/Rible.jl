module NC
    using LinearAlgebra
    const I3 = Matrix(1.0I,3,3)
    # Is Float64 Okay?
    const r̄i = [0.0,0.0,0.0]
    const r̄j = [1.0,0.0,0.0]
    const ū =  [0.0,1.0,0.0]
    const v̄ =  [0.0,0.0,1.0]
    # Local Coordinates
    function X̄(r̄p)
        ret = Matrix{eltype(r̄p)}(undef,3,3)
        ret[:,1] .= r̄p .- r̄i
        ret[:,2] .= ū
        ret[:,3] .= v̄
        ret
    end

    const invX̄ = inv(X̄(r̄j))

    c(r̄p) = invX̄*(r̄p-r̄i)
    function C(c)
        ret = Matrix{eltype(c)}(undef,3,3*4)
        ret[:,1:3] = (1-c[1])*I3
        ret[:,4:6] = c[1]*I3
        ret[:,7:9] = c[2]*I3
        ret[:,10:12] = c[3]*I3
        ret
    end


    function inertia2Ji(inertia)
        Ji = similar(inertia)
        Ji .= -copy(inertia)
        Ji[1,1] = (inertia[2,2] + inertia[3,3] - inertia[1,1])/2
        Ji[2,2] = (inertia[1,1] + inertia[3,3] - inertia[2,2])/2
        Ji[3,3] = (inertia[1,1] + inertia[2,2] - inertia[3,3])/2
        return Ji
    end
    # Mass matrix
    # only once
    function form_mass_matrix(m,r̄G,Ji)
        a = invX̄*(r̄G-r̄i)
        z = invX̄*Ji*transpose(invX̄)
        M = vcat(
            hcat((m-2m*a[1]+z[1,1])*I3, (m*a[1]-z[1,1])*I3, (m*a[2]-z[1,2])*I3, (m*a[3]-z[1,3])*I3),
            hcat(   (m*a[1]-z[1,1])*I3, z[1,1]*I3, z[1,2]*I3, z[1,3]*I3),
            hcat(   (m*a[2]-z[2,1])*I3, z[2,1]*I3, z[2,2]*I3, z[2,3]*I3),
            hcat(   (m*a[3]-z[3,1])*I3, z[3,1]*I3, z[3,2]*I3, z[3,3]*I3)
        )
    end

    #Constraint
    function Φ(q)
        q_reshaped = reshape(q,3,4)
        ri,rj,u,v = [q_reshaped[:,k] for k = 1:4]
        [
        (rj-ri)⋅(rj-ri) - 1.0^2,
        (rj-ri)⋅u,
        (rj-ri)⋅v,
        u⋅v,
        u⋅u - 1.0^2,
        v⋅v - 1.0^2
        ]
    end

    function Φq(q)
        q_reshaped = reshape(q,3,4)
        ri,rj,u,v = [q_reshaped[:,k] for k = 1:4]
        ret = zeros(eltype(q), 6, 12)
        ret[1,1:3] = transpose(2(ri-rj))
        ret[1,4:6] = transpose(2(rj-ri))

        ret[2,1:3] = transpose(-u)
        ret[2,4:6] = transpose( u)
        ret[2,7:9] = transpose(rj-ri)

        ret[3,1:3] = transpose(-v)
        ret[3,4:6] = transpose( v)
        ret[3,10:12] = transpose(rj-ri)

        ret[4,7:9] = transpose(v)
        ret[4,10:12] = transpose(u)

        ret[5,7:9] = transpose(2u)

        ret[6,10:12] = transpose(2v)
        ret
    end
end
