to_local_coords(::QC,c) = c

function to_position(::QC,x,c)
    ro = @view x[1:3]
    q = @view x[4:7]
    ro .+ Rmat(q)*c
end

function to_transformation(::QC,x,c)
    q = @view x[4:7]
    T = eltype(q)
    I3 = SMatrix{3,3}(one(T)*I)
    hcat(
        I3,
        ∂Rη∂q(q,c)
    )
end

function ∂Cẋ∂x(x,ẋ,c)
    η = SA[c[1],c[2],c[3]]
    q̇ = @view ẋ[4:7]
    q̇0 = q̇[1]
    v̇ = SA[q̇[2],q̇[3],q̇[4]]
    T = eltype(ẋ)
    O3 = @SMatrix zeros(T,3,3)
    hcat(
        O3,
        2hcat(
             q̇0*c+skew(v̇)*η,
            -q̇0*skew(η)+kron(transpose(η),v̇)-skew(v̇)*skew(η)
        )
    )
end

function ∂Cᵀf∂x(x,f,c)
    η = SA[c[1],c[2],c[3]]
    T = eltype(x)
    O37 = @SMatrix zeros(T,3,7)
    O43 = @SMatrix zeros(T,4,3)
    hessians = ∂²Rη∂qᵀ∂q(η)
    vcat(
        O37,
        hcat(
            O43,
            sum(
                hessians[k]*f[k]
                for k = 1:3
            )
        )
    )
end

function cartesian_frame2coords(::QC,origin_frame)
    (;position,velocity,axes,angular_velocity) = origin_frame
    rotation_matrix = RotMatrix(axes.X)
    q = QuatRotation(rotation_matrix).q |> vec
    Ω = inv(rotation_matrix)*angular_velocity
    q̇ = localAngular2quatVel(q,Ω)
    x = vcat(position,q)
    ẋ = vcat(velocity,q̇)
    x, ẋ
end

function find_rotation(::QC,x)
    q = @view x[4:7]
    # RotMatrix(QuatRotation(q))
    q0, q1, q2, q3 = q
    q0², q1², q2², q3² = q0^2, q1^2, q2^2, q3^2
    SA[
        2(q0²+q1²)-1   2(q1*q2-q0*q3) 2(q1*q3+q0*q2);
        2(q1*q2+q0*q3) 2(q0²+q2²)-1   2(q2*q3-q0*q1);
        2(q1*q3-q0*q2) 2(q2*q3+q0*q1) 2(q0²+q3²)-1;
    ]
end

function find_local_angular_velocity(::QC,x,ẋ)
    q = @view x[4:7] 
    q̇ = @view ẋ[4:7]
    quatVel2localAngular(q,q̇)
end

function find_angular_velocity(qcs::QC,x,ẋ)
    R = find_rotation(qcs,x)
    Ω = find_local_angular_velocity(qcs,x,ẋ)
    R*Ω
end