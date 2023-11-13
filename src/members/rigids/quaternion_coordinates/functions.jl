to_local_coords(::QC,c) = c

function to_transformation(::QC,c)
    ĉ = skew(c)
    function inner_C(x)
        q = @view x[4:7]
        T = eltype(q)
        I3 = SMatrix{3,3}(one(T)*I)
        hcat(
            I3,
            -2Rmat(q)*ĉ*Lmat(q)
        )
    end
end

function make_∂Cẋ∂x(c)
    ĉ = skew(c)
    function ∂Cẋ∂x(x,ẋ)
        q = @view x[4:7]
        q̇ = @view ẋ[4:7]
        R = Rmat(q)
        ĉL̇ = ĉ*Lmat(q̇)
        η = ĉL̇*q
        T = eltype(q)
        O3 = @SMatrix zeros(T,3,3)
        hcat(
            O3,
            2(R*ĉL̇ + ∂Rη∂q(q,η))
        )        
    end
end

function make_∂Cᵀf∂x(c)
    ĉᵀ = transpose(skew(c))
    function ∂Cᵀf∂x(x,f)
        q = @view x[4:7]
        ĉᵀRᵀf = ĉᵀ*transpose(Rmat(q))*f
        T = eltype(q)
        O37 = @SMatrix zeros(T,3,7)
        O43 = @SMatrix zeros(T,4,3)
        vcat(
            O37,
            hcat(
                O43,
                -2(Lᵀmat(q)*ĉᵀ*∂Rᵀη∂q(q,f) + ∂Lᵀη∂q(ĉᵀRᵀf))
            )
        )        
    end
end


function make_∂Cẋ∂x_forwarddiff(C,nc,nx)
    function ∂Cẋ∂x(x,ẋ)
        function Cẋ(_x)
            C(_x)*ẋ
        end
        out = zeros(eltype(ẋ),nc,nx)
        ForwardDiff.jacobian!(out,Cẋ,x)
    end
end

function cartesian_frame2coords(::QC,origin_position,R,origin_velocity,ω)
    Rmat = RotMatrix(R)
    q = QuatRotation(Rmat).q |> vec
    Ω = inv(Rmat)*ω
    q̇ = localAngular2quatVel(q,Ω)
    x = vcat(origin_position,q)
    ẋ = vcat(origin_velocity,q̇)
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