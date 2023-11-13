
function Lmat(q::AbstractVector)
    q0,q1,q2,q3 = q
    SA[
        -q1  q0  q3 -q2;
        -q2 -q3  q0  q1;
        -q3  q2 -q1  q0;
    ]
end

function Lᵀmat(q::AbstractVector)
    q0,q1,q2,q3 = q
    SA[
        -q1 -q2 -q3;
         q0 -q3  q2;
         q3  q0 -q1;
        -q2  q1  q0;
    ]
end

function ∂Lᵀη∂q(η::AbstractVector)
    η1,η2,η3 = η
    SA[
         0  -η1 -η2 -η3;
        η1    0  η3 -η2;
        η2  -η3   0  η1;
        η3   η2 -η1   0;
    ]
end

function quatVel2localAngular(q::AbstractVector,q̇::AbstractVector)
    Ω = 2Lmat(q)*q̇
end

function localAngular2quatVel(q::AbstractVector,Ω::AbstractVector)
    q̇ = Lᵀmat(q)*Ω./2
end

function Rmat(q::AbstractVector)
    q0,q1,q2,q3 = q
    v = SA[q1,q2,q3]
    q0² = q0^2
    q1² = q1^2
    q2² = q2^2
    q3² = q3^2
    Symmetric(SA[
        q0² + q1² - q2² - q3²     2q1*q2        2q1*q3;
           0        q0² - q1² + q2² - q3²       2q2*q3;
           0           0          q0² - q1² - q2² + q3²;
    ]) + 
    skew(2q0.*v)
end

function ∂Rη∂q(q::AbstractVector,η::AbstractVector)
    q0,q1,q2,q3 = q
    v = SA[q1,q2,q3]
    Q = q0*I + skew(v)
    Qη = Q*η    
    2hcat(
        Qη,
        transpose(v)*η*I - skew(Qη)
    )
end

function ∂Rᵀη∂q(q::AbstractVector,η::AbstractVector)
    q0,q1,q2,q3 = q
    v = SA[q1,q2,q3]
    Q = q0*I - skew(v)
    Qη = Q*η    
    2hcat(
        Qη,
        transpose(v)*η*I + skew(Qη)
    )
end