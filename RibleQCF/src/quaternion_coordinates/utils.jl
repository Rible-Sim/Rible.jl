
function Lmat(q::AbstractVector)
    q0,q1,q2,q3 = q
    SA[
        -q1  q0  q3 -q2;
        -q2 -q3  q0  q1;
        -q3  q2 -q1  q0;
    ]
end

Lᵀmat(q::AbstractVector) = transpose(Lmat(q))

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

"""
    Rmat(q::AbstractVector)

rotation matrix
"""
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

Rmat(q::Quaternion) = Rmat(vec(q))

function ∂Rη∂q(q::AbstractVector,η::AbstractVector)
    q0,q1,q2,q3 = q
    v = SA[q1,q2,q3]
    Qη = q0*η - skew(η)*v
    2hcat(
        Qη,
        transpose(η)*v*I - skew(Qη)
    )
end

∂Rη∂q(q::Quaternion,η::AbstractVector) = ∂Rη∂q(vec(q),η)

function ∂²Rη∂qᵀ∂q(η::AbstractVector)
    H = -skew(η)
    hessians = [
        begin
            rb = @SMatrix [
                ifelse(
                    (i==k),
                    η[j],
                    0
                ) + 
                ifelse(
                    (j==k),
                    η[i],
                    0
                )
                for i = 1:3, j = 1:3
            ]
            2vcat(
                hcat(η[k],transpose(H[k,:])),
                hcat(H[k,:],rb-η[k]I)
            )
        end
        for k = 1:3
    ]
end

function ∂Rᵀf∂q(q::AbstractVector,f::AbstractVector)
    q0,q1,q2,q3 = q
    v = SA[q1,q2,q3]
    Q = q0*I - skew(v)
    Qf = Q*f    
    2hcat(
        Qf,
        transpose(v)*f*I + skew(Qf)
    )
end

∂Rᵀf∂q(q::Quaternion,f::AbstractVector) = ∂Rᵀf∂q(vec(q),f)

function Pmat(q::AbstractVector)
    q0,q1,q2,q3 = q
    SA[
        q0 -q1 -q2 -q3;
        q1  q0 -q3  q2;
        q2  q3  q0 -q1;
        q3 -q2  q1  q0;
    ]
end

Pmat(q::Quaternion) = Pmat(vec(q))

function Mmat(q::AbstractVector)
    q0,q1,q2,q3 = q
    SA[
        q0 -q1 -q2 -q3;
        q1  q0  q3 -q2;
        q2 -q3  q0  q1;
        q3  q2 -q1  q0;
    ]
end

Mmat(q::Quaternion) = Mmat(vec(q))

const Inv_mat = Diagonal(SA[1,-1,-1,-1])