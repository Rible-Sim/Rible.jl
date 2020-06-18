struct MLCP{T}
    A::Matrix{T}
    B::Matrix{T}
    C::Matrix{T}
    D::Matrix{T}
    a::Vector{T}
    b::Vector{T}
    n::Int
    m::Int
    u::Vector{T}
    v::Vector{T}
end

function pgs(problem::MLCP{T}) where T
    @unpack A,B,C,D,a,b,u,v,n,m = problem
    itermax = 1000
    pgsExplicit = true
    tol = 1e-7

    diagA = zeros(T,n)
    diagB = zeros(T,m)
    eps2 = eps(T)^2
    for i = 1:n
        if abs2(A[i,i]) < eps2
            @error "Vanishing diagonal term"
            return nothing
        else
            diagA[i] = 1.0/A[i,i]
        end
    end
    for i = 1:m
        if abs2(B[i,i]) < eps2
            @error "Vanishing diagonal term"
            return nothing
        else
            diagB[i] = 1.0/B[i,i]
        end
    end

    # start iteration
    iter = 0
    err = mlcp_compute_error(problem)
    while (iter < itermax && err < tol)
        iter += 1
        for i = 1:n
            u[i] = 0
            u[i] = -(a[i] + A[i,:]⋅u + C[i,:]⋅v)*diagA[i]
        end
        for i = 1:m
            v[i] = 0
              vi = -(b[i] + D[:,i]⋅u + B[:,i]⋅v)*diagB[i]
            v[i] = ifelse( vi < 0, 0.0, vi)
        end
    end

    if err > tol
        @warn """
        No convergence of PGS after $iter iterations. The residue is : $err"""
    end

    return nothing
end

function mlcp_compute_error(A,B,C,D,a,b,u,v)
    n = length(u)
    m = length(v)
    we = a + A*u + C*v
    wi = b + D*u + B*v
    ei = zero(wi)
    for j = 1:m
        vj  = v[j]
        wij = wi[j]
        if wij <= 0
            if vj > 0
                ei[j] = wij^2
            else #vj <= 0
                ei[j] = vj^2+wij^2
            end
        else #wij > 0
            if vj <=0
                ei[j] = vj^2
            else #vj > 0
                ei[j] = min(vj,wij)^2
            end
        end
    end

    mlcp_error = sqrt(sum(we.^2) + sum(ei))
end

function mlcp_compute_error(problem::MLCP{T}) where T
    @unpack a,b,A,B,C,D,u,v,m,n = problem
    we = a + A*u + C*v
    wi = b + D*u + B*v
    ei = zero(wi)
    for j = 1:m
        vj  = v[j]
        wij = wi[j]
        if wij <= 0
            if vj > 0
                ei[j] = wij^2
            else #vj <= 0
                ei[j] = vj^2+wij^2
            end
        else #wij > 0
            if vj <=0
                ei[j] = vj^2
            else #vj > 0
                ei[j] = min(vj,wij)^2
            end
        end
    end

    mlcp_error = sqrt(sum(we.^2) + sum(ei))
end

using PATHSolver

M = [0  0 -1 -1 ;
     0  0  1 -2 ;
     1 -1  2 -2 ;
     1  2 -2  4 ]

q = [2; 2; -2; -6]

myfunc(x) = M*x + q

n = 4
lb = zeros(n)
ub = 100*ones(n)

options(convergence_tolerance=1e-2, output=:yes, time_limit=3600)

z, f = solveLCP(myfunc, M, lb, ub)
