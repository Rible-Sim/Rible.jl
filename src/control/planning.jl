struct Chart{T,ϕcT,make_ψcT}
    nq::Int
    nX::Int
    xc::Vector{T}
    Uc::Matrix{T}
    ϕc::ϕcT
    make_ψc::make_ψcT
    ε::T
    cosα::T
    ρ::T
    σ::T
    neighbor_idx::Vector{Int}
    neighbor_vectors::Vector{Vector{T}}
end

function build_Jac_F(st)
    A = build_A(st)
    Φ = build_Φ(st)
    nq = st.ncoords
    function inner_F(x)
        q = @view x[1:nq]
        q̇ = @view x[nq+1:2nq]
        vcat(Φ(q), A(q) * q̇)
    end
    function inner_Jac_F(x)
        q = @view x[1:nq]
        q̇ = @view x[nq+1:2nq]
        A_matrix = A(q)
        ∂Aq̇∂q_matrix = ∂Aq̇∂q(st, q̇)
        Jac1 = hcat(A_matrix, zero(A_matrix))
        Jac2 = hcat(∂Aq̇∂q_matrix, A_matrix)
        Jac_F = vcat(Jac1, Jac2)
    end
    inner_F, inner_Jac_F
end

function Chart(st, qc::AbstractVector{T}, q̇c::AbstractVector{T};
                ε=zero(T), α=zero(T), ρ=one(T), σ=2one(T)) where T
    nc = get_ncstr(st)
    F, Jac_F = build_Jac_F(st)
    xc = vcat(qc, q̇c)
    Q, R = qr(transpose(Jac_F(xc)))
    Uc = Q[:, size(R, 1)+1:end]
    nX = size(Uc, 2)
    nq = length(qc)
    UcT = transpose(Uc)
    ϕc(x) = UcT * (x - xc)
    function make_ψc(y)
        function inner_ψc!(ψ,x)
            ψ[1:2nc] .= F(x)
            ψ[2nc+1:end] .= UcT * (x - xc) - y
        end
        function inner_Jac_ψc!(J,x)
            J[1:2nc,:] .= Jac_F(x)
            J[2nc+1:end,:] .= UcT
        end
        inner_ψc!, inner_Jac_ψc!
    end
    neighbor_idx = Int[]
    neighbor_vectors = Vector{Vector{T}}()
    @assert σ > ρ
    Chart(nq, nX,
        xc, Uc,
        ϕc, make_ψc,
        ε, cos(α), ρ, σ,
        neighbor_idx,neighbor_vectors)
end
