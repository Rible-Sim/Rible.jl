function static_kinematic_determine(
        ℬᵀ;
        atol = 0.0
    )
    ndof,ncables = size(ℬᵀ)
    U,S,V = svd(ℬᵀ; full = true, alg = LinearAlgebra.QRIteration())
    @show size(U),size(V),size(S)    
    rtol = min(size(ℬᵀ)...)*eps(real(float(one(eltype(ℬᵀ)))))
    tol = max(atol, rtol*S[1])
    @show tol,S[1]
    rank_ℬᵀ = count(x -> x > tol, S)
    S_nonzero = S[begin:rank_ℬᵀ]
    dsi = ncables - rank_ℬᵀ
    dki = ndof - rank_ℬᵀ
    @show ndof,ncables,rank_ℬᵀ,dsi,dki
    A = vcat(
        ℬᵀ,
        Matrix(-1I,ncables,ncables)
    )

    hr = hrep(A, zeros(ndof+ncables),  BitSet(1:ndof))
    ph = polyhedron(hr, lib)
    vr = vrep(ph)
    @assert npoints(vr) == 1
    @show nrays(vr)
    self_stress_states = reduce(hcat,[ray.a for ray in rays(vr)])
    # self_stress_states = V[:,rank_ℬᵀ+1:rank_ℬᵀ+dsi]
    stiffness_directions = U[:,rank_ℬᵀ+1:rank_ℬᵀ+dki]
    self_stress_states,stiffness_directions
end

function classical_gram_schmidt(matrix)
    # orthogonalises the columns of the input matrix
    num_vectors = size(matrix)[2]
    orth_matrix = copy(matrix)
    for vec_idx = 1:num_vectors
        for span_base_idx = 1:(vec_idx-1)
            # substrack sum
            orth_matrix[:, vec_idx] -= dot(orth_matrix[:, span_base_idx], orth_matrix[:, vec_idx])*orth_matrix[:, span_base_idx]
        end
        # normalise vector
        orth_matrix[:, vec_idx] = orth_matrix[:, vec_idx]/norm(orth_matrix[:, vec_idx])
    end
    return orth_matrix
end

function modified_gram_schmidt(matrix)
    # orthogonalises the columns of the input matrix
    num_vectors = size(matrix)[2]
    orth_matrix = copy(matrix)
    for vec_idx = 1:num_vectors
        orth_matrix[:, vec_idx] = orth_matrix[:, vec_idx]/norm(orth_matrix[:, vec_idx])
        for span_base_idx = (vec_idx+1):num_vectors
            # perform block step
            orth_matrix[:, span_base_idx] -= dot(orth_matrix[:, span_base_idx], orth_matrix[:, vec_idx])*orth_matrix[:, vec_idx]
        end
    end
    return orth_matrix
end

function optimize_maximum_stiffness(mat𝒦ps,vec𝒦m,vecI,A,b,nx)

    model = COSMO.Model()

    constraint1 = COSMO.Constraint(
        hcat(
            mat𝒦ps,zero(vecI),-vecI
        ), 
        vec𝒦m, 
        # COSMO.PsdConeTriangle,
        COSMO.PsdCone
    )

    constraint2 = COSMO.Constraint(
        A, 
        b, 
        COSMO.ZeroSet
    )

    constraints = [constraint1,constraint2]

    custom_settings = COSMO.Settings(
        verbose = true, 
        eps_abs = 1e-7, 
        eps_rel = 1e-7
    )

    COSMO.assemble!(
        model, 
        zeros(nx,nx), 
        -Matrix(1.0I,nx,nx)[:,end], 
        constraints,
        settings = custom_settings
    )

    COSMO.optimize!(model)
end

function optimize_minimum_stiffness(mat𝒦ps,vec𝒦m,vecI,A,b,nx,)

    model = COSMO.Model()

    constraint1 = COSMO.Constraint(
        hcat(
            mat𝒦ps,zero(vecI)
        ), 
        vec𝒦m, 
        # COSMO.PsdConeTriangle,
        COSMO.PsdCone
    )

    constraint2 = COSMO.Constraint(
        A, 
        b, 
        COSMO.ZeroSet
    )

    # constraint3 = COSMO.Constraint(
    #     Matrix(1.0I,nx,nx)[:,size(mat𝒦ps,2)+1]', 
    #     [-σmax], 
    #     COSMO.Nonnegatives
    # )

    constraints = [constraint1,constraint2,]

    custom_settings = COSMO.Settings(
        verbose = true, 
        eps_abs = 1e-7, 
        eps_rel = 1e-7
    )

    COSMO.assemble!(
        model, 
        zeros(nx,nx), 
        -Matrix(1.0I,nx,nx)[:,end], 
        constraints,
        settings = custom_settings
    )
    # COSMO.warm_start_primal!(model, x_0)
    COSMO.optimize!(model)
end

# function super_stability(bot)
#     q = TR.get_q(bot.tg)
#     # q̌ = TR.get_q̌(bot.tg)
#     Ǎ = TR.make_A(bot.tg)(q)
#     Ň_ = TR.nullspace(Ǎ)
#     Ň = modified_gram_schmidt(Ň_)
#     Q̃ = TR.build_Q̃(bot.tg)
#     L̂ = TR.build_L̂(bot.tg)

#     # Left hand side
#     Q̃L̂ = Q̃*L̂

#     Bᵀ = -Q̃L̂
#     ℬᵀ = transpose(Ň)*Bᵀ

#     S,D = static_kinematic_determine(ℬᵀ)
#     k = TR.get_cables_stiffness(bot.tg)

#     ns = size(S,2)
#     for i = 1:ns
#         f = S[:,i]
#         # equivalent μ
#         # μ = l .- (f./k)
#         λ = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*f
#         # @show f,λ
#         Ǩa = - TR.∂Aᵀλ∂q̌(bot.tg,λ)
#         𝒦a = transpose(Ň)*Ǩa*Ň |> sparse
#         D_𝒦a = ldl(𝒦a).D.diag |> sort

#         Ǩm, Ǩg = TR.build_Ǩm_Ǩg!(bot.tg,q,f,k)

#         𝒦m = transpose(Ň)*Ǩm*Ň |> sparse
#         𝒦g = transpose(Ň)*Ǩg*Ň |> sparse
#         D_𝒦g = ldl(𝒦g).D.diag |> sort

#         𝒦ga = 𝒦g.+ 𝒦a |> sparse
#         D_𝒦ga = ldl(𝒦ga).D.diag |> sort
#         println("the $i self-stress state: ")
#         @show f,count(>(1e-11), f)
#         @show D_𝒦a
#         @show D_𝒦g
#         @show D_𝒦ga
#     end
# end

