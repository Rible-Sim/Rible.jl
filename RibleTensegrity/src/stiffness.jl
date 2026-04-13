function static_kinematic_determine(
        ℬᵀ;
        atol = 0.0
    )
    num_of_dof,ncables = size(ℬᵀ)
    U,S,V = svd(ℬᵀ; full = true, alg = LinearAlgebra.QRIteration())
    # @show size(U),size(V),size(S)    
    rtol = min(size(ℬᵀ)...)*eps(real(float(one(eltype(ℬᵀ)))))
    tol = max(atol, rtol*S[1])
    # @show tol,S[1]
    rank_ℬᵀ = count(x -> x > tol, S)
    S_nonzero = S[begin:rank_ℬᵀ]
    dsi = ncables - rank_ℬᵀ
    dki = num_of_dof - rank_ℬᵀ
    # @show num_of_dof,ncables,rank_ℬᵀ,dsi,dki
    A = vcat(
        ℬᵀ,
        Matrix(-1I,ncables,ncables)
    )

    hr = hrep(A, zeros(num_of_dof+ncables),  BitSet(1:num_of_dof))
    ph = polyhedron(hr, lib)
    vr = vrep(ph)
    @assert npoints(vr) == 1
    # @show nrays(vr)
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

function optimize_maximum_stiffness(mat𝒦ps,vec𝒦m,vecI,A,b,nx; verbose=false)

    model = COSMO.Model()

    cstr1 = COSMO.Constraint(
        hcat(
            mat𝒦ps,zero(vecI),-vecI
        ), 
        vec𝒦m, 
        # COSMO.PsdConeTriangle,
        COSMO.PsdCone
    )

    cstr2 = COSMO.Constraint(
        A, 
        b, 
        COSMO.ZeroSet
    )

    cstr3 = COSMO.Constraint(
        Matrix(1.0I,nx,nx)[:,end-1] |> Diagonal, 
        zeros(nx), 
        COSMO.Nonnegatives
    )

    cstr = [cstr1,cstr2,cstr3]

    custom_settings = COSMO.Settings(
        verbose = verbose, 
        eps_abs = 1e-10, 
        eps_rel = 1e-10,
        max_iter = Int(1e5)
    )

    COSMO.assemble!(
        model, 
        zeros(nx,nx), 
        -Matrix(1.0I,nx,nx)[:,end], 
        cstr,
        settings = custom_settings
    )

    COSMO.optimize!(model)
end

function optimize_zero_stiffness(mat𝒦ps,vec𝒦m,vecI,A,b,nx,x_0; verbose=false)

    model = COSMO.Model()

    cstr1 = COSMO.Constraint(
        hcat(
            mat𝒦ps,zero(vecI)
        ), 
        vec𝒦m, 
        # COSMO.PsdConeTriangle,
        COSMO.PsdCone
    )

    cstr2 = COSMO.Constraint(
        A, 
        b, 
        COSMO.ZeroSet
    )

    cstr = [cstr1,cstr2]

    custom_settings = COSMO.Settings(
        verbose = verbose, 
        eps_abs = 1e-3, 
        # eps_rel = 1e-10,
        max_iter = 50000,
        rho = 1e-6,
        # sigma = 1e-3,
        # adaptive_rho = false,
        alpha = 0.5,
        # scaling = 0,
        # accelerator = COSMO.EmptyAccelerator
    )

    COSMO.assemble!(
        model, 
        zeros(nx,nx), 
        -Matrix(1.0I,nx,nx)[:,end], 
        cstr,
        settings = custom_settings
    )
    COSMO.warm_start_primal!(model, x_0)
    COSMO.optimize!(model)
end

function optimize_zero_stiffness_Clarabel(mat𝒦ps,vec𝒦m,vecI,Aeq,beq,nx,x_0; verbose=false)

    solver = Clarabel.Solver()

    A = sparse(
        [
            hcat(
                mat𝒦ps,zero(vecI)
            );
            Aeq; 
        ]
    )
    b = [
        vec𝒦m;
        beq;
    ]

    # Clarabel.jl cone specification
    cones = [
        Clarabel.PSDTriangleConeT(3),
        Clarabel.ZeroConeT(length(beq)), 
    ]

    settings = Clarabel.Settings(
        verbose = verbose, 
        time_limit = 5,
        tol_gap_abs = 1e-12, #absolute duality gap tolerance
        tol_gap_rel = 1e-12, #relative duality gap tolerance
        tol_feas = 1e-12, #feasibility check tolerance (primal and dual)
        tol_infeas_abs = 1e-12, #absolute infeasibility tolerance (primal and dual)
        tol_infeas_rel = 1e-12, #relative infeasibility tolerance (primal and dual)
    )

    Clarabel.setup!(
        solver, 
        sparse(zeros(nx,nx)), 
        -Matrix(1.0I,nx,nx)[:,end],
        A, 
        b,
        cones, 
        settings
    )

    solution = Clarabel.solve!(solver)
end

# function super_stability(bot)
#     q = RB.get_coords(bot.structure)
#     # q̌ = RB.get_free_coords(bot.structure)
#     Ǎ = RB.make_cstr_jacobian(bot.structure)(q)
#     Ň_ = RB.nullspace(Ǎ)
#     Ň = modified_gram_schmidt(Ň_)
#     Q̃ = RB.build_Q(bot.structure)
#     L̂ = RB.build_L̂(bot.structure)

#     # Left hand side
#     Q̃L̂ = Q̃*L̂

#     Bᵀ = -Q̃L̂
#     ℬᵀ = transpose(Ň)*Bᵀ

#     S,D = static_kinematic_determine(ℬᵀ)
#     k = RT.get_cables_stiffness(bot.structure)

#     ns = size(S,2)
#     for i = 1:ns
#         f = S[:,i]
#         # equivalent μ
#         # μ = l .- (f./k)
#         λ = inv(Ǎ*transpose(Ǎ))*Ǎ*Bᵀ*f
#         # @show f,λ
#         Ǩa = RB.cstr_forces_jacobian(bot.structure,λ)
#         𝒦a = transpose(Ň)*Ǩa*Ň |> sparse
#         D_𝒦a = ldl(𝒦a).D.diag |> sort

#         Ǩm, Ǩg = RB.build_Ǩm_Ǩg!(bot.structure,q,f,k)

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

