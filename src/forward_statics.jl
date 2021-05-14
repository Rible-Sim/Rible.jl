
function forward(tg,startsols_input,
                    start_parameters_input,
                    target_parameters_input,
                    F = reshape(build_G!(tg),:,1))
    @var q[1:tg.ncoords]
    @var a[1:tg.nstrings]
    @var λ[1:tg.nconstraint]
    @var u[1:tg.nstrings]
    @var g[1:size(F,2)]
    polyq = 1.0q
    polya = 1.0a
    polyλ = 1.0λ
    polyu = 1.0u

    q0,a0,λ0 = startsols_input
    u0,g0 = start_parameters_input
    u1,g1 = target_parameters_input

    Φ = build_Φ(tg)
    A = build_A(tg)
    Q̃ = build_Q̃(tg)
    U = build_U(tg,polya,polyu)
    S = build_S(tg,polya,polyq)

    F = [-transpose(A(polyq))*polyλ  + F*g + Q̃*U*polyq;
        S;
        Φ(polyq)]
    # @show maximum(abs.(to_number.(subs(F, q=>q0, a=>a0, λ=>λ0, u=>u0, g=>g0))))
    # reset_forces!(tg)
    # update_strings_apply_forces!(tg)
    # l = [s.state.length for s = tg.strings]
    # u1 = u0
    # g0 = 1.0
    # g1 = 0.0

    start_parameters = vcat(u0,g0)
    target_parameters = vcat(u1,g1)
    startsols = [[q0; a0; λ0]]
    vg = [q,a,λ]
    # FFF = [subs(f, u=>l,g=>0.0) for f in F]
    Fsys = System(F; variable_groups=vg, parameters = [u;g])

    # ParameterHomotopy(Fsys,start_parameters,target_parameters),startsols
    # [subs(f, q=>q0,λ=>λ0,a=>a0,u=>u0,g=>1.0) for f in F]
    result = HomotopyContinuation.solve(Fsys, startsols; start_parameters, target_parameters, threading = false)
    # result = HomotopyContinuation.solve(Fsys, startsols)
end

function deform_forward(tgstruct,startsols_input,
                        start_parameters_input,
                        target_parameters_input)
    @var q[1:tgstruct.ncoords]
    @var a[1:tgstruct.nstrings]
    @var λ[1:tgstruct.nconstraint]
    @var d[1:tgstruct.nconstraint]
    @var u[1:tgstruct.nstrings]
    @var g
    polyq = 1.0q
    polya = 1.0a
    polyλ = 1.0λ
    polyd = 1.0d
    polyu = 1.0u

    q0,a0,λ0 = startsols_input
    d0,u0,g0 = start_parameters_input
    d1,u1,g1 = target_parameters_input

    Φ = build_Φ(tgstruct)
    A = build_A(tgstruct)
    Q̃ = build_Q̃(tgstruct)
    U = build_U(tgstruct,polya,polyu)
    S = build_S(tgstruct,polya,polyq)

    G = build_G!(tgstruct)

    F = [-transpose(A(polyq))*polyλ  + g*G + Q̃*U*polyq;
        S;
        Φ(polyq,polyd)]

    @debug "$(maximum(abs.(to_number.(subs(F, q=>q0, a=>a0, λ=>λ0, d=>d0, u=>u0, g=>g0)))))"

    start_parameters = vcat(d0,u0,g0)
    target_parameters = vcat(d1,u1,g1)
    startsols = [[q0; a0; λ0]]
    vg = [q,a,λ]

    Fsys = System(F; variable_groups=vg, parameters = [d;u;g])

    result = HomotopyContinuation.solve(Fsys, startsols; start_parameters, target_parameters, threading = false)
end

function all_forward(tgstruct,startsols_input,
                        start_parameters_input,
                        target_parameters_input)
    @var q[1:tgstruct.ncoords]
    @var a[1:tgstruct.nstrings]
    @var λ[1:tgstruct.nconstraint]
    @var d[1:tgstruct.nconstraint]
    @var k[1:tgstruct.nstrings]
    @var u[1:tgstruct.nstrings]
    @var g
    polyq = 1.0q
    polya = 1.0a
    polyλ = 1.0λ
    polyd = 1.0d
    polyk = 1.0k
    polyu = 1.0u

    q0,a0,λ0 = startsols_input
    d0,k0,u0,g0 = start_parameters_input
    d1,k1,u1,g1 = target_parameters_input

    Φ = build_Φ(tgstruct)
    A = build_A(tgstruct)
    Q̃ = build_Q̃(tgstruct)
    U = build_U(tgstruct,polya,polyu,polyk)
    S = build_S(tgstruct,polya,polyq)

    G = build_G!(tgstruct)

    F = [-transpose(A(polyq))*polyλ  + g*G + Q̃*U*polyq;
        S;
        Φ(polyq,polyd)]

    @debug "$(maximum(abs.(to_number.(subs(F, q=>q0, a=>a0, λ=>λ0, d=>d0, k=>k0, u=>u0, g=>g0)))))"

    start_parameters = vcat(d0,k0,u0,g0)
    target_parameters = vcat(d1,k1,u1,g1)
    startsols = [[q0; a0; λ0]]
    vg = [q,a,λ]

    Fsys = System(F; variable_groups=vg, parameters = [d;k;u;g])

    result = HomotopyContinuation.solve(Fsys, startsols; start_parameters, target_parameters, threading = false)
end
