function build_T(tgstruct)
    nq_per_body = 4
    body2q = tgstruct.connectivity.body2q
    nbody = tgstruct.nbody
    nq = body2q[end][end]
    T = spzeros(Int,nq,nbody*nq_per_body)
    for (rbid,rb) in enumerate(tgstruct.rigidbodies)
        pindex = body2q[rbid]
        js = (rbid-1)*nq_per_body
        for (j,index) in enumerate(pindex)
            T[index,js+j] = 1
        end
    end
    T
end
# T = build_T(manipulator)
# Array(T)
# @code_warntype build_T(manipulator)
function build_C(tgstruct)
    nq_per_body = 4
    @unpack nbody,npoints,ndim = tgstruct
    C = spzeros(nbody*nq_per_body,npoints*ndim)
    js = 0
    for (rbid,rb) in enumerate(tgstruct.rigidbodies)
        is = (rbid-1)*nq_per_body
        for (apid,Cp) in enumerate(rb.state.cache.Cp)
            jlen,_ = size(Cp)
            C[is+1:is+nq_per_body,js+1:js+jlen] .= transpose(Cp)
            # @show is+1:is+nq_per_body, js+1:js+jlen
            js += jlen
        end
    end
    C
end
# C = build_C(manipulator)
# Array(C)
# @code_warntype build_C(manipulator)
# manipulator.connectivity.string2ap[1]

function build_D(tgstruct)
    @unpack nstring,npoints,ndim = tgstruct
    @unpack body2q,string2ap = tgstruct.connectivity
    D = spzeros(Int,npoints*ndim,nstring*ndim)
    D_raw = spzeros(Int,npoints,nstring)
    iss = [0]
    for (rbid,rb) in enumerate(tgstruct.rigidbodies)
        push!(iss,iss[end]+rb.prop.naps)
    end

    for (sid,ap) in enumerate(string2ap)
        D_raw[iss[ap[1].rbid]+ap[1].apid,sid] = 1
        D_raw[iss[ap[2].rbid]+ap[2].apid,sid] = -1
    end
    D .= kron(D_raw,Matrix(1I,2,2))
end

function build_Q̃(tgstruct)
    Q̃=build_T(tgstruct)*build_C(tgstruct)*build_D(tgstruct)
end

function fvector(tgstruct)
    @unpack strings = tgstruct
    ret = vcat([str.state.tension*str.state.direction for str in strings]...)
end

function iksolve(prob;ftol=1e-14)
    @unpack funcs,q0,u0,λ0,nq,nu,nλ = prob
    A,F! = funcs
    F = zeros(nq)
    function IK_R!(R,x)
        u = @view x[   1:nu]
        λ = @view x[nu+1:nu+nλ]
        F!(F,u)
        R .= transpose(A(q0))*λ - F
    end
    initial_x = vcat(u0,λ0)
    ik_result = nlsolve(IK_R!,initial_x,ftol=ftol)
    @info "Convergence: $(converged(ik_result)), Iterations: $(ik_result.iterations), ftol=$ftol"
    u_result = ik_result.zero[   1:nu]
    λ_result = ik_result.zero[nu+1:nu+nλ]
    u_result,λ_result
end
# D = build_D(manipulator)
# Array(D)

# Q̃=build_T(manipulator)*build_C(manipulator)*build_D(manipulator)
# Array(Q̃)
