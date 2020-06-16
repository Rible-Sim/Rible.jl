function build_T(st2d)
    nq_per_body = 4
    body2q = st2d.connectivity.body2q
    nbody = st2d.nbody
    nq = body2q[end][end]
    T = spzeros(Int,nq,nbody*nq_per_body)
    for (rbid,rb) in enumerate(st2d.rigidbodies)
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
function build_C(st2d)
    nq_per_body = 4
    @unpack nbody,npoints,ndim = st2d
    C = spzeros(nbody*nq_per_body,npoints*ndim)
    js = 0
    for (rbid,rb) in enumerate(st2d.rigidbodies)
        is = (rbid-1)*nq_per_body
        for (apid,Cp) in enumerate(rb.state.auxs.Cp)
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
# manipulator.connectivity.string2bp[1]

function build_D(st2d)
    @unpack nstring,npoints,ndim = st2d
    @unpack body2q,string2bp = st2d.connectivity
    D = spzeros(Int,npoints*ndim,nstring*ndim)
    D_raw = spzeros(Int,npoints,nstring)
    iss = [0]
    for (rbid,rb) in enumerate(st2d.rigidbodies)
        push!(iss,iss[end]+rb.prop.number_aps)
    end

    for (sid,bp) in enumerate(string2bp)
        D_raw[iss[bp[1].rbid]+bp[1].apid,sid] = 1
        D_raw[iss[bp[2].rbid]+bp[2].apid,sid] = -1
    end
    D .= kron(D_raw,Matrix(1I,2,2))
end

function build_Q̃(st2d)
    Q̃=build_T(st2d)*build_C(st2d)*build_D(st2d)
end

function fvector(st2d)
    @unpack strings = st2d
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
