function get_ΔLs(ll,l0,L)
    Diagonal((ll.-l0)./ll)*L
end


function get_rotation_matrix(ω)
    normω =norm(ω)
    𝛚 = TR.NCF.skew(ω)
    R = I + sin(normω)*𝛚./normω + (1-cos(normω))*(𝛚^2)./normω
end

function update_N!(N,nc,Δn,R)
    for i in 1:size(N,1)
        N[i,:] = transpose(R*(N[i,:]-transpose(nc))) + nc + Δn
    end
end

function initialize_tsff(tg)
    rbs = tg.rigidbodies
    ss = tg.strings
    cnt = tg.connectivity
    @unpack nbodies,nstrings = tg
    @unpack string2ap = cnt
    nnodes_by_body = [rb.prop.naps for rb in rbs]
    N_by_body = [reduce(vcat,[transpose(rp) for rp in rb.state.rps]) for rb in rbs]
    C_by_body = [zeros(Int,nstrings,nnodes) for nnodes in nnodes_by_body]
    for sstring in ss
        (;id)= sstring
        a,b = string2ap[id]
        C_by_body[a.rbid][id,a.apid] = -1
        C_by_body[b.rbid][id,b.apid] =  1
    end
    𝐤 = [sstring.k for sstring in ss]
    K = Diagonal(𝐤)
    l0 = TR.get_strings_restlen(tg)
    nc = [transpose(rb.state.ro[:]) for rb in rbs]
    nbodies,N_by_body,C_by_body,nnodes_by_body,K,l0,nc
end

function update_FtTt!(Ft,Tt,nb,N_by_body,C_by_body,nnodes_by_body,K,l0,nc)
    L_by_body = [C*N for (C,N) = zip(C_by_body,N_by_body)]
    L = sum(L_by_body)
    ll = [norm(L[i,:]) for i in 1:size(L,1)]
    # @show ll
    ΔLs = get_ΔLs(ll,l0,L)
    for rbid = 1:nb
        F = -transpose(C_by_body[rbid])*K*ΔLs
        N = N_by_body[rbid]
        T = reduce(vcat,transpose((N[i,:]-transpose(nc[rbid]))×F[i,:]) for i in 1:nnodes_by_body[rbid])
        Ft[rbid] = sum(F,dims=1)
        Tt[rbid] = sum(T,dims=1)
    end
end

function tsff(tg;e0=1e-4,maxiters=10000)
    rbs = tg.rigidbodies
    m = 0.1
    J = 10
    c = 100
    Δt = 0.0001
    nb,N_by_body,C_by_body,nnodes_by_body,K,l0,nc = initialize_tsff(tg)
    vt = [zeros(1,3) for i in 1:nb]
    wt = [zeros(1,3) for i in 1:nb]
    Ft = [zeros(1,3) for i in 1:nb]
    Tt = [zeros(1,3) for i in 1:nb]
    errs = Vector{Float64}()
    Ns = [N_by_body]
    for iteration = 1:maxiters
        update_FtTt!(Ft,Tt,nb,N_by_body,C_by_body,nnodes_by_body,K,l0,nc)
        err = max(maximum(norm.(Ft)),maximum(norm.(Tt)))
        push!(errs,err)
        # @show Ft
        # @show err
        if  err < e0
            break
        else
            for rbid in 2:nb
                # @show "before",vt[rbid], wt[rbid]
                vt[rbid] .= (Ft[rbid] .+ (m/Δt-c/2).*vt[rbid])./(m/Δt+c/2)
                wt[rbid] .= (Tt[rbid] .+ (J/Δt-c/2).*wt[rbid])./(J/Δt+c/2)
                # @show "after",vt[rbid], wt[rbid]
                Δn = Δt.*vt[rbid]
                ΔΦ = Δt.*wt[rbid]
                R = get_rotation_matrix(ΔΦ)
                # @show R
                update_N!(N_by_body[rbid],nc[rbid],Δn,R)
                nc[rbid] .+= Δn
                # @show N_by_body[rbid]
            end
            push!(Ns,deepcopy(N_by_body))
        end
    end
    errs,Ns
end
