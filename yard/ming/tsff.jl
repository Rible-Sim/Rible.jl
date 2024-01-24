function get_ΔLs(ll,l0,L)
    Diagonal((ll.-l0)./ll)*L
end


function get_rotation_matrix(ω)
    normω =norm(ω)
    𝛚 = RB.skew(ω)
    R = I + sin(normω)*𝛚./normω + (1-cos(normω))*(𝛚^2)./normω
end

function update_N!(N,nc,Δn,R)
    for i in 1:size(N,1)
        N[i,:] = transpose(R*(N[i,:]-transpose(nc))) + nc + Δn
    end
end

function initialize_tsff(st)
    rbs = st.rigidbodies
    ss = st.strings
    cnt = st.connectivity
    @unpack nbodies,nstrings = st
    @unpack string2ap = cnt
    nnodes_by_body = [body.prop.naps for rb in rbs]
    N_by_body = [reduce(vcat,[transpose(rp) for rp in body.state.loci_states]) for rb in rbs]
    C_by_body = [zeros(Int,nstrings,num_of_loci) for num_of_loci in nnodes_by_body]
    for sstring in ss
        (;id)= sstring
        a,b = string2ap[id]
        C_by_body[a.bodyid][id,a.apid] = -1
        C_by_body[b.bodyid][id,b.apid] =  1
    end
    𝐤 = [sstring.k for sstring in ss]
    K = Diagonal(𝐤)
    l0 = RB.get_strings_restlen(st)
    nc = [transpose(body.state.ro[:]) for rb in rbs]
    nbodies,N_by_body,C_by_body,nnodes_by_body,K,l0,nc
end

function update_FtTt!(Ft,Tt,nb,N_by_body,C_by_body,nnodes_by_body,K,l0,nc)
    L_by_body = [C*N for (C,N) = zip(C_by_body,N_by_body)]
    L = sum(L_by_body)
    ll = [norm(L[i,:]) for i in 1:size(L,1)]
    # @show ll
    ΔLs = get_ΔLs(ll,l0,L)
    for bodyid = 1:nb
        F = -transpose(C_by_body[bodyid])*K*ΔLs
        N = N_by_body[bodyid]
        T = reduce(vcat,transpose((N[i,:]-transpose(nc[bodyid]))×F[i,:]) for i in 1:nnodes_by_body[bodyid])
        Ft[bodyid] = sum(F,dims=1)
        Tt[bodyid] = sum(T,dims=1)
    end
end

function tsff(st;e0=1e-4,maxiters=10000)
    rbs = st.rigidbodies
    m = 0.1
    J = 10
    c = 100
    Δt = 0.0001
    nb,N_by_body,C_by_body,nnodes_by_body,K,l0,nc = initialize_tsff(st)
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
            for bodyid in 2:nb
                # @show "before",vt[bodyid], wt[bodyid]
                vt[bodyid] .= (Ft[bodyid] .+ (m/Δt-c/2).*vt[bodyid])./(m/Δt+c/2)
                wt[bodyid] .= (Tt[bodyid] .+ (J/Δt-c/2).*wt[bodyid])./(J/Δt+c/2)
                # @show "after",vt[bodyid], wt[bodyid]
                Δn = Δt.*vt[bodyid]
                ΔΦ = Δt.*wt[bodyid]
                R = get_rotation_matrix(ΔΦ)
                # @show R
                update_N!(N_by_body[bodyid],nc[bodyid],Δn,R)
                nc[bodyid] .+= Δn
                # @show N_by_body[bodyid]
            end
            push!(Ns,deepcopy(N_by_body))
        end
    end
    errs,Ns
end
