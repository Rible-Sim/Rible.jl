# Mass matrices
Ī2J̄(::Union{NC2D2C,NC3D3C,NC2D4C,NC2D6C,NC3D6C},Ī)  = Ī
Ī2J̄(::NC3D12C,Ī) = 1/2*tr(Ī)*I-Ī

function Īg2z(nmcs,m,Īₘ,r̄ₘ)
    (;r̄i,invX̄) = nmcs.data
    cₘ = to_local_coords(nmcs,r̄ₘ)
    J̄g = Ī2J̄(nmcs,Īₘ)
    J̄o = J̄g + m*r̄ₘ*transpose(r̄ₘ)
    J̄i = J̄o - m*r̄i*transpose(r̄ₘ) - m*r̄ₘ*transpose(r̄i) + m*r̄i*transpose(r̄i)
    z = invX̄*J̄i*transpose(invX̄)
    Symmetric(z)
end


function build_M̄_std(nmcs::NC, m::T, Īₘ, r̄ₘ) where T
    nld = get_num_of_local_dims(nmcs)
    a = to_local_coords(nmcs, r̄ₘ)
    z = Īg2z(nmcs, m, Īₘ, r̄ₘ)
    M̄_std = zeros(T, 1 + nld, 1 + nld)
    M̄_std[1, 1] = m
    M̄_std[2:1+nld, 1] = m * a
    M̄_std[1, 2:1+nld] = M̄_std[2:1+nld, 1]
    M̄_std[2:1+nld, 2:1+nld] .= z
    M̄_std
end

function build_M̄(nmcs::NC, m::T, Īₘ, r̄ₘ) where T
    M̄_std = build_M̄_std(nmcs, m, Īₘ, r̄ₘ)
    cv = get_conversion_core(nmcs)
    M̄ = transpose(cv) * M̄_std * cv
end

function build_M̄_and_B(nmcs::NC{N,M,T,L,NCOORDS,NCOORDS2}, m::T, Īₘ, r̄ₘ) where {N,M,T,L,NCOORDS,NCOORDS2}
    num_of_dim = get_num_of_dims(nmcs)
    M̄_std = build_M̄_std(nmcs, m, Īₘ, r̄ₘ)
    cv = get_conversion_core(nmcs)
    M̄ = transpose(cv) * M̄_std * cv
    B̄ = M̄_std[[1], :] * cv
    Id = IMatrix(num_of_dim)
    B = kron(B̄, Id)
    M̄, B
end

## Mass matrices: standard
function make_M(nmcs::NC, m::T, Īₘ, r̄ₘ) where {T} # ami (area moment of inertia tensor)
    num_of_dim = get_num_of_dims(nmcs)
    ncoords = get_num_of_coords(nmcs)
    M̄ = build_M̄(nmcs, m, Īₘ, r̄ₘ)
    M = SMatrix{ncoords,ncoords}(kron(M̄, IMatrix(num_of_dim)))
end
