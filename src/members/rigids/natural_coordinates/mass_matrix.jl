
# Mass matrices
Ī2J̄(::Union{NC2D2C,NC3D3C,NC2D6C,NC3D6C},Ī)  = Ī
Ī2J̄(::NC3D12C,Ī) = 1/2*tr(Ī)*I-Ī

function Īg2az(nmcs,m,Īg,mass_center)
    (;r̄i,invX̄) = nmcs
    a = invX̄*(mass_center-r̄i)
    J̄g = Ī2J̄(nmcs,Īg)
    J̄o = J̄g + m*mass_center*transpose(mass_center)
    J̄i = J̄o - m*r̄i*transpose(mass_center) - m*mass_center*transpose(r̄i) + m*r̄i*transpose(r̄i)
    z = invX̄*J̄i*transpose(invX̄)
    #@assert issymmetric(z)
    a,Symmetric(z)
end

## Mass matrices: standard
function make_M(cf::CoordinateFunctions,m::T,Īg,mass_center) where {T} # ami (area moment of inertia tensor)
    (;nmcs) = cf
    ndim = get_num_of_dims(nmcs)
    nld = get_num_of_local_dims(nmcs)
    ncoords = get_num_of_coords(nmcs)
    a,z = Īg2az(nmcs,m,Īg,mass_center)
    M_raw = zeros(T,1+nld,1+nld)
    M_raw[1,1] = m
    M_raw[2:1+nld,1] = m*a
    M_raw[1,2:1+nld] = M_raw[2:1+nld,1]
    M_raw[2:1+nld,2:1+nld] .= z
    M_std = kron(M_raw,IMatrix(ndim))
    cv = nmcs.conversion_to_std
    M = SMatrix{ncoords,ncoords}(transpose(cv)*M_std*cv)
end
