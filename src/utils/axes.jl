struct Axes{N,T,L}
    X::SMatrix{N,N,T,L}
end


function Base.getproperty(a::Axes{2}, p::Symbol) 
    if (p === :n) || (p === :normal)
        return a.X[:,2]
    elseif (p === :t) || (p === :tangent)
        return  a.X[:,1]
    else # fallback to getfield
        return getfield(a, p)
    end
end

function Base.getproperty(a::Axes{3}, p::Symbol) 
    if (p === :n) || (p === :normal)
        return a.X[:,3]
    elseif (p === :t) || (p === :tangent)
        return  a.X[:,1]
    elseif (p === :b) || (p === :bitangent)
        return  a.X[:,2]
    else # fallback to getfield
        return getfield(a, p)
    end
end

Axes(X::Axes) = X

# planar_orthonormal_axes(n)
function Axes(n::SVector{2})
    n /= norm(n)
    Axes(
        SMatrix{2,2}(
            -n[2], n[1], 
             n[1], n[2]
        )
    )
end

# spatial_orthonormal_axes(n)
function Axes(n::SVector{3})
    normal = SVector{3}(n) 
    normal /= norm(normal)
    tangent, bitangent = HouseholderOrthogonalization(normal)
    SMatrix{3,3}(
        tangent[1],   tangent[2],   tangent[3],
        bitangent[1], bitangent[2], bitangent[3],
        normal[1],    normal[2],    normal[3],
    ) |> Axes
end

function (Base.:*)(R::StaticArray{Tuple{N,N}},axes::Axes{N}) where N
    Axes(SMatrix{N,N}(R*axes.X))
end

function (Base.:*)(axes::Axes{N},R::StaticArray{Tuple{N,N}}) where N
    Axes(SMatrix{N,N}(axes.X*R))
end

function (Base.:*)(axes1::Axes{N},axes2::Axes{N}) where N
    Axes(axes1.X*axes2.X)
end

# work around
function (Base.:*)(n::StaticArray{Tuple{1,N}},axes::Axes{N}) where N
    Axes(SVector{N}(n)).X*axes
end

function (Base.:*)(axes::Axes{N},n::SVector{N}) where N
    axes.X*n
end
