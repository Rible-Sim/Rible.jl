"""
Natural Coordinates Types

This module defines the concrete types for natural coordinate systems.
"""

"""
Data for local natural coordinates.
$(TYPEDEF)
"""
struct LNCData{N,M,T,L}
    r̄i::SArray{Tuple{N},T,1,N}
    X̄::SArray{Tuple{N,M},T,2,L}
    invX̄::SArray{Tuple{M,N},T,2,L}
end

function LNCData(r̄i::SVector,X̄::SMatrix)
    LNCData(r̄i,X̄,(pinv(X̄)))
end

const I2_Bool = IMatrix(2)
const I3_Bool = IMatrix(3)

"""
Natural Coordinates

`N`: dimension of space
`M`: local dimension of natural coordinates
`T`: floating point type
`L`: number of coordinates
$(TYPEDEF)
"""
struct NC{N,M,T,L,NCOORDS,NCOORDS2} <: AbstractNonminimalCoordinates{N,T}
    np::Int64
    nv::Int64
    data::LNCData{N,M,T,L}
    conversion_to_std::SMatrix{NCOORDS,NCOORDS,Int64,NCOORDS2}
    conversion_to_X::SMatrix{NCOORDS,NCOORDS,Int64,NCOORDS2}
    hessians::Vector{Symmetric{Int64, SparseMatrixCSC{Int64, Int64}}}
    num_of_cstr::Int64
    cstr_idx::Vector{Int64}
end

"""
2D Natural Coordinates.
$(TYPEDEF)
"""
const NC2D{M,T,L,NCOORDS,NCOORDS2} = NC{2,M,T,L,NCOORDS,NCOORDS2}

"""
3D Natural Coordinates.
$(TYPEDEF)
"""
const NC3D{M,T,L,NCOORDS,NCOORDS2} = NC{3,M,T,L,NCOORDS,NCOORDS2}

"""
2D 2-component (point-like) Natural Coordinates.
$(TYPEDEF)
"""
const NC2D2C{T} = NC{2,0,T,0,2,4}

"""
2D 4-component (bar-like) Natural Coordinates.
$(TYPEDEF)
"""
const NC2D4C{T} = NC{2,1,T,2,4,16}

"""
2D 6-component (rigid body) Natural Coordinates.
$(TYPEDEF)
"""
const NC2D6C{T} = NC{2,2,T,4,6,36}

"""
3D 3-component (point-like) Natural Coordinates.
$(TYPEDEF)
"""
const NC3D3C{T} = NC{3,0,T,0,3,9}

"""
3D 6-component (bar-like) Natural Coordinates.
$(TYPEDEF)
"""
const NC3D6C{T} = NC{3,1,T,3,6,36}

"""
3D 12-component (rigid body) Natural Coordinates.
$(TYPEDEF)
"""
const NC3D12C{T} = NC{3,3,T,9,12,144}
