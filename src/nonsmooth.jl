struct FischerBurmeister{T}
    ϵ::T
end
function (f::FischerBurmeister)(x,y)
    √(x^2+y^2+2f.ϵ^2) - (x+y)
end
