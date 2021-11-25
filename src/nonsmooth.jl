struct FischerBurmeister{T}
    ϵ::T
    X₁::T
    X₂::T
end
FischerBurmeister() = FischerBurmeister(1e-14,1.0,1.0)
function (f::FischerBurmeister)(x,y)
    √((x/f.X₁)^2+(y*f.X₂)^2+2f.ϵ^2) - (x/f.X₁+y*f.X₂)
end
