using SymEngine
using LinearAlgebra
c1,c2,c3 = symbols("c1 c2 c3")
I3 = Matrix(1I,3,3)

Ciuvw_raw = [1 c1 c2 c3]
Ciuvw = kron(Ciuvw_raw,I3)
CCiuvw = expand.(transpose(Ciuvw)*Ciuvw)
CCiuvw == kron(transpose(Ciuvw_raw)*Ciuvw_raw,I3)

Cijvw_raw = [(1-c1) c1 c2 c3]
Cijvw = kron(Cijvw_raw,I3)
transpose(Cijvw_raw)*Cijvw_raw
CCijvw = expand.(transpose(Cijvw)*Cijvw)


Cijkw_raw = [(1-c1-c2) c1 c2 c3]
CCijkw_raw=transpose(Cijkw_raw)*Cijkw_raw

Cijkl_raw = [(1-c1-c2-c3) c1 c2 c3]
CCijkl_raw = transpose(Cijkl_raw)*Cijkl_raw
expand.(CCijkl_raw[1,1])
