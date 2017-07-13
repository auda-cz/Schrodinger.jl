"""
    rotation(θ,n=(1,0,0))

Generate a qubit rotation operator about an axis defined by the vector \$\\vec{n}\$. Note that the `n` vector will be normalized, allowing for inputs like `(1,1,0)`. The rotation operator is defined as
```math
\\hat{R}_{\\vec{n}}(θ) = \\exp\\left(−iθ\\vec{n}⋅\\vec{σ}/2\\right) = \\cos\\frac{θ}{2}𝟙 - i\\sin\\frac{θ}{2}(n_xσ_x + n_yσ_y + n_zσ_z).
```
"""
function rotation(θ::Real,n::NTuple{3,Int}=(1,0,0))
    R = Matrix{Complex128}(2,2)
    nx,ny,nz = n
    a = 1/√(nx^2+ny^2+nz^2)
    nx,ny,nz = a*nx,a*ny,a*nz
    c = cos(0.5θ); s = sin(0.5θ)
    R[1,1] = Complex(c,-nz*s)
    R[2,1] = Complex(ny*s,-nx*s)
    R[1,2] = Complex(-ny*s,-nx*s)
    R[2,2] = Complex(c,nz*s)
    return Operator(R,(2,))
end

const H = Operator(Float64[1/√2 1/√2; 1/√2 -1/√2],(2,),true)
const S = Operator(sparse(Complex128[1 0; 0 1im]),(2,),false)
const T = Operator(sparse(Complex128[1 0; 0 Complex(1/√2,1/√2)]),(2,),false)

const cNOT = Operator(sparse(Float64[1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0]),(2,2),true)
const rcNOT = Operator(sparse(Float64[0 1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]),(2,2),true)
const cZ = Operator(sparse(Float64[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 -1]),(2,2),true)
end
