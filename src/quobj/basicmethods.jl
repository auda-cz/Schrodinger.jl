import Base: complex, length, size, LinAlg.checksquare, getindex, setindex!,
     diag, full, norm, trace, normalize!, scale!, ishermitian, issymmetric, isdiag,
     similar, copy, hash, isequal, ==, isapprox, show

# Special QuObject methods
data(A::QuObject) = A.data
dims(A::QuObject) = A.dims
isnormalized(A::QuVector) = norm(A) ≈ 1.0
isnormalized(A::QuMatrix) = trace(A) ≈ 1.0
function dimsmatch(A::QuObject,B::QuObject)
    A.dims==B.dims || throw(DimensionMismatch("subspace dimensions do not match"))
    return nothing
end
dense(x::Ket) = Ket(full(x.data),x.dims)
dense(x::Bra) = Bra(full(x.data),x.dims)
dense(A::Operator) = Operator(full(A.data),A.dims)


# Translate basic Base array methods to QuObjects
length(A::QuObject) = length(A.data)
size(A::QuObject, d) = size(A.data, d)
size(A::QuObject) = size(A.data)
checksquare(A::QuMatrix) = checksquare(A.data)
getindex(A::QuObject, idx...) = getindex(A.data,idx...)
setindex!(A::QuObject, idx...) = setindex!(A.data,v,idx...)
diag(A::QuMatrix,k::Int=0) = diag(A.data,k)
full(A::QuObject) = full(A.data)
complex(x::Ket) = Ket(complex(x.data),x.dims)
complex(x::Bra) = Bra(complex(x.data),x.dims)
complex(A::Operator) = Operator(complex(A.data),A.dims)
norm(x::QuVector,n::Int=2) = norm(x.data,n)
trace(A::QuMatrix) = trace(A.data)
normalize!(x::QuVector) = (normalize!(x.data,2);x)
normalize!(A::QuMatrix) = (scale!(A.data,1/trace(A.data));A)
scale!(A::QuObject,b::Number) = (scale!(A.data,b);A)
ishermitian{T,D}(A::Operator{Herm,T,D}) = true
ishermitian{T,D}(A::Operator{NonHerm,T,D}) = false
issymmetric(A::QuMatrix) = issymmetric(A.data)
isdiag(A::QuMatrix) = isdiag(A.data)
similar{T<:QuObject}(A::T) = T(similar(A.data),A.dims)
copy{T<:QuObject}(A::T) = T(copy(A.data),A.dims)
hash{T<:QuObject}(A::T,h::UInt) = hash(hash(A.data,hash(A.dims,hash(T))),h)
isequal{T<:QuObject}(A::T,B::T) = isequal(A.dims,B.dims)&&isequal(A.data,B.data)
==(A::Ket,B::Ket) = isequal(A.dims,B.dims)&&(A.data==B.data)
==(A::Bra,B::Bra) = isequal(A.dims,B.dims)&&(A.data==B.data)
==(A::Operator,B::Operator) = isequal(A.dims,B.dims)&&(A.data==B.data)
isapprox(A::Ket,B::Ket) =  isequal(A.dims,B.dims)&&isapprox(A.data,B.data)
isapprox(A::Bra,B::Bra) =  isequal(A.dims,B.dims)&&isapprox(A.data,B.data)
isapprox(A::Operator,B::Operator) = isequal(A.dims,B.dims)&&isapprox(A.data,B.data)

# Show methods
function show{T<:QuMatrix}(io::IO, A::T)
    n = size(A,1)
    dim = join(A.dims,'⊗')
    print(io, "$n×$n $T with space dimensions $dim")
end
function show{T<:QuMatrix}(io::IO, ::MIME"text/plain", A::T)
    n = size(A,1)
    dim = join(A.dims,'⊗')
    println(io, "$n×$n $T with space dimensions $dim:")
    Base.showarray(io, A.data, false, header=false)
end
show(io::IO, ψ::QuVector) = print(io, braket(ψ))
function show(io::IO, ::MIME"text/plain", ψ::QuVector)
    n = size(ψ,1)
    dim = join(ψ.dims,'⊗')
    println(io, "$n-d $(typeof(ψ)) with space dimensions $dim:")
    println(io, braket(ψ))
end

# Bra-ket printing
function braket(ψ::QuVector, N::Int = 5)
    # N is the max number of kets/bras to show
    idx = find(x->abs2(x)>0.01,ψ.data)
    isempty(idx) && return "0"
    perm = sortperm(ψ.data[idx], by=abs2, rev=true)
    idx = length(perm)>N ? idx[perm[1:N]] : idx[perm]
    coeffs = map(full(ψ.data[idx])) do x
        @sprintf("%.2f∠%d°", abs(x), rad2deg(angle(x)))
    end
    labels = join.(reverse.(braketlabels.(idx.-1,[ψ.dims])),[","])
    return prettybraket(ψ,coeffs,labels)*(length(perm)>N?" +...":"")
end

prettybraket(::Ket,coeffs,labels) = join(string.(coeffs,["|"],labels,["⟩"])," + ")
prettybraket(::Bra,coeffs,labels) = join(string.(coeffs,["⟨"],labels,["|"])," + ")

function braketlabels{N}(n::Integer,bases::SDims{N})
    l = zeros(Int,N)
    for i in 1:N
        l[i] = rem(n, bases[N-i+1])
        n = div(n, bases[N-i+1])
    end
    return l
end
