abstract type QuantumState{D} end

struct Ket{D} <: QuantumState{D}
    data::Vector{ComplexF64}
    function Ket(v::Vector{T}) where {T <: Number}
        v = convert(Vector{ComplexF64}, v)
        new{length(v)}(v)
    end
end

Ket(args::T...) where {T <: Number} = Ket(collect(args))

struct Bra{D} <: QuantumState{D}
    data::Vector{ComplexF64}
    function Bra(v::Vector{T}) where {T <: Number}
        v = convert(Vector{ComplexF64}, v)
        new{length(v)}(v)
    end
end

Bra(args::T...) where {T <: Number} = Bra(collect(args))

abstract type AbstractOperator{D} end

struct Operator{D} <: AbstractOperator{D}
    data::Matrix{ComplexF64}
    function Operator(m::Matrix{T}) where {T <: Number}
        @assert size(m, 1) == size(m, 2)
        m = convert(Matrix{ComplexF64}, m)
        new{size(m, 1)}(m)
    end
end

Ket(d::Int) = Ket(zeros(ComplexF64, d))
Ket(b::Bra) = Ket(deepcopy(data(b)))
Ket(k::Ket) = Ket(deepcopy(data(k)))
Bra(d::Int) = Bra(zeros(ComplexF64, d))
Bra(k::Ket) = Bra(deepcopy(data(k)))
Bra(b::Bra) = Bra(deepcopy(data(b)))
Operator(d::Int) = Operator(zeros(ComplexF64, d, d))
Operator(op::Operator) = Operator(deepcopy(data(op)))
Operator(k::QuantumState) = density(k)

vec(ket::Ket) = data(ket)
vec(bra::Bra) = transpose(bra.data)
mat(op::Operator) = data(op)

data(qs::QuantumState) = qs.data
dims(::QuantumState{D}) where {D} = D
data(op::Operator) = op.data
dims(::AbstractOperator{D}) where {D} = D

# Mutliplication of Quantum States
*(k::Ket, b::Bra) = Operator(vec(k) * vec(b))
*(b::Bra, k::Ket) = vec(b) * vec(k)
*(op::Operator, k::Ket) = Ket(mat(op) * vec(k))
*(b::Bra, op::Operator) = Bra(transpose(vec(b) * mat(op)))
*(op1::Operator, op2::Operator) = Operator(mat(op1) * mat(op2))

# Addition of Quantum States
+(k1::Ket, k2::Ket) = Ket(data(k1) + data(k2))
+(b1::Bra, b2::Bra) = Bra(data(b1) + data(b2))
+(op1::Operator, op2::Operator) = Operator(data(op1) + data(op2))

# Subtraction of Quantum States
-(k1::Ket, k2::Ket) = Ket(data(k1) - data(k2))
-(b1::Bra, b2::Bra) = Bra(data(b1) - data(b2))
-(op1::Operator, op2::Operator) = Operator(data(op1) - data(op2))

# Scalar Multiplication of Quantum States
*(c::Number, k::Ket) = Ket(c * data(k))
*(k::Ket, c::Number) = Ket(c * data(k))
*(c::Number, b::Bra) = Bra(c * data(b))
*(b::Bra, c::Number) = Bra(c * data(b))
*(c::Number, op::Operator) = Operator(c * data(op))
*(op::Operator, c::Number) = Operator(c * data(op))

# Scalar Division of Quantum States
/(k::Ket, c::Number) = Ket(data(k) / c)
/(b::Bra, c::Number) = Bra(data(b) / c)
/(op::Operator, c::Number) = Operator(data(op) / c)

# Addition and Assignement of Quantum States 
+₌(k1::Ket, k2::Ket) = (data(k1) .+= data(k2); return k1)
+₌(b1::Bra, b2::Bra) = (data(b1) .+= data(b2); return b1)
+₌(op1::Operator, op2::Operator) = (data(op1) .+= data(op2); return op1)

# Subtraction and Assignement of Quantum States
-₌(k1::Ket, k2::Ket) = (data(k1) .-= data(k2); return k1)
-₌(b1::Bra, b2::Bra) = (data(b1) .-= data(b2); return b1)
-₌(op1::Operator, op2::Operator) = (data(op1) .-= data(op2); return op1)

# Scalar Multiplication and Assignement of Quantum States
*₌(k::Ket, c::Number) = (data(k) .*= c; return k)
*₌(c::Number, k::Ket) = (data(k) .*= c; return k)
*₌(b::Bra, c::Number) = (data(b) .*= c; return b)
*₌(c::Number, b::Bra) = (data(b) .*= c; return b)
*₌(op::Operator, k::Ket) = (data(k) .=  data(op) * data(k); return k)
*₌(b::Bra, op::Operator) = (data(b) .= data(b) * data(op); return b)
*₌(op1::Operator, op2::Operator) = (data(op1) .*= data(op2); return op1)


# Scalar Division and Assignement of Quantum States
/₌(k::Ket, c::Number) = (data(k) ./= c; return k)
/₌(b::Bra, c::Number) = (data(b) ./= c; return b)
/₌(op::Operator, c::Number) = (data(op) ./= c; return op)

# Normalize Quantum States
normalize(k::Ket) = Ket(normalize(data(k)))
normalize(b::Bra) = Bra(normalize(data(b)))
normalize(op::Operator) = Operator(data(op) / tr(data(op)))

normalize!(k::Ket) = (normalize!(data(k)); return k)
normalize!(b::Bra) = (normalize!(data(b)); return b)
normalize!(op::Operator) = (data(op) .= data(op) / tr(data(op)); return op)

# Outer Product
density(k::Ket) = kron(k, dag(k))
density(b::Bra) = kron(dag(b), b)

# Conjugate
conj(k::Ket) = Ket(conj(data(k)))
conj(b::Bra) = Bra(conj(data(b)))
conj(op::Operator) = Operator(conj(data(op)))

conj!(qs::QuantumState) = (data(qs) .= conj(data(qs)); return qs)
conj!(op::Operator) = (data(op) .= conj(data(op)); return op)

# Adjoint
adjoint(k::Ket) = Bra(conj(data(k)))
adjoint(b::Bra) = Ket(conj(data(b)))
adjoint(op::Operator) = Operator(Array(adjoint(data(op))))

adjoint!(op::Operator) = (data(op) .= adjoint(data(op)); return op)


# Dagger
dag(k::Ket) = adjoint(k)
dag(b::Bra) = adjoint(b)
dag(op::Operator) = adjoint(op)


dag!(op::Operator) = adjoint!(op)


# Kronecker Product
kron(k::Ket, b::Bra) = Operator(kron(vec(k), vec(b)))
kron(b::Bra, k::Ket) = kron(vec(b), vec(k))
kron(k1::Ket, k2::Ket) = Ket(kron(vec(k1), vec(k2)))
kron(b1::Bra, b2::Bra) = Bra(kron(vec(b1), vec(b2)))
kron(op1::Operator, op2::Operator) = Operator(kron(data(op1), data(op2)))


# Trace
tr(op::Operator) = tr(data(op))

# Hilbert Space
hilbertspace(::QuantumState{D}) where {D} = D^2
hilbertspace(::Operator{D}) where {D} = D^2

# Get Index
getindex(k::Ket, i::Int) = data(k)[i]
getindex(b::Bra, i::Int) = data(b)[i]

getindex(op::Operator, i::Int) = data(op)[i]
getindex(op::Operator, i::Int, j::Int) = data(op)[i, j]


# inner product
inner(k1::Ket, k2::Ket) = dag(k1) * k2 
innner(k::Ket) = inner(k, k)


# outer product
outer(k1::Ket{D}, k2::Ket{D}) where {D} = k1 * dag(k2)
outer(k::Ket) = density(k)