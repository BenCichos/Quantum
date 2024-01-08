abstract type Measurement{D} end

struct Projection{D} <: Measurement{D}
    ϕ::QuantumState{D}
    ρ̂::Operator{D}
    
    function Projection(ϕ::QuantumState{D}, ρ̂::Operator{D}) where {D}
        new{D}(ϕ, ρ̂)
    end
end

projector(m::Projection) = m.ϕ
operator(m::Projection) = m.ρ̂

dims(::Projection{D}) where {D} = D
hilbertspace(::Projection{D}) where {D} = D^2

measurement(ϕ::Ket{D}, ρ̂::Operator{D}) where {D} = Projection(ϕ, ρ̂)
measurement(ϕ::Bra{D}, ρ̂::Operator{D}) where {D} = Projection(ϕ, ρ̂)
measurement(ϕ::Ket{D}, ψ::Bra{D}) where {D} = Projection(ϕ, density(ψ))
measurement(ϕ::Bra{D}, ψ::Ket{D}) where {D} = Projection(ϕ, density(ψ))
measurement(ϕ::Ket{D}, ψ::Ket{D}) where {D} = Projection(ϕ, density(ψ))
measurement(ϕ::Bra{D}, ψ::Bra{D}) where {D} = Projection(ϕ, density(ψ))
evaluate(m::Projection) = expectation(projector(m), operator(m))


expectation(ϕ::Ket{D}, ρ̂::Operator{D}) where {D} = real(dag(ϕ) * ρ̂ * ϕ)
expectation(ϕ::Bra{D}, ρ̂::Operator{D}) where {D} = real(ϕ * ρ̂ * dag(ϕ))


fidelity(ρ̂::Operator{D}, σ̂::Operator{D}) where {D} = real(tr(ρ̂ * σ̂))
fidelity(ϕ::QuantumState{D}, ψ::QuantumState{D}) where {D} = real(tr(density(ϕ) * density(ψ)))


# Pauli Matrices Two Level Quantum Systems
PAULI0, PAULI1, PAULI2, PAULI3 = Operator([1 0; 0 1]), Operator([0 1; 1 0]), Operator([0 -im; im 0]), Operator([1 0; 0 -1])
PAULIMATRICES = [PAULI0, PAULI1, PAULI2, PAULI3]

# Pauli Matrices 2ᴺ Level Quantum Systems
function paulimatrix(index::Int; dims::Int=2)
    @assert ispow2(dims)
    dims == 2 && return PAULIMATRICES[index]
    (quotient, modulus) = fldmod1(index, dims)
    return kron(paulimatrix(quotient; dims=2), paulimatrix(modulus; dims=dims÷2))
end

paulimatrices(dims::Int) = map(i -> paulimatrix(i; dims=dims), 1:(dims^2))

paulimeasurement(ϕ::Ket) = map(pauli -> measurement(ϕ, pauli), paulimatrices(dims(ϕ)))
paulireconstruct(measurements::Vector{Projection{D}}) where {D} =  1/D * sum(evaluate.(measurements) .* paulimatrices(D))