struct NoisyState{D} <: QuantumState{D}
    state::QuantumState{D}
    noise::ContinuousUnivariateDistribution

    function NoisyState(state::QuantumState{D}, noise::ContinuousUnivariateDistribution) where {D}
        new{D}(state, noise)
    end
end

NoisyState(state::QuantumState{D}) where {D} = NoisyState(state, Normal())
NoisyKet(v::Vector{T}, noise::ContinuousUnivariateDistribution) where {T <: Number} = NoisyState(Ket(v), noise)
NoisyKet(v::Vector{T}) where {T <: Number} = NoisyState(Ket(v))
NoisyKet(args::T...) where {T <: Number} = NoisyState(Ket(args...))
NoisyBra(v::Vector{T}, noise::ContinuousUnivariateDistribution) where {T <: Number} = NoisyState(Bra(v), noise)
NoisyBra(v::Vector{T}) where {T <: Number} = NoisyState(Bra(v))
NoisyBra(args::T...) where {T <: Number} = NoisyState(Bra(args...))

vec(ns::NoisyState) = vec(ns.state)
data(ns::NoisyState) = data(ns.state)
dims(::NoisyState{D}) where {D} = D

# Mutliplication of Quantum States







struct NoisyOperator{D} <: AbstractOperator{D}
    operator::AbstractOperator{D}
    noise::ContinuousUnivariateDistribution

    function NoisyOperator(operator::AbstractOperator{D}, noise::ContinuousUnivariateDistribution) where {D}
        new{D}(operator, noise)
    end
end

NoisyOperator(operator::AbstractOperator{D}) where {D} = NoisyOperator(operator, Normal())
