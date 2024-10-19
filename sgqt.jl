# Self-Guided Quantum Tomography

# Algorithm Parameters
β(i::Int; b::Float64, t::Float64=1/6) =  b / ((i + 1)^t)
γ(i::Int; a::Float64, A::Float64=0.0, s::Float64=1.0) = a / ((i + 1 + A)^s)


# Direct Fidelity Estimation
function estimate_fidelity(ρ̂::Operator{4}, σ̂::Operator{4}, Sᵢ::Vector{Operator{4}}) 
    return 1/length(Sᵢ) * sum(ϵ -> fidelity(ρ̂, ϵ)*fidelity(ϵ, σ̂), Sᵢ)
end

# Self-Guided Quantum Tomography - Two Level Quantum System
function sgqt_walk(ρ̂::Operator{2}, ϕ::Ket{2}, i::Int, paramsᵦ::NamedTuple, paramsᵧ::NamedTuple)
    Δᵢ = rand(Ket{2})
    βᵢ = β(i; paramsᵦ...)
    γᵢ = γ(i; paramsᵧ...)

    ϕ₊, ϕ₋ = normalize!(ϕ + (Δᵢ*βᵢ)), normalize!(ϕ - (Δᵢ*βᵢ))
    E₊, E₋ = evaluate(measurement(ϕ₊, ρ̂)), evaluate(measurement(ϕ₋, ρ̂))
    gᵢ = (E₊ - E₋) / 2βᵢ
    return normalize(ϕ + γᵢ*gᵢ*Δᵢ)
end

function sgqt(ρ̂::Operator{2}, iterations::Int; paramsᵦ::NamedTuple=(;b=0.1, t=1/6), paramsᵧ::NamedTuple=(;a=3.0, A=0.0, s=1.0))
    ϕᵢ = rand(Ket{2})
    ψs = Ket{2}[]
    sizehint!(ψs, iterations)
    for i in 1:iterations
        push!(ψs, ϕᵢ)
        ϕᵢ = sgqt_walk(ρ̂, ϕᵢ, i, paramsᵦ, paramsᵧ)
    end
    ψs
end


sgqt(ϕ::Ket{2}, iterations::Int; paramsᵦ::NamedTuple=(;b=0.1, t=1/6), paramsᵧ::NamedTuple=(;a=3.0, A=0.0, s=1.0)) = sgqt(density(ϕ), iterations, paramsᵦ=paramsᵦ, paramsᵧ=paramsᵧ)

function sgqt_walk(ρ̂::Operator{4}, ϕ::Ket{4}, measures::AbstractVector, cardinality::Int, i::Int, paramsᵦ::NamedTuple, paramsᵧ::NamedTuple, weights::Vector{Float64})
    Δᵢ = rand(Ket{4})
    Sᵢ = sample(measures, Weights(weights), cardinality, replace=false)
    βᵢ = β(i; paramsᵦ...)
    γᵢ = γ(i; paramsᵧ...)

    σ̂₊, σ̂₋ = normalize!(density(ϕ + (Δᵢ*βᵢ))), normalize!(density(ϕ - (Δᵢ*βᵢ)))
    F̃₊, F̃₋ = estimate_fidelity(ρ̂, σ̂₊, Sᵢ), estimate_fidelity(ρ̂, σ̂₋, Sᵢ)
    gᵢ = (F̃₊ - F̃₋)/2βᵢ
    return normalize(ϕ + γᵢ*gᵢ*Δᵢ)
end

function sgqt(ρ̂::Operator{4}, iterations::Int; measures::AbstractVector=paulimatrices(4), cardinality::Int=8, paramsᵦ::NamedTuple=(;b=0.1, t=1/6), paramsᵧ::NamedTuple=(;a=3.0, A=0.0, s=1.0), weights::Vector{Float64}=fill(1/16, 16))
    ϕᵢ = rand(Ket{4})
    ψs = Ket{4}[]
    sizehint!(ψs, iterations)
    for i in 1:iterations
        push!(ψs, ϕᵢ)
        ϕᵢ = sgqt_walk(ρ̂, ϕᵢ, measures, cardinality, i, paramsᵦ, paramsᵧ, weights)
    end
    ψs
end

sgqt(ϕ::Ket{4}, iterations::Int; measures::AbstractVector=paulimatrices(4), cardinality::Int=8, paramsᵦ::NamedTuple=(b=0.1, t=1/6), paramsᵧ::NamedTuple=(a=3.0, A=0.0, s=1.0), weights::Vector{Float64}=fill(1/16, 16)) = sgqt(density(ϕ), iterations; measures=measures, cardinality=cardinality, paramsᵦ=paramsᵦ, paramsᵧ=paramsᵧ, weights=weights)