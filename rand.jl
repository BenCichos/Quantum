rand(rng::AbstractRNG, ::SamplerType{Ket{N}}) where {N} = Ket(normalize!(rand(rng, ComplexF64, N)))
rand(rng::AbstractRNG, ::SamplerType{Bra{N}}) where {N} = Bra(normalize!(rand(rng, ComplexF64, N)))

randket(d::Int) = rand(Ket{d})
randbra(d::Int) = rand(Bra{d})