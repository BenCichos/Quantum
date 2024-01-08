rand(::Type{Ket{N}}) where {N} = Ket(normalize!(rand(ComplexF64, N)))
rand(::Type{Bra{N}}) where {N} = Bra(normalize!(rand(ComplexF64, N)))
randket(d::Int) = Ket(normalize!(rand(ComplexF64, d)))
randbra(d::Int) = Bra(normalize!(rand(ComplexF64, d)))