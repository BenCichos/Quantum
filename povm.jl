sicpovm_projectors = begin
    x = sqrt(2 + sqrt(5))
    sic_povm = [
        Ket([x, 1, 1, 1]),
        Ket([x, 1, -1, -1]),
        Ket([x, -1, 1, -1]),
        Ket([x, -1, -1, 1]),
        Ket([im, x, 1, -im]),
        Ket([im, x, -1, im]),
        Ket([-im, x, 1, im]),
        Ket([-im, x, -1, -im]),
        Ket([im, im, x, -1]),
        Ket([im, -im, x, 1]),
        Ket([-im, im, x, 1]),
        Ket([-im, -im, x, -1]),
        Ket([im, 1, -im, x]),
        Ket([im, -1, im, x]),
        Ket([-im, 1, im, x]),
        Ket([-im, -1, -im, x])]
    1 / sqrt(5 + sqrt(5)) .* sic_povm
end


sicpovm_matrices = map(Operator, sicpovm_projectors)