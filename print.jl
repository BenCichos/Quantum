function println(k::Ket)
    ket_data = data(k)
    string_rep = map(eachindex(ket_data)) do i
        zero_index = i - 1
        iszero(zero_index) ? "($(ket_data[i])) | $(string(i-1, base=2)) \\rangle " : " + ($(ket_data[i])) | $(string(i-1, base=2)) \\rangle "
    end
    display(latexstring(join(string_rep)))
end