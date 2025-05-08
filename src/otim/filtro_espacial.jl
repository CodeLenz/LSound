#
# Filtro de vizinhança espacial
#
function Filtro(ne,vizinhos::Vector{Vector{T1}},pesos::Vector{Vector{T2}},vetor::Vector{T2}) where {T1,T2}

    # Aloca o vetor saída filtrada
    filtrado = zeros(T2,ne)

    # Loop em cada elemento
    for ele = 1:ne

        # Recupera os vizinhos e os pesos do elemento
        viz  = vizinhos[ele]
        peso = pesos[ele]

        # Parte de cima
        cima = sum( vetor[viz].*peso )

        # Parte de baixo
        baixo = sum(peso)

        # Vetor filtrado
        filtrado[ele] = cima/baixo

    end # ele

    # Retorna o vetor filtrado
    return filtrado

end