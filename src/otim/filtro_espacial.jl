#
# Filtro de vizinhança espacial
#
# filtra vetor, levando em consideração vizinhos e pesos, que foram 
# obtidos somente para os elementos de projeto da malha 
#
#
function Filtro(vizinhos::Vector{Vector{T1}},pesos::Vector{Vector{T2}},
                vetor::Vector{T2},elements_design::Vector) where {T1,T2}

    # Aloca o vetor saída filtrada
    filtrado = copy(vetor)

    # Loop em cada elemento
    for i in LinearIndices(vizinhos)

        # Recupera os vizinhos e os pesos do elemento
        viz  = vizinhos[i]
        peso = pesos[i]

        # Parte de cima
        cima = sum( vetor[viz].*peso )

        # Parte de baixo
        baixo = sum(peso)

        # Elemento que está sendo filtrado
        ele = elements_design[i] 

        # Vetor filtrado
        filtrado[ele] = cima/baixo

    end # ele

    # Retorna o vetor filtrado
    return filtrado

end