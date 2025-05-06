#
# Filtro de vizinhança espacial
#
function Filtro(ne,vizinhos,pesos,vetor)

    # Aloca o vetor saída filtrada
    filtrado = zeros(ne)

    # Loop em cada elemento
    for ele = 1:ne

        # Recupera os vizinhos e os pesos do elemento
        viz  = vizinhos[ele,:]
        peso = pesos[ele,:]

        # Parte de cima
        cima = sum( vetor[viz].*peso )

        # Parte de baixo
        baixo = sum(pesos)

        # Vetor filtrado
        filtrado[ele] = cima/baixo

    end # ele

    # Retorna o vetor filtrado
    return filtrado

end