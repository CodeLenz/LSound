#
# Retorna um vetor de vetores, com os vizinhos de cada elemento 
#
function Vizinhanca(ne,centroides,raio_filtro,elements_design::Vector)

    # Aloca a lista de saída para os vizinhos
    vizinhos = Vector{Vector{Int64}}(undef,ne)

    # Aloca a lista de pesos
    pesos = Vector{Vector{Float64}}(undef,ne)

    # Loop pelos elementos da malha
    for ele in elements_design 

        # Centroide deste elemento
        cele = centroides[ele,:]

        # Aloca um vetor para armazenarmos os vizinhos deste elemento
        vele = Int64[]

        # Aloca um vetor para armazenarmos os pesos
        pele = Float64[]

        # Loop por todos os elementos da malha
        for viz = 1:ne

            # Centroide do vizinho
            cviz = centroides[viz,:]

            # Distância entre os centróides
            dist = norm(cviz.-cele)

            # Se essa distância for menor do que o raio_filtro,
            # armazena como vizinho
            if dist < raio_filtro
               push!(vele, viz)
               push!(pele, 1 - dist/raio_filtro)
            end

        end # viz

        # Armazena vele na linha ele de vizinhos
        vizinhos[ele] = copy(vele)

        # Armazena pele na linha ele de pesos
        pesos[ele] = copy(pele)

    end # ele

    # Retorna o vetor de vetores com os vizinhos do elemento
    # e também os pesos
    return vizinhos, pesos

end