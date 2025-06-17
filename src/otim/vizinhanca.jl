#
# Retorna um vetor de vetores, com os vizinhos de cada elemento 
# considerando SOMENTE os elementos de projeto
#
function Vizinhanca(ne,centroides,raio_filtro,elements_design::Vector)

    # Número de variáveis de projeto
    np = length(elements_design)

    # Aloca a lista de saída para os vizinhos
    vizinhos = Vector{Vector{Int64}}(undef,np)

    # Aloca a lista de pesos
    pesos = Vector{Vector{Float64}}(undef,np)

    # Número mínimo e máximo de vizinhos de 
    # um elemento da malha
    n_min_viz = np
    n_max_viz = 0

    # Loop pelos elementos da malha
    contador = 1
    for ele in elements_design 

        # Centroide deste elemento
        cele = centroides[ele,:]

        # Aloca um vetor para armazenarmos os vizinhos deste elemento
        vele = Int64[]

        # Aloca um vetor para armazenarmos os pesos
        pele = Float64[]

        # Loop por todos os elementos da malha
        for viz in elements_design

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

        # Número de vizinhos deste elemento 
        nviz = length(vele)

        # Verifica número de vizinhos
        n_min_viz = min(n_min_viz,nviz) 
        n_max_viz = max(n_max_viz,nviz)

        # Armazena vele na linha ele de vizinhos
        vizinhos[contador] = copy(vele)

        # Armazena pele na linha ele de pesos
        pesos[contador] = copy(pele)

        # Acumula o contador
        contador += 1

    end # ele

    # Mostra o número mínimo e máximo de vizinhos
    println("Número mínimo de vizinhos na malha: ", n_min_viz)
    println("Número máximo de vizinhos na malha: ", n_max_viz)

    # Retorna o vetor de vetores com os vizinhos do elemento
    # e também os pesos. Esses vetores de vetores 
    # contém somente os elementos de projeto
    return vizinhos, pesos

end