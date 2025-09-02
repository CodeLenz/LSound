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

#
# Find all elements that share edges with the current element
#
function NeighborEdges(ne,connect,elements_design)

    # Vector of vectors
    vizinhos = Vector{Vector{Int64}}(undef,ne)

    # Local vector 
    local_vector = Int64[]

    # Loop pelos elementos de projeto
    for ele in elements_design

        # Limpa o vetor local 
        empty!(local_vector)

        # Nós deste elemento 
        nos_ele = connect[ele,3:end]

        # Tipo de elemento 
        etype = connect[ele,1]

        # Dependendo do tipo de elemento, temos diferentes requisitos para 
        # comparação do número de nós que são necessários para definir uma 
        # face. Não é tão simples se considerarmos a pirâmide de 5 nós, mas 
        # para os outros podemos fazer direto 
        #
        # Válido para triângulo e para quadrado (etypes 2 e 3)
        n_lado = 2
        if etype==4
            n_lado  = 3
        elseif  etype==5
            n_lado = 4
        elseif etype>5
            error("fix this code...") 
        end

        # Número de nós do elemento -2 (para comparação)
        ncompara = length(nos_ele)-n_lado

        # Loop pelos outros elementos de projeto 
        for viz in elements_design

            # pulamos o próprio elemento 
            if ele!=viz

               # Nós deste elemento 
               nos_viz = connect[viz,3:end]

               # Verifica se temos ao menos 2 nós em comum. Vou fazer isso com o
               # comando setfiff(a,b), que retorna todos os elementos que estão em a 
               # mas não estão em b. Se 
               compara = setdiff(nos_ele,nos_viz)

               # Se tivemos (ao menos) menos 2 nós em compara
               # então eles são em comum com o candidato a vizinho
               # Guarda para vizinho de ele 
               if length(compara)<=ncompara
                  push!(local_vector,viz)
               end

            end


        end #viz

        # Armazena no vetor de vetores com os vizinhos
        vizinhos[ele] = copy(local_vector)
        
    end #ele

    # Retorna 
    return vizinhos

end