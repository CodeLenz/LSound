#
# BESO as in Zhi Hao Zuo, Yi Min Xie, 
# A simple and compact Python code for complex 3D topology optimization, 2015
#
#
# x    -> vetor de variáveis de projeto
# D    -> vetor com os sensitivity indexes
# Vlim -> Volume limite para esta iteração do BESO
# V    -> Vetor com o volume (não parametrizado) de cada elemento
# tol  -> critério de parada do método da biseção
# xmin -> valor mínimo ("vazio") 
# xmax -> valor máximo ("cheio")
#
#
function BESO(x::Vector{T1}, D::Vector{T1}, V::Vector{T1}, Vlim::Float64, elements_design::Vector, tol=1E-6, xmin=1E-3, xmax=0.99) where T1

    # Valores limites para a iteração 
    D_min = minimum(D[elements_design])
    D_max = maximum(D[elements_design])
 
    # Copia o vetor atual para um novo vetor
    xn = copy(x)

    # Inicializa o volume fora do loop para podermos 
    # recuperar depois
    volume = 0.0

    # Contador de iterações do loop de biseção
    # Se o contador for zero, então não fizemos nenhuma modificação 
    # na malha
    contador = 0

    # Loop de atualização, mantendo a restrição de volume
    # vamos usar um método de biseção, como no critério de
    # ótimo
    while (D_max - D_min)/D_max > tol

        # Atualiza o contador 
        contador += 1

        # Valor de corte, interpolando pelos valores limites
        corte = (D_min + D_max)/2 

        # Inicializa o volume atual da estrutura 
        volume = 0.0

        # Para cada variável de projeto
        for i in elements_design #LinearIndices(x)

            # Se o indicador for abaixo do valor de 
            # referência, tira material do elemento i
            if D[i] < corte

                xn[i] = xmin  

            # Do contrário, adiciona material no elemento i 
            elseif D[i] > corte

                xn[i] = xmax 

            end #if

            # Soma o volume
            volume += xn[i]*V[i]
        
        end #i

        # Se o volume atual estiver acima do valor limite
        # movemos o corte para o D_min, do contrário para 
        # o D_max
        if volume > Vlim
            D_min = corte
        else
            D_max = corte
        end

    end # while
    
    # Retorna o novo vetor de variáveis de projeto e o contador
    return xn, contador

end


#
# Essa rotina está acessando elementos que não são de projeto 
# e também está deixando elementos com densidades intermediárias <- problemão
#
function BESO2(ne::Int64, x::Vector{T1}, D::Vector{T1}, V::Vector{T1}, Vlim::Float64, elements_design::Vector, ar=0.05; tol=1E-6, xmin=1E-3, xmax=0.99) where T1

    # Devolve uma lista com os elementos 
    # em ordem crescente de índices de sensibilidade (menor para o maior)
    eles_D_sort = sortperm(D)

    # Precisamos pegar somente os que são de projeto 
    # mas mantendo a sequência numérica dos elementos
    eles_D_sortp = indexin(elements_design,eles_D_sort)
 
    # Número de elementos em eles_D_sortp
    nep = length(eles_D_sortp)

    # Copia o vetor atual para um novo vetor
    xn = copy(x)

    # Inicializa o volume fora do loop para podermos 
    # recuperar depois
    volume = 0.0

    # Contador de iterações do loop de biseção
    # Se o contador for zero, então não fizemos nenhuma modificação 
    # na malha
    contador = 1

    # O número de elementos adicionados/removidos deve ser igual 
    # ou próximo a
    ne_mod = ceil(Int64,nep*ar)

    # Meio da lista
    meio = floor(Int64,nep/2)

    @show ne_mod, meio

    # Numero efetivo de elementos adicionados/removidos
    AR = 0

    # Loop pelas extremidades de eles_D_sort
    for i=1:floor(Int64,ne_mod/2)

        # Para os elementos que estão do meio para o 
        # começo da lista modificamos para vazio
        if xn[eles_D_sortp[meio+i]]!=xmin
           AR += 1
        end 
        xn[eles_D_sortp[meio+i]] = xmin

        # Os que estão no do meio para o final da lista
        # passam para cheio
        if xn[eles_D_sortp[meio-i]]!=xmax
            AR += 1
        end  
        xn[eles_D_sortp[meio-i]] = xmax

    end

    # Volume atualizado
    volume = sum(V.*xn)

    @show AR, AR/ne
    @show volume, Vlim
    
    # Se o volume estiver acima do volume target, podemos 
    # corrigir alterando as variáveis de projeto 
    
    # Retorna o novo vetor de variáveis de projeto e o contador
    return xn, contador

end



