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
function BESO(x::Vector{T1}, D::Vector{T1}, V::Vector{T1}, Vlim::Float64,  tol=1E-6, xmin=1E-3, xmax=1.0) where T1

    # Valores limites para a iteração 
    D_min = minimum(D)
    D_max = maximum(D)
 
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
        for i in LinearIndices(x)

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


