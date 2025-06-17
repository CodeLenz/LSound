#
# Devolve o volume associado aos elementos cheios e o volume associado aos 
# elementos vazios
#
function Volume_cheios_vazios(x::Vector,V::Vector,elementos_projeto::Vector,xmin::Float64,xmax::Float64)

    # Volume cheio
    cheio = 0.0

    # Volume vazio
    vazio = 0.0

    # Loop pelos elementos de projeto
    for ele in elementos_projeto
        if x[ele] ≈ xmax 
           cheio += V[ele]
        elseif x[ele] ≈ xmin
           vazio += V[ele]
        end
    end

    # Retorna volume cheio e volume vazio
    return cheio, vazio

end


#
# Determina os valores de corte para remoção/adição
#
#
# TODO
#
# Como cada elemento pode ter um volume diferente, ne_adicionar
# e ne_remover deveriam ser calculados pelo número de elementos 
# que podemos acessar nos vetores decrecente e crescente, até 
# bater o volume em questão (adicionar/remover)
#
# Fazer o somatorio de volume no caso em cada situação para comparação
#
#
function Elementos_a_modificar(x::Vector, D::Vector,elementos_projeto::Vector,V::Vector,cheio::Float64,Vlim::Float64,xmin::Float64,xmax::Float64)

    # Números de sensibilidade somente dos elementos de projeto
    Dprojeto = D[elementos_projeto]

    # Valores de corte
    flag_adicao = false
    ne_adicionar = 0
    adicao  = 0.0

    flag_remocao = false
    ne_remover = 0
    remocao = 0.0

    # Listas de modificação  --> lista booleana ?  
    lista_adicao  = falses(length(x))
    lista_remocao = falses(length(x))

    #
    #                                           Adição
    #

    # Organiza do maior para o menor
    #
    ordem_decrescente = sortperm(Dprojeto,rev=true)

    # O volume a adicionar é 
    adicionar = Vlim-cheio

    println("Vlim-cheio = ", adicionar)

    # Se temos algo a adicionar (adicionar>0) então temos que ver quantos elementos 
    # podemos modificar para preencher essa diferenca (adicionar)
    if adicionar>0

        # Indica que precisamos adicionar material 
        flag_adicao = true

        # loop por todos os elementos de projeto, na ordem_decrescente
        # vendo se a variável de projeto do elemento é xmin, soma o volume  
        # e marca esse elemento para adição
        for ele in elementos_projeto[ordem_decrescente]   

            # Se o valor de D está acima do valor de corte para adição, adicionamos material 
            if x[ele] ≈ xmin

                # Adiciona o volume do elemento no total 
                adicao += V[ele]
                
                # Se passarmos do valor de volume a adicionar,
                # temos que sair do loop
                if adicao>=adicionar
                    break
                end

                # Número de elementos a adicionar
                ne_adicionar += 1

                # Marca o elemento para adição
                lista_adicao[ele] = true 

              
            end # if xmin   

        end # ele

    end # if adicionar > 0 

    
    #             Remoção

    # Organiza do  menor para o maior
    #
    ordem_crescente = sortperm(Dprojeto)

    # O volume a remover é o complementar ao da adição
    remover = cheio - Vlim

    println("cheio - Vlim = ", remover)

    # Se temos algo a remover (remover>0) então temos que ver quantos elementos 
    # podemos modificar para preencher essa diferenca (remover)
    if remover>0

        # Indica que precisamos remover material 
        flag_remocao = true

        # loop por todos os elementos de projeto, na ordem_crescente
        # vendo se a variável de projeto do elemento é xmax, soma o volume  
        # e marca esse elemento para remocao
        for ele in elementos_projeto[ordem_crescente]  

            # Se o elemento tem xmax, podemos marcar para passar para xmin
            if x[ele] ≈ xmax

                # Acrescenta o volume do elemento no total para remoção
                remocao += V[ele]

                # Se passarmos do valor de volume a remover,
                # temos que sair do loop
                if remocao>=remover
                    break
                end

                # Marca o elemento para remoação
                lista_remocao[ele] = true 

                # Número de elementos a remover
                ne_remover += 1
                
            end #if xmax   
            
        end #ele

    end # if remover

    @show ne_adicionar, ne_remover

    # Retorna o valor de corte para adição/remoção
    return lista_adicao, lista_remocao

end

#
# Mais uma tentativa de BESO
#
function BESO3(x::Vector{T1}, D::Vector{T1}, V::Vector{T1}, Vlim::Float64, elements_design::Vector; xmin=1E-3, xmax=0.99) where T1

    # Copia o vetor de variáveis de projeto para um vetor de trabalho
    xn = copy(x)

    # Primeiro calculamos os volumes de cheio de de vazio
    cheio,_  = Volume_cheios_vazios(xn,V,elements_design,xmin,xmax)
    println("Volume cheio inicial: ", cheio)

    # Agora verificamos quais são os valores de corte para adição e remoção de material 
    lista_adicao, lista_remocao = Elementos_a_modificar(x,D,elements_design,V,cheio,Vlim, xmin, xmax)

    # println("Flag de adição: ", flag_adicao, " Adição: ", adicao)
    # println("Flag de remoção: ", flag_remocao, " Remoção: ", remocao)
    # println("Lista de adição: ", lista_adicao)
    # println("Lista de remoção: ", lista_remocao)

    # Muda as densidades relativas dos elemento 
    xn[lista_adicao]   .= xmax
    xn[lista_remocao]  .= xmin

    # retorna o novo vetor 
    return xn

end