
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
#
function Valores_corte(D::Vector,elementos_projeto::Vector,V::Vector,cheio::Float64,Vlim::Float64)

    # Uma roubadinha marota aqui...volume médio dos elementos
    Vmedio = mean(V[elementos_projeto])

    # Números de sensibilidade somente dos elementos de projeto
    Dprojeto = D[elementos_projeto]

    # Valores de corte
    flag_adicao = false
    adicao  = 0.0

    flag_remocao = false
    remocao = 0.0

    #             Adição

    # Organiza do maior para o menor
    #
    # ordem_decrecente = sortperm(Dprojeto,rev=true)
    #
    # toda a vez que quiser acessar Dprojeto, pode usar
    # a lista ordem_decrecente direto. Mas para acessar, 
    # por exemplo, V tem que usar elementos_projeto[alguma posição em ordem crescente]
    #
    # decrecente = Dprojeto[ordem_decrecente]
    #
    decrecente = sort(Dprojeto, rev=true)

    # O volume a adicionar é 
    adicionar = Vlim-cheio

    # Se temos algo a adicionar (adicionar>0) então temos que ver quantos elementos 
    # podemos modificar para preencher essa diferenca (adicionar)
    if adicionar>0

        # Indica que precisamos adicionar material 
        flag_adicao = true

        # loop por todos os elementos de projeto, na ordem_decrecente
        # vendo se a variável de projeto do elemento é xmin, soma o volume  
        # e marca esse elemento para adição

        # O número de elementos associado a essa diferença será 
        # arredonda para cima para evitarmos que não seja removido
        # nenhum elemento (por questão de tamanho da malha)
        ne_adicionar = ceil(Int64,adicionar/Vmedio)

        # Valor de corte para adição seria 
        adicao = decrecente[ne_adicionar]

    end

    
    #             Remoção

    # Organiza do  menor para o maior
    #
    # ordem_crescente = sortperm(Dprojeto)
    #
    # toda a vez que quiser acessar Dprojeto, pode usar
    # a lista ordem_crescente direto. Mas para acessar, 
    # por exemplo, V tem que usar elementos_projeto[alguma posição em ordem crescente]
    #
    # crescente = Dprojeto[ordem_crescente]
    #

    # Organiza do menor para o maior
    crescente = sort(Dprojeto)

    # O volume a remover é o complementar ao da adição
    remover = cheio - Vlim

    # Se temos algo a remover (remover>0) então temos que ver quantos elementos 
    # podemos modificar para preencher essa diferenca (remover)
    if remover>0

        # Indica que precisamos remover material 
        flag_remocao = true

        # loop por todos os elementos de projeto, na ordem_crescente
        # vendo se a variável de projeto do elemento é xmax, soma o volume  
        # e marca esse elemento para remocao


        # O número de elementos associado a essa diferença será 
        # arredonda para cima para evitarmos que não seja adicionaro
        # nenhum elemento (por questão de tamanho da malha)
        ne_remover = ceil(Int64,remover/Vmedio)

        # Valor de corte para remocao seria 
        remocao = crescente[ne_remover]

    end

    # Retorna o valor de corte para adição/remoção
    return flag_adicao, adicao, flag_remocao, remocao

end

#
# Mais uma tentativa de BESO
#
function BESO3(x::Vector{T1}, D::Vector{T1}, V::Vector{T1}, Vlim::Float64, elements_design::Vector; xmin=1E-3, xmax=0.99) where T1

    # Copia o vetor de variáveis de projeto para um vetor de trabalho
    xn = copy(x)

    # Primeiro calculamos os volumes de cheio de de vazio
    cheio, _  = Volume_cheios_vazios(x,V,elements_design,xmin,xmax)

    # Agora verificamos quais são os valores de corte para adição e remoção de material 
    flag_adicao, adicao, flag_remocao, remocao = Valores_corte(D,elements_design,V,cheio,Vlim)

    # Acima, já recebemos as listas de modificação
    # xn[lista_adicao]  .= xmax
    # xn[lista_remocao] .= xmin 
      
    # E, com isso, podemos iterar nos elementos de projeto, modificando a variável de projeto
    for ele in elements_design

        # Se o valor de D está abaixo do valor de corte para remoção, removemos
        if flag_remocao && D[ele]<remocao && x[ele] ≈ xmax
            xn[ele] = xmin
        end    

        # Se o valor de D está acima do valor de corte para adição, socamos material 
        if flag_adicao && D[ele]>adicao && x[ele] ≈ xmin
            xn[ele] = xmax
        end    

    end #ele

    # retorna o novo vetor 
    return xn

end