#
# Retorna um vetor com o volume (ou área se 2D) de cada elemento 
# da malha
#
function Volumes(ne,connect,coord)

    # Inicializa o vetor de saída
    V = zeros(ne)

    # A primeira coluna em coonect tem o tipo de elemento
    # 2 -> triangular (linear)
    # 3 -> quadrangular (linear)
    #
    #          3D
    # 4 -> Tetrahedra (linear)
    # 5 -> hexaedra (linear)
    # 7 -> pyramid (linear)

    # Loop em cada elemento, identificado o tipo e chamando a função 
    # correta
    for ele = 1:ne 

        # Tipo de elemento
        et = connect[ele,1]

        # Descobre nos, X e Y para este elemento
        nos, X = Nos_Coordenadas(ele,et,coord,connect) 

        # Calcula a área ou o volume de cada elemento
        if et==3
            V[ele] = Area_bi4(X)
        elseif et==2
            V[ele] = Area_tri3(X)
        elseif et==4
            V[ele] = Volume_tet4(X)
        elseif et==5
            V[ele] = Volume_hex8(X)
        elseif et==7
            V[ele] = Volume_pyr5(X)    
        else
            error("Volumes::Elemento não definido")
        end

    end #ele
    
    # Retorna o vetor com os volumes de cada um dos elementos da malha
    return V

end

#
# Rotina que calcula o volume da próxima iteração, baseado no conceito de 
# er (evolutionary rate, ou taxa de evolução)
#
function Calcula_volume_er(er,volume_atual,Vast)

    # Volume a ser utilizado como limite para esta iteração
    # Este é o principal parâmetro de controle para a estabilidade
    # do processo de otimização
    if volume_atual > Vast
        vol = max(Vast, volume_atual*(1.0 - er))
    elseif volume_atual < Vast
        vol = min(Vast,volume_atual*(1 + er))
    else
        vol = Vast
    end 
 
    # Caso tenhamos um volume nulo, decorrente de uma 
    # inicialização com todos os γ=γ_min, devemos utilizar
    # vol como sendo algum valor pré-determinado, pois o 
    # if das linhas anteriores vai retornar zero
    if vol==0
       vol = er*Vast
    end

    # Retorna o volume a ser atendido pelo BESO
    return vol

end 