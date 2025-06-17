#
# Retorna uma matriz com as posições dos centróides de
# cada elemento da malha
#
function Centroides(ne,connect,coord,elements_design)

    # Inicializa a matriz de saída
    C = zeros(ne,3)

    # A primeira coluna em coonect tem o tipo de elemento
    # 2 -> triangular (linear)
    # 3 -> quadrangular (linear)
    #
    #          3D
    # 4 -> tetrahedra (linear)
    # 5 -> hexaedra (linear)
    # 7 -> pyramid (linear)

    # Loop em cada elemento, identificado o tipo 
    for ele in elements_design

        # Tipo de elemento
        et = connect[ele,1]

        # Descobre nos, X e Y para este elemento
        nos, X = Nos_Coordenadas(ele,et,coord,connect) 

        # Média em cada uma das dimensões
        C[ele,:] .= mean(X,dims=1)'

    end #ele
    
    # Retorna a matriz com os centróides
    return C

end