

#
# Monta as matrizes globais K e M
#
function Monta_KM(nn,ne,coord,connect,materials)

    # Aloca as matrizes globais
    K = spzeros(nn,nn)
    M = spzeros(nn,nn)

    # Loop pelos elementos
    for ele=1:ne

        # Material
        mat = connect[ele,2] 

        # Velocidade
        c = materials[mat,2]

        # Monta as matrizes dos elementos
        Ke,Me,nos = KMe(ele,c,coord,connect)

        # Sobreposição das locais nas globais
        K[nos,nos] += Ke
        M[nos,nos] += Me

    end

    return K, M

end
