

#
# Monta as matrizes globais K e M
#
function Monta_KM(nn,ne,coord,connect,materials)

    # Aloca vetores para a montagem eficiente 
    # das matrizes esparsas
    I = Int64[]
    J = Int64[]
    VK = Float64[]
    VM = Float64[]
 
    # Loop pelos elementos
    for ele=1:ne

        # Material
        mat = connect[ele,2] 

        # Velocidade
        c = materials[mat,2]

        # Tipo de elemento
        et = connect[ele,1]

        # Descobre nos, X e Y para este elemento
        nos, X, Y = Nos_Coordenadas(ele,et,coord,connect) 

        # Monta as matrizes dos elementos
        if et==3
           Ke, Me = KMe_bi4(ele,c,X,Y)
        elseif et==2
           Ke, Me = KMe_tri3(ele,c,X,Y)
        else
            error("Elemento não definido")
        end
 
        # Sobreposição das locais nas globais
        for i in LinearIndices(nos)
            ni = nos[i]
            for j in LinearIndices(nos)
                push!(I,ni)
                push!(J,nos[j])
                push!(VK,Ke[i,j])
                push!(VM,Me[i,j])
            end
        end

    end

    # Retorna as matrizes globais
    return sparse(I,J,VK), sparse(I,J,VM)

end
