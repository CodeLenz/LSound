

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

        # Monta as matrizes dos elementos
        Ke,Me,nos = KMe(ele,c,coord,connect)

        # Sobreposição das locais nas globais
        for i=1:4
            ni = nos[i]
            for j=1:4
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
