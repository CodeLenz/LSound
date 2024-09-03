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



function Matriz_C(nn,damping,materials,coord,connect)

    # Aloca os vetores para a matriz de amortecimento
    I = Int64[]
    J = Int64[]
    VC = Float64[]

    # Loop pelo dicionário de damping
    for d in damping

        # Recover data from Dictionary
        valor    = d["value"]
        elements = d["elements"] 

        # Agora precisamos fazer um loop sobre os elementos
        # e suas arestas
        for e in axes(elements,1)

            # Element and edge
            ele  = elements[e,1]
            edge = elements[e,2]

            # Tipo de elemento
            et = connect[ele,1]

            # Material for this element
            mat = connect[ele,2]

            # Find material density
            # ρ = materials[mat,1]
            # c = materials[mat,2]

            # Cte de amortecimento (passando o Y_n diretamente)
            damp = valor #ρ/valor

            # Find nodes and coordinates
            nos,X,Y = Nos_Coordenadas(ele,et,coord,connect)

            # Calcula a matriz local do elemento
            if et==3
                Ce = Damping_local_bi4(edge,damp,X,Y)
            elseif et==2
                Ce = Damping_local_tri3(edge,damp,X,Y)
            else
                error("Tipo de elemento não definido")
            end

            # Loop para informar o amortecimento
            for i in LinearIndices(nos)
                ni = nos[i]
                for j in LinearIndices(nos)
                    nj = nos[j]
                    push!(I,ni)
                    push!(J,nj)
                    push!(VC,Ce[i,j])    
                end #j
            end #i
        end #e
    end # dict (d)

   # Devolve a matriz esparsa de amortecimento
   return sparse(I,J,VC,nn,nn)

 end