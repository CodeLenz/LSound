
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

            # Material for this element
            mat = connect[ele,2]

            # Find material density
            # ρ = materials[mat,1]
            # c = materials[mat,2]

            # Cte de amortecimento (passando o Y_n diretamente)
            damp = valor #ρ/valor

            # Find nodes and coordinates
            nos,X,Y = Nos_Coordenadas(ele,coord,connect)

            # Calcula a matriz local do elemento
            Ce = Damping_local(edge,damp,X,Y)

            # Loop para informar o amortecimento
            for i=1:4
                ni = nos[i]
                for j=1:4
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