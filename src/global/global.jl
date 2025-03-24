#
# Monta as matrizes globais K e M
#
function Monta_KM(nn,ne,coord,connect,materials)  
    #  function Monta_KM(nn,ne,coord,connect,materials,ρ = 1.204)
    # Qual é a posição na entrada de dados em que está o ρ  (como o c) ?    
    # A 101325 Pa e 20 °C = 293,15 K, ρ = 101325 / (287,058 ⋅ 293,15) = 1,204 kg/m³

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

        # Find material density
        ρ = materials[mat,1]

        # Velocidade
        c = materials[mat,2]
        
        # Calcula o módulo de compressibilidade κ
        # bulk modulus
        κ = ρ*c^2

        # Calcula as inversas 
        iρ = 1/ρ
        iκ = 1/κ    
        
        # Tipo de elemento
        et = connect[ele,1]

        # Descobre nos, X e Y para este elemento
        nos, X = Nos_Coordenadas(ele,et,coord,connect) 

        # Monta as matrizes dos elementos
        if et==3
           Ke, Me = KMe_bi4(iρ,iκ,X)
        elseif et==2
           Ke, Me = KMe_tri3(iρ,iκ,X)
        elseif et==4
            Ke, Me = KMe_tet4(iρ,iκ,X)   
        elseif et==5
           Ke, Me = KMe_hex8(iρ,iκ,X)
        elseif et==7
            Ke, Me = KMe_pyr5(iρ,iκ,X)    # 7  5-node pyramid.
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



function Matriz_C(nn,damping,coord,connect)

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

            # Cte de amortecimento (passando o Y_n diretamente)
            damp = 1/valor 

            # Find nodes and coordinates
            nos,X = Nos_Coordenadas(ele,et,coord,connect)

            # Calcula a matriz local do elemento
            if et==3
                Ce = Damping_local_bi4(edge,damp,X)
            elseif et==2
                Ce = Damping_local_tri3(edge,damp,X)
            elseif et==4
                Ce = Damping_local_tet4(edge,damp,X)
            elseif et==5
                Ce = Damping_local_hex8(edge,damp,X)
            elseif et==7
                Ce = Damping_local_pyr5(edge,damp,X)   # 7  5-node pyramid.
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