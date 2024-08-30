#
# Versão lendo direto do .msh (gmsh)
#
function Analise(meshfile::String;nev=4,Tf=1.0,Δt=1E-6,metodo=:Newmark,output=true)

    # Evita chamar um .geo
    if occursin(".geo",meshfile)
        error("Chamar com .msh..")
    end

    # Verifica se o método é válido
    metodo in [:Modal, :Bathe, :Newmark] || error("Métodos disponíveis são :Modal, :Bathe e :Newmark")

    # Le dados da malha
    nn, coord, ne, connect, materials, nodes_open, velocities, damping = Parsemsh_Daniele(meshfile)

    # Teste até arrumarmos os cálculos com os triângulos
    if any(connect[:,1].==2)
        println("################# CUIDADO ::: Elemento triangular está sendo implementado ")
        println("################# CUIDADO ::: resultados ainda não estão OK ")
    end

    # Calcula as matrizes globais
    K,M = Monta_KM(nn,ne,coord,connect,materials)

    # DOFs livres do problema
    livres = setdiff(collect(1:nn),nodes_open)
    
    # Inicializa um arquivo de pós-processamento do gmsh
    nome = "$(metodo).pos"
    etype = connect[:,1]
    Lgmsh_export_init(nome,nn,ne,coord,etype,connect[:,3:end])

    ###############################################################
    #                           Modal
    ###############################################################
    if metodo===:Modal

        println("Análise Modal com ", nev, " autovalores")

        # Chama a rotina de solução
        freq, X = Modal(K,M,livres,nev)
 
        # Exporta os modos para visualização no gmsh
        # Exporta os modos 1 até nev
        if output
            for  i in LinearIndices(freq)

                # Inicializa um vetor com todos os gls
                vv = zeros(nn)

                # Copia o modo para as posições livres
                vv[livres] = X[:,i]

                # Adiciona ao arquivo
                Lgmsh_export_nodal_scalar(nome,vv,"Modo $i, freq $(freq[i]) Hz")

            end # i 
        end

        return freq, X

    end # Modal

    ##############################################################
    #                       Transiente
    ##############################################################
    println("Análise transiente: usando o método $(metodo)")

    # E a de amortecimento
    C = Matriz_C(nn,damping,materials,coord,connect)

    # Faz a jogadinha para chamar os integradores no tempo
    P = zeros(nn)
    F(t) = Vetor_P!(t,nn,materials,velocities,coord,connect,P)

    # Chama o integrador
    if metodo===:Newmark
        tempos, MP = Newmark(M, C, K, F, livres, Δt, Tf)#,U0=U0)
    elseif metodo===:Bathe
        tempos, MP = B1B2Bathe2d(M, C, K, F, livres, Δt, Tf)#,U0=U0)
    end

    # Inicializa um vetor com todos os gls, pois se tivermos paredes 
    # abertas (Open), não calculamos estes gls
    vv = zeros(nn)

    # Exporta os tempos
    if output
        println("Escrevendo os dados transientes para o arquivo .pos")
        for  i=1:1:size(MP,2)

            # Copia o deslocamento para as posições livres
            vv .= MP[:,i]

            # Adiciona ao arquivo
            Lgmsh_export_nodal_scalar(nome,vv,"Pressão")

        end
    end

    return tempos, MP

end

