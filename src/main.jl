#
# Versão lendo direto do .msh (gmsh)
#
function Analise(meshfile::String,metodo=:Modal;nev=4,Tf=1.0,Δt=1E-6,freqs=[],U0=[],V0=[],output=true)

    # Evita chamar um .geo
    if occursin(".geo",meshfile)
        error("Chamar com .msh..")
    end

    # Verifica se o método é válido
    metodo in [:Modal, :Bathe, :Newmark, :Harmonic] || error("Métodos disponíveis são :Modal, :Bathe, :Newmark e :Harmonic")

    # Le dados da malha
    nn, coord, ne, connect, materials, nodes_open, velocities, damping, nodes_probe = Parsemsh_Daniele(meshfile)

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
    #                       Hamonic e Transiente
    ##############################################################
    
    # E a de amortecimento
    C = Matriz_C(nn,damping,materials,coord,connect)

    # Faz a jogadinha para chamar os integradores no tempo
    P = zeros(nn)
    F(t) = Vetor_P!(t,nn,materials,velocities,coord,connect,P)


    ##############################################################
    #                       Harmonic
    ##############################################################
    if metodo===:Harmonic

        # Verificamos se existem frequências sendo informadas
        if isempty(freqs)
            error("Analise Harmonica:: freqs deve ser um vetor não vazio")
        end

        # Número de frequências
        nω = length(freqs)
        
        # Número de nós a monitorar
        np = length(nodes_probe) 


        @show nodes_probe
        
        # Se não temos nós a monitorar, então monitoramos 
        # todos os livres
        if isempty(nodes_probe)
            nodes_probe = livres
            np = length(livres)
        end
        
        # Aloca matriz com os valores a serem monitorados
        monitor = zeros(ComplexF64,np,nω)

        # pre-aloca vetor de resposta
        U = zeros(ComplexF64, nn)

        # Loop pelas frequências
        contador = 1
        @showprogress "Harmonic... "  for f in freqs

            # Converte a freq para rad/s
            ω = 2*pi*f

            # Monta a matriz de rigidez dinâmica
            Kd = K[livres,livres] .+ im*ω*C[livres,livres] .- (ω^2)*M[livres,livres]

            # Monta o vetor de forças, que depende da frequência  
            Vetor_P!(0.0,nn,materials,velocities,coord,connect,P,ω=ω)

            # Soluciona 
            U[livres] .= Kd\P

            # Armazena os resultados na matriz de monitoramento
            monitor[:,contador] .= U[nodes_probe]

            # Incrementa o contador
            contador += 1

        end #f

        return monitor
        
    end


    ##############################################################
    #                       Transiente
    ##############################################################
    
    # Chama o integrador
    if metodo===:Newmark
        tempos, MP = Newmark(M, C, K, F, livres, Δt, Tf,U0=U0,V0=V0)
    elseif metodo===:Bathe
        tempos, MP = B1B2Bathe2d(M, C, K, F, livres, Δt, Tf,U0=U0,V0=V0)
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

