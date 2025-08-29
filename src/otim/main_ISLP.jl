#
# TODO
#
# raio_filtro -> pensar em uma maneira de entrar com este dado
#

#
# Rotina principal
#
"""
 Otim(meshfile::String,freqs=[])

 Basic input:

 meshfile -> arquivo de entrada (.msh)

 freqs -> vetor com as frquências (em Hz)

 vA -> Vetor com as sensibilidades para cada frequência

 verifica_derivada -> bool 

 Inputs -> freqs, a vector with the frequencies to sweep

"""
function Otim_ISLP(meshfile::String,freqs::Vector, vA::Vector;verifica_derivada=false)
    
    # Evita chamar um .geo
    occursin(".geo",meshfile) && error("Chamar com .msh..")
    
    # Define os nomes dos arquivos de entrada (yaml) e de saída
    # (pos) em função do nome de meshfile
    arquivo_yaml = meshfile[1:end-3]*"yaml"
    arquivo_pos  = meshfile[1:end-3]*"pos"

    # Define o nome dos arquivos contendo a distribuição inicial e final (otimizada)
    # das variáveis de projeto. Esses arquivos serão utilizados posteriormente em 
    # Processa_FRF
    arquivo_γ_ini = meshfile[1:end-3]*"_γ_ini.dat"
    arquivo_γ_fin = meshfile[1:end-3]*"_γ_opt.dat"

    # Verificamos se existem frequências sendo informadas
    isempty(freqs) && error("Analise Harmonica:: freqs deve ser um vetor não vazio")

    # Evita passar as frequências como algo diferente de um Vetor de floats
    isa(freqs,Vector{Float64}) || error("freqs deve ser um vetor de floats")
    
    # Verifica se os arquivos de entrada existem
    isfile(meshfile) || error("Otim:: arquivo de entrada $meshfile não existe")

    # Número de frequências
    nf = length(freqs)

    # Verifica se A foi definido
    if isempty(vA)
       vA = ones(nf)
    else 
       length(vA)==nf || error("Otim:: dimensão de vA deve ser nf")
    end

    # Arquivo .yaml
    isfile(arquivo_yaml) || error("Otim:: arquivo de entrada $(arquivo_yaml) não existe")

    # Le dados da malha
    nn, coord, ne, connect, materials, nodes_open, velocities, nodes_pressure, pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed = Parsemsh_Daniele(meshfile)

    # Lista com os elementos que são de projeto
    elements_design = setdiff(1:ne,sort!(elements_fixed))

    # Le os dados do arquivo yaml
    raio_filtro, niter, nhisto, ϵ1, ϵ2, vf, parametrizacao, γ_min, γ_max, partida = Le_YAML(arquivo_yaml)

    # Seleciona as rotinas de parametrização de material de acordo com 
    # a opção 
    @show parametrizacao

    # Vetores de funções 
    vetor_fρ  = [fρ_pereira, fρ_duhring]
    vetor_dfρ = [dfρ_pereira, dfρ_duhring]
    vetor_fκ  = [fκ_pereira, fκ_duhring]
    vetor_dfκ = [dfκ_pereira, dfκ_duhring]

    # 1 para PEREIRA e 2 para Duhring
    ponteiro_parametrizacao = 2
 
    if parametrizacao=="PEREIRA"
         println("Utilizando a parametrização de PEREIRA")
         #fρ(γ)  = fρ_pereira(γ) #,ψ, ρ_ar = ρ_ar, ρ2 = ρ_solido)
         #dfρ(γ) = dfρ_pereira(γ)
         #fκ(γ)  = fκ_pereira(γ)
         #dfκ(γ) = dfκ_pereira(γ)

    elseif parametrizacao=="DUHRING"
         println("Utilizando a parametrização de DUHRING")
         ponteiro_parametrizacao = 2

         #fρ(γ)  = fρ_duhring(γ)
         #dfρ(γ) = dfρ_duhring(γ)
         #fκ(γ)  = fκ_duhring(γ)
         #dfκ(γ) = dfκ_duhring(γ)
    end # parametrizacao
     
    # Agora que queremos otimizar o SPL, vamos precisar OBRIGATÓRIAMENTE de nodes_target,
    # que vai funcionar como nodes_probe aqui
    isempty(nodes_target) && error("Otim:: nodes_target deve ter ao menos um nó informado")

    # Vamos colocar nodes_target em ordem crescente
    sort!(nodes_target)

    # Precisamos de um material
    isempty(materials) && error("Analise:: at least one material is necessary")

    # Não precisamos do centróide e de vizinhança se for para verificar a derivada
    if !verifica_derivada
    
         # TODO 
         # Calcular centróides e vizinhos somente de elementos de projeto
         # 

         # Calcula a matriz com os centróides de cada elemento da malha
         println("Determinando os centróides dos elementos")
         @time centroides = Centroides(ne,connect,coord,elements_design)

         # TODO 
         # Ver cálculo automático de raio se raio_filtro for nulo
         #
         
         # Obtém os vizinhos de cada elemento da malha
         println("Determinando a vizinhança para um raio de $(raio_filtro)")
         @time vizinhos, pesos = Vizinhanca(ne,centroides,raio_filtro,elements_design)

    end # verifica_derivada

    # Vamos inicializar o vetor de variáveis de projeto.
    # γ = 0 --> ar
    # γ = 1 --> sólido
    #
    # Não podemos começar com todas as posições nulas, pois 
    # isso vai fazer com que a atualização de volume seja 
    # 0*(1+er) = sempre zero.
    # Então, podemos começar com um padrão que seja fisicamente
    # adequado para o problema em questão.
    println("Inicializando o vetor de variáveis de projeto")
    γ = partida*ones(ne)
    println("Ponto de partida = ", partida )
    println()

    # Fixa os valores prescritos de densidade relativa
    Fix_γ!(γ,elements_fixed,values_fixed)
    
    # Grava o arquivo com a distribuição inicial de densidades
    writedlm(arquivo_γ_ini,γ)

    # Vamos avisar que a análise de sensibilidade ainda não está consideranto
    # pressões impostas diretamente
    #if !isempty(pressures) 
    #  error("Aplicação de Pressure é válida para análise, mas não estamos considerando na otimização")
    #end

    # Concatena nodes_open e nodes_pressure
    nodes_mask = sort(vcat(nodes_open,nodes_pressure))
    
    # Posições que não precisam ser calculadas no sistema de equações 
    livres = setdiff(collect(1:nn),nodes_mask)
   
    # Inicializa um arquivo de pós-processamento do gmsh
    etype = connect[:,1]
    Lgmsh_export_init(arquivo_pos,nn,ne,coord,etype,connect[:,3:end])
    
    # Grava a topologia inicial 
    Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter 0")
   
    # Sweep na topologia inicial, para comparação com a otimizada
    MP,_ =  Sweep(nn,ne,coord,connect,γ,vetor_fρ[ponteiro_parametrizacao],vetor_fκ[ponteiro_parametrizacao],freqs,livres,velocities,pressures)

   
    # Exporta por frequência
    for i=1:nf

      # frequência
      f = freqs[i]

      # Exporta
      Lgmsh_export_nodal_scalar(arquivo_pos,abs.(MP[:,i]),"Pressure in $f Hz [abs] - initial topology")

    end

    # Verifica derivadas
    if verifica_derivada

      # Derivada utilizando o procedimento analítico
      MP,K,M =  Sweep(nn,ne,coord,connect,γ,vetor_fρ[ponteiro_parametrizacao],vetor_fκ[ponteiro_parametrizacao],freqs,livres,velocities,pressures)

      # Calcula a derivada da função objetivo em relação ao vetor γ
      dΦ = Derivada(ne,nn,γ,connect,coord,K,M,livres,freqs,pressures,vetor_dfρ[ponteiro_parametrizacao],vetor_dfκ[ponteiro_parametrizacao],nodes_target,MP,elements_design,vA) 

      println("Verificando as derivadas utilizando diferenças finitas centrais...")

      # Derivada numérica
      dnum = Verifica_derivada(γ,nn,ne,coord,connect,vetor_fρ[ponteiro_parametrizacao],vetor_fκ[ponteiro_parametrizacao],freqs,livres,velocities,pressures,nodes_target,elements_design,vA)

      # Relativo
      rel = (dΦ.-dnum)./dnum
      
      # Exporta para a visualização no Gmsh
      Lgmsh_export_element_scalar(arquivo_pos,dΦ,"Analitica")
      Lgmsh_export_element_scalar(arquivo_pos,dnum,"Numerica")
      Lgmsh_export_element_scalar(arquivo_pos,rel,"relativa")

      # Retorna as derivadas
      return dΦ, dnum 

    end

    ########################################################################################################
    ################################ Começo do loop principal de otimização topológica #####################
    ########################################################################################################

    # Volume of each element
    V = Volumes(ne,connect,coord)

    # Total volume sem a parametrização 
    # mas só dos elementos de projeto
    volume_full_projeto = sum(V[elements_design])

    # Target volume
    Vast = vf*volume_full_projeto

    # Sensitivity index in the current iteration
    SN = zeros(ne)    

    # Sensitivity index in the last iterations
    # Aqui vamos  simplificar um pouco a lógica, sacrificando memória
    ESED_F_ANT = zeros(ne,niter)

    # Valor médio da sensibilidade (entre duas iterações)
    ESED_F_media = zeros(ne)

    # Monitora o histórico de volume e de SPL ao longo das 
    # iterações 
    historico_V   = zeros(niter)
    historico_SLP = zeros(niter)

    #############################  Main loop ###########################
    for iter = 1:niter

        # Volume atual da estrutura
        volume_atual = 0.0
        for ele in elements_design
            if γ[ele]≈γ_max
               volume_atual += V[ele]
            end
        end

        # Armazena o volume no histório de volumes
        historico_V[iter] = volume_atual

        # Faz o sweep. A matriz MP tem dimensão nn × nf, ou seja, 
        # cada coluna é o vetor P para uma frequência de excitação
        MP,K,M =  Sweep(nn,ne,coord,connect,γ,vetor_fρ[ponteiro_parametrizacao],vetor_fκ[ponteiro_parametrizacao],freqs,livres,velocities,pressures)

        # Calcula o SPL para esta iteração 
        objetivo = Objetivo(MP,nodes_target,vA)

        println("Iteração       ", iter)
        println("Objetivo       ", objetivo)
        println("Volume atual   ", volume_atual)
        println("Volume target  ", Vast)
        println()

        # Armazena no historico de SPL
        historico_SLP[iter] = objetivo

        # Calcula a derivada da função objetivo em relação ao vetor γ
        dΦ = Derivada(ne,nn,γ,connect,coord,K,M,livres,freqs,pressures,vetor_dfρ[ponteiro_parametrizacao],vetor_dfκ[ponteiro_parametrizacao],nodes_target,MP,elements_design,vA) 
  
        # ESED - Normaliza a derivada do objetivo
        SN[elements_design] .= dΦ[elements_design] ./ V[elements_design]
        
        # Filtro de vizinhança espacial
        ESED_F =  Filtro(vizinhos,pesos,SN,elements_design)

        # Guarda na coluna de ESED_F_media
        ESED_F_ANT[elements_design,iter] .= ESED_F[elements_design]

        # Mean value using the last iterations
        # Aqui temos que ter um cuidado muito importante. O número de colunas 
        # para calcularmos a média deve ser o menor entre iter e nhisto
        pini = max(1,iter-nhisto)
        pfin = max(iter,iter-nhisto) 
        
        # Valor médio 
        ESED_F_media[elements_design] .= mean(ESED_F_ANT[elements_design,pini:pfin],dims=2)
     
        # Visualiza as ESEDS...
        Lgmsh_export_element_scalar(arquivo_pos,SN,"SN")  
        Lgmsh_export_element_scalar(arquivo_pos,ESED_F,"ESED_F")  
        Lgmsh_export_element_scalar(arquivo_pos,ESED_F_media,"ESED_media")  

         #
         # Esquema de relaxação da restrição de volume Eq. 4
         #
         #
         # Um vetor b = [\Delta g ]
         #
         
         # Calcular as variações ΔV^k para essa iteração
         # lembrando que só temos uma restrição (de volume)

         # Lado direito da restrição de volume linearizada em Δγ
         ΔV = Vast - volume_atual

         # Lógica para relaxar a restrição de volume 
         #
         # Parâmetros para comparação 
         α = (1-ϵ1)*volume_atual
         β = (1+ϵ2)*volume_atual

        if Vast < α
           ΔV = -ϵ1*volume_atual
        elseif Vast > β
           ΔV = ϵ2*volume_atual
        end

        # Como só temos uma restrição, 
        b = [ΔV]

         # 
         # Chama a rotina de LP passando  
         # os vetores 
         # c = ESED_F_media
         # A = [ Gradiente da restrição de volume '  ] 
         A = Matrix(transpose(V[elements_design]))

         # Vetor de coeficientes da função objetivo
         c = ESED_F_media[elements_design]

         # Vetor de variáveis de projeto no ponto de ótimo
         # LP(n, c, A, b)
         Δγ =  LP(c, A , b, γ[elements_design])

         # Incrementa os γ
         γ[elements_design] .+= Δγ

         #
         # Cuidado que agora as variáveis de projeto vão estar sempre em 0/1
         #
         # Só para garantir...
         #
         for ele in elements_design
             γ[ele] = round(γ[ele])
         end

         # Grava a topologia para visualização 
         Lgmsh_export_element_scalar(arquivo_pos,γ,"Iter $iter")

    end # iterações externas

    println("Final da otimização, executando a análise SWEEP na topologia otimizada")

   # Roda o sweep na topologia otimizada e exporta para visualização 
   MP,_ =  Sweep(nn,ne,coord,connect,γ,vetor_fρ[ponteiro_parametrizacao],vetor_fκ[ponteiro_parametrizacao],freqs,livres,velocities,pressures)

   # Exporta por frequência
   for i=1:nf

      # frequência
      f = freqs[i]

      # Exporta
      Lgmsh_export_nodal_scalar(arquivo_pos,abs.(MP[:,i]),"Pressure in $f Hz [abs]")

   end

   # Grava um arquivo com os γ finais
   writedlm(arquivo_γ_fin,γ)

   # Retorna o histórico de volume e também o da função objetivo 
   println("Retornando históricos de V e de SLP")
   return historico_V, historico_SLP

end # main_otim

