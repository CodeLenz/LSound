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

 verifica_derivada -> bool 

 Inputs -> freqs, a vector with the frequencies to sweep

"""
function Otim(meshfile::String,freqs::Vector;verifica_derivada=false)
    
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

    # Arquivo .yaml
    isfile(arquivo_yaml) || error("Otim:: arquivo de entrada $(arquivo_yaml) não existe")

    # Le dados da malha
    nn, coord, ne, connect, materials, nodes_open, velocities, nodes_pressure, pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed = Parsemsh_Daniele(meshfile)

    # Lista com os elementos que são de projeto
    elements_design = setdiff(1:ne,sort!(elements_fixed))

    # Le os dados do arquivo yaml
    raio_filtro, niter, nhisto, er, vf, parametrizacao, γ_min, γ_max = Le_YAML(arquivo_yaml)

    # Seleciona as rotinas de parametrização de material de acordo com 
    # a opção 
    @show parametrizacao

    # Vetores de funções 
    vetor_fρ  = [fρ_pereira, fρ_duhring]
    vetor_dfρ = [dfρ_pereira, dfρ_duhring]
    vetor_fκ  = [fκ_pereira, fκ_duhring]
    vetor_dfκ = [dfκ_pereira, dfκ_duhring]


    # 1 para PEREIRA e 2 para Duhring
    ponteiro_parametrizacao = 1
 

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
    end
     
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


    end

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
    #println("Utilizando a fração de volume como ponto de partida")
    γ = γ_min*ones(ne) #+ 1E-2*randn(ne)

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

    # Número de frequências
    nf = length(freqs)
   
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
      dΦ = Derivada(ne,nn,γ,connect,coord,K,M,livres,freqs,pressures,vetor_dfρ[ponteiro_parametrizacao],vetor_dfκ[ponteiro_parametrizacao],nodes_target,MP,elements_design) 

      println("Verificando as derivadas utilizando diferenças finitas centrais...")

      # Derivada numérica
      dnum = Verifica_derivada(γ,nn,ne,coord,connect,vetor_fρ[ponteiro_parametrizacao],vetor_fκ[ponteiro_parametrizacao],freqs,livres,velocities,pressures,nodes_target,elements_design)

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

    # Define o γn para as atualizações das variáveis de projeto
    γn = zeros(ne)

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

    # Define vol fora do loop para podermos recuperar depois
    vol = 0.0

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

        
        println("Iteração       ", iter)
        #println("Objetivo       ", objetivo)
        println("Volume atual   ", volume_atual)
        #println("Volume próxima ", vol)
        println("Volume target  ", Vast)
        println()

        # Calcula a derivada da função objetivo em relação ao vetor γ
        dΦ = Derivada(ne,nn,γ,connect,coord,K,M,livres,freqs,pressures,vetor_dfρ[ponteiro_parametrizacao],vetor_dfκ[ponteiro_parametrizacao],nodes_target,MP,elements_design) 
  
        # Zera a derivada dos elementos fixos
        Fix_D!(dΦ,elements_fixed)

        # Como podemos ter variação de sinal na derivada, devemos tomar cuidado 
        # com a lógica dos esquemas que funcionam para compliance (derivadas sempre
        # negativas)

        # ESED - Normaliza a derivada do objetivo
        #        e corrige o sinal para a definição do Índice de Sensibilidade
        SN[elements_design] .= -dΦ[elements_design] ./ V[elements_design]
        
        # Filtro de vizinhança espacial
        ESED_F =  Filtro(vizinhos,pesos,SN,elements_design)

        # Zera os valores fixos
        # Fix_D!(ESED_F,elements_fixed)

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

        # Update the relative densities
        # Método baseado na biseção
        # γn, niter_beso = BESO(γ, ESED_F_media, V, vol, elements_design,γ_min=γ_min,γ_max=γ_max)

        # Loop, variando o er para que o objetivo da próxima iteração seja menor
        # do que o objetivo atual
        #
        # Mudando o er, mudamos o vol e, com isso, mudamos o número de elementos que 
        # podem ser colocados ou retirados. 
        println("Entrando no LS com objetivo ", objetivo)

        objetivo_slp = 0.0

        # Fator de redução do passo 
        τ = 0.005 # --> QUAL SERIA O VALOR CORRETO DE UTILIZAR AQUI?
        
        # Constante 'c'
        c = 1E-4

        # Passo inicial
        α_0 = er

        for ls=1:10
 
            # Lógica de atualização do er
            # ARMIJO BACKTRACKING LINE SEARCH

            # Volume a ser utilizado como limite para esta iteração
            # Este é o principal parâmetro de controle para a estabilidade
            # do processo de otimização
            vol = Calcula_volume_er(er,volume_atual,Vast)

            # Faz o sweep. A matriz MP tem dimensão nn × nf, ou seja, 
            # cada coluna é o vetor P para uma frequência de excitação
            MP,K,M =  Sweep(nn,ne,coord,connect,γ,vetor_fρ[ponteiro_parametrizacao],vetor_fκ[ponteiro_parametrizacao],freqs,livres,velocities,pressures)

            # BESO Clássico
            γn .= BESO3(γ, ESED_F_media,V,vol,elements_design,xmin=γ_min,xmax=γ_max)
            
            # Calcula o SPL para esta iteração 
            objetivo_slp = Objetivo(MP,nodes_target)

            println("Objetivo  LS     ", objetivo_slp)

            # comparação do objetivo_ls com o objetivo atual

            # "Gradiente" para a condição de ARMIJO
            # Será SN[elements_design] ou dΦ[elements_design] ?
            # gradiente = -dΦ[elements_design]  # ???
            # Direção de descida 
            # descida = - gradiente     # --> deve ser oposto ao gradiente...somente sinal?

            # Valor de t
            # t = ?  --> t = -c*m   --> c * gradiente' * descida   (considerar o gradiente no ponto)

            # Condição de ARMIJO é satisfeita?
            contador = 0
             
            while true # Pode usar assim ? ou o correto seria a condição de "ARMIJO"?
                # Calcula a nova objetivo após o passo
                novo_objetivo = Objetivo(MP, nodes_target)

                # Verificar a condição de ARMIJO
                if objetivo_slp - novo_objetivo >= α_0 * c * gradiente' * descida
                    println("Condição de ARMIJO satisfeita. Passo aceito = ", α_0)
               
                    # Condição for satisfeita ?? sair do loop
                break 
                else
                    # Caso contrário, reduzir o passo multiplicando por τ e tentar novamente
                    α_0 *= τ
                    contador += 1
                    println("Condição de ARMIJO não satisfeita. Reduzindo ", α_0, " na iteração ", contador)
                end # if

                # Limitar o número de tentativas ??
                if contador > 50 
                    println("Limite de tentativas atingido. Parando a busca!")
                    break
                end

            end # while

        end #ls

        # Atualiza o objetivo
        objetivo = objetivo_slp 

        # Armazena no historico de SPL
        historico_SLP[iter] = objetivo


        # Grava a topologia para visualização 
        Lgmsh_export_element_scalar(arquivo_pos,γn,"Iter $iter")

        # Se niter_beso for nula, então o problema stagnou
        if norm(γn-γ,Inf) ≈ 0 
           γ .= γn
           println("BESO não atualizou as variáveis")
           break
        end

        # Atualiza o γ para a próxima iteração 
        γ .= γn

        
    end # iterações externas

    println("Final da otimização, executando a análise SWEEP na topologia otimizada")

   # Roda o sweep na topologia otimizada e exporta para visualização 
   MP,_ =  Sweep(nn,ne,coord,connect,γ,vetor_fρ[ponteiro_parametrizacao],vetor_fκ[ponteiro_parametrizacao],freqs,livres,velocities,pressures)

   # Número de frequências
   nf = length(freqs)
    
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
   return historico_V, historico_SLP

end # main_otim

#
# Rotina que calcula o volume da próxima iteração, baseado no conceito de 
# er (evolutionary rate, ou taxa de evolução)
#
function Calcula_volume_er(er,volume_atual,Vast)

    # Volume a ser utilizado como limite para esta iteração
    # Este é o principal parâmetro de controle para a estabilidade
    # do processo de otimização
    if volume_atual > Vast
        vol = max(Vast, volume_atual*(1.0 - er))
    elseif volume_atual < Vast
        vol = min(Vast,volume_atual*(1 + er))
    else
        vol = Vast
    end 
 
    # Caso tenhamos um volume nulo, decorrente de uma 
    # inicialização com todos os γ=γ_min, devemos utilizar
    # vol como sendo algum valor pré-determinado, pois o 
    # if das linhas anteriores vai retornar zero
    if vol==0
       vol = er*Vast
    end

    # Retorna o volume a ser atendido pelo BESO
    return vol

end