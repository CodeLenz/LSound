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

 scale -> [1.0 ; 1.0 ; 1.0] scale to apply to the geometry (1.0 meter)


 Inputs -> freqs, a vector with the frequencies to sweep


"""
function Otim(meshfile::String,freqs::Vector;verifica_derivada=false)
    
    # Evita chamar um .geo
    occursin(".geo",meshfile) && error("Chamar com .msh..")
    
    # Define os nomes dos arquivos de entrada (yaml) e de saída
    # (pos) em função do nome de meshfile
    arquivo_yaml = meshfile[1:end-3]*"yaml"
    arquivo_pos  = meshfile[1:end-3]*"pos"

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
    raio_filtro, niter, nhisto, er, vf, parametrizacao = Le_YAML(arquivo_yaml)

    # Seleciona as rotinas de parametrização de material de acordo com 
    # a opção 
    if parametrizacao=="PEREIRA"
         println("Utilizando a parametrização de PEREIRA")
         fρ(γ)  = fρ_pereira(γ) #,ψ, ρ_ar = ρ_ar, ρ2 = ρ_solido)
         dfρ(γ) = dfρ_pereira(γ)
         fκ(γ)  = fκ_pereira(γ)
         dfκ(γ) = dfκ_pereira(γ)
    #elseif parametrizacao=="DUHRING"
    #     println("Utilizando a parametrização de DUHRING")
    #     fρ(γ)  = fρ_duhring(γ)
    #     dfρ(γ) = dfρ_duhring(γ)
    #     fκ(γ)  = fκ_duhring(γ)
    #     dfκ(γ) = dfκ_duhring(γ)
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
    
         # Calcula a matriz com os centróides de cada elemento da malha
         println("Determinando os centróides dos elementos")
         @time centroides = Centroides(ne,connect,coord)

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
    println("Utilizando a fração de volume como ponto de partida")
    γ = vf*ones(ne) #+ 1E-2*randn(ne)
    
    # Fixa os valores prescritos de densidade relativa
    Fix_γ!(γ,elements_fixed,values_fixed)
    
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
    @time MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities,pressures)

    # Número de frequências
    nf = length(freqs)
   
    # Exporta por frequência
    for i=1:nf

      # frequência
      f = freqs[i]

      # Exporta
      Lgmsh_export_nodal_scalar(arquivo_pos,abs.(MP[:,i]),"Pressão em $f Hz [abs] - topolgia inicial")

    end

    # Verifica derivadas
    if verifica_derivada

      # Derivada utilizando o procedimento analítico
      MP,K,M =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities,pressures)

      # Calcula a derivada da função objetivo em relação ao vetor γ
      dΦ = Derivada(ne,nn,γ,connect,coord,K,M,livres,freqs,pressures,dfρ,dfκ,nodes_target,MP,elements_design) 

      println("Verificando as derivadas utilizando diferenças finitas centrais...")

      # Derivada numérica
      dnum = Verifica_derivada(γ,nn,ne,coord,connect,fρ,fκ,freqs,livres,velocities,pressures,nodes_target,elements_design)

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

    # Sensitivity index in the last iteration
    ESED_F_ANT = zeros(ne)

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
        volume_atual = sum(γ[elements_design].*V[elements_design])

        # Armazena o volume no histório de volumes
        historico_V[iter] = volume_atual

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

        # Faz o sweep. A matriz MP tem dimensão nn × nf, ou seja, 
        # cada coluna é o vetor P para uma frequência de excitação
        MP,K,M =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities,pressures)

        # Calcula o SPL para esta iteração 
        objetivo = Objetivo(MP,nodes_target)

        println("Iteração       ", iter)
        println("Objetivo       ",objetivo)
        println("Volume atual   ", volume_atual)
        println("Volume próxima ", vol)
        println("Volume target  ", Vast)
        println()

        # Armazena no historico de SPL
        historico_SLP[iter] = objetivo

        # Calcula a derivada da função objetivo em relação ao vetor γ
        dΦ = Derivada(ne,nn,γ,connect,coord,K,M,livres,freqs,pressures,dfρ,dfκ,nodes_target,MP,elements_design) 
  
        # Zera a derivada dos elementos fixos
        # Fix_D!(dΦ,elements_fixed)

        # Como podemos ter variação de sinal na derivada, devemos tomar cuidado 
        # com a lógica dos esquemas que funcionam para compliance (derivadas sempre
        # negativas)

        # ESED - Normaliza a derivada do objetivo
        #        e corrige o sinal para a definição do Índice de Sensibilidade
        SN = -dΦ ./ V
        
        # Filtro de vizinhança espacial
        ESED_F =  Filtro(ne,vizinhos,pesos,SN,elements_design)

        # Zera os valores fixos
        # Fix_D!(ESED_F,elements_fixed)

        # Mean value using the last iteration
        if iter > 1
           ESED_F_media .= (ESED_F .+ ESED_F_ANT)./2
        else
           ESED_F_media .= ESED_F
        end

        # Store the value for the next iteration
        ESED_F_ANT .= ESED_F

        # Update the relative densities
        γn, niter_beso = BESO(γ, ESED_F_media, V, vol, elements_design)

        # Garante que os elementos fixos não tenham sido alterados
        # Fix_γ!(γn,elements_fixed,values_fixed)

        # Se niter_beso for nula, então o problema stagnou
        if niter_beso==0
           println("BESO não atualizou as variáveis")
           break
        end

        # Grava a topologia para visualização 
        Lgmsh_export_element_scalar(arquivo_pos,γn,"Iter $iter")

        # Atualiza o γ para a próxima iteração 
        γ .= γn
     
    end # iterações externas

    println("Final da otimização, executando a análise SWEEP na topologia otimizada")

    # Roda o sweep na topologia otimizada e exporta para visualização 
    @time MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities,pressures)

    # Número de frequências
    nf = length(freqs)
    
    # Exporta por frequência
    for i=1:nf

        # frequência
        f = freqs[i]

        # Exporta
        Lgmsh_export_nodal_scalar(arquivo_pos,abs.(MP[:,i]),"Pressão em $f Hz [abs]")

    end

    # Retorna o histórico de volume e também o da função objetivo 
    return historico_V, historico_SLP

end # main_otim
