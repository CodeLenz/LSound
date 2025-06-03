"""
 Otim(meshfile::String,freqs=[])

 Basic input:

 meshfile -> arquivo de entrada (.msh)

 Inputs -> freqs, a vector with the frequencies to sweep

"""
function Processa_FRF(meshfile::String,freqs::Vector, γ_opt::Vector)
    
    # Evita chamar um .geo
    occursin(".geo",meshfile) && error("Chamar com .msh..")
    
    # Define os nomes dos arquivos de entrada (yaml) e de saída
    # (pos) em função do nome de meshfile
    arquivo_yaml = meshfile[1:end-3]*"yaml"

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
    
    # Le os dados do arquivo yaml
    raio_filtro, niter, nhisto, er, vf, parametrizacao = Le_YAML(arquivo_yaml)

    # Lista com os elementos que são de projeto
    elements_design = setdiff(1:ne,sort!(elements_fixed))

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
    isempty(nodes_target) && error("Processa_FRF:: nodes_target deve ter ao menos um nó informado")

    # Vamos colocar nodes_target em ordem crescente
    sort!(nodes_target)

    # Precisamos de um material
    isempty(materials) && error("Processa_FRF:: at least one material is necessary")
    
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

    # Concatena nodes_open e nodes_pressure
    nodes_mask = sort(vcat(nodes_open,nodes_pressure))
    
    # Posições que não precisam ser calculadas no sistema de equações 
    livres = setdiff(collect(1:nn),nodes_mask)
   
    # Sweep na topologia inicial, para comparação com a otimizada
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities,pressures)

    # Número de frequências
    nf = length(freqs)

    # Cria um vetor para armazenar o SPLn em cada frequência considerada na otimização
    historico_SLPn_inicial = zeros(length(freqs))
   
    # Processa o SPLn por frequência
    for i=1:nf

      # Calcula o valor e armazena
      historico_SLPn_inicial[i] = SPLn(MP[:,i],20E-6)

    end

    # Aloca o vetor de saída - SPLn para cada frequência
    historico_SLPn_final = zeros(nf)

    
    # Roda o sweep na topologia otimizada e exporta para visualização 
    MP,_ =  Sweep(nn,ne,coord,connect,γ_opt,fρ,fκ,freqs,livres,velocities,pressures)
    
    # Calcula o SLPn em cada uma das frequências 
    for i=1:nf

      historico_SLPn_final[i] = SPLn(MP[:,i],20E-6)

   end

   # Retorna os históricos
   return  freqs, historico_SLPn_inicial, historico_SLPn_final

end 
