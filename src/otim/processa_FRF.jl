"""
 Processa_FRF(meshfile::String,freqs=[],γ_opt)

 Calcula o SPLn para as frequências em freqs, tanto para a topologia 
 inicial quanto para a otimizada (γ_opt) 

 Basic input:

 meshfile -> arquivo de entrada (.msh)

 freqs ->  vector with the frequencies to sweep 
 
 γ_opt -> optimal distribution of γ

"""
function Processa_FRF(meshfile::String,freqs::Vector)
    
    # Evita chamar um .geo
    occursin(".geo",meshfile) && error("Chamar com .msh..")
    
    # Define os nomes dos arquivos de entrada (yaml) e de saída
    # (pos) em função do nome de meshfile
    arquivo_yaml = meshfile[1:end-3]*"yaml"

    # Arquivos que contém as distribuições inicial/otimizada de γ
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
    
    # Le os dados do arquivo yaml
    raio_filtro, niter, nhisto, ϵ1, ϵ2, vf, parametrizacao, γ_min, γ_max, partida= Le_YAML(arquivo_yaml)

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
    #  println("Utilizando a parametrização de DUHRING")
    #  fρ(γ)  = fρ_duhring(γ)
    #  dfρ(γ) = dfρ_duhring(γ)
    #  fκ(γ)  = fκ_duhring(γ)
    #  dfκ(γ) = dfκ_duhring(γ)
    end
     
    # Agora que queremos otimizar o SPL, vamos precisar OBRIGATÓRIAMENTE de nodes_target,
    # que vai funcionar como nodes_probe aqui
    isempty(nodes_target) && error("Processa_FRF:: nodes_target deve ter ao menos um nó informado")

    # Vamos colocar nodes_target em ordem crescente
    sort!(nodes_target)

    # Precisamos de um material
    isempty(materials) && error("Processa_FRF:: at least one material is necessary")
    
    # Lê as variáveis de projeto iniciais
    γ = vec(readdlm(arquivo_γ_ini))  
   
    # Concatena nodes_open e nodes_pressure
    nodes_mask = sort(vcat(nodes_open,nodes_pressure))
    
    # Posições que não precisam ser calculadas no sistema de equações 
    livres = setdiff(collect(1:nn),nodes_mask)
   
    # Sweep na topologia inicial, para comparação com a otimizada
    MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities,pressures)

    # Número de frequências
    nf = length(freqs)

    # Cria um vetor para armazenar o SPLn em cada frequência considerada na otimização
    FRF_SLPn_inicial = zeros(length(freqs))
   
    # Processa o SPLn por frequência
    for i=1:nf

      # Calcula o valor e armazena
      FRF_SLPn_inicial[i] = SPLn(MP[:,i],20E-6)

    end

    # Aloca o vetor de saída - SPLn para cada frequência
    FRF_SLPn_final = zeros(nf)

    # Le as variáveis de projeto finais (otimizadas)
    γ_opt = vec(readdlm(arquivo_γ_fin))  
    
    # Roda o sweep na topologia otimizada e exporta para visualização 
    MP,_ =  Sweep(nn,ne,coord,connect,γ_opt,fρ,fκ,freqs,livres,velocities,pressures)
    
    # Calcula o SLPn em cada uma das frequências 
    for i=1:nf

      FRF_SLPn_final[i] = SPLn(MP[:,i],20E-6)

    end

   # Retorna o vetor com as frequências, a varredura inicial e a final
   return  freqs, FRF_SLPn_inicial, FRF_SLPn_final

end 