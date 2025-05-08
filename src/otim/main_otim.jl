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

 scale -> [1.0 ; 1.0 ; 1.0] scale to apply to the geometry (1.0 meter)


 Inputs -> freqs, a vector with the frequencies to sweep


"""
function Otim(meshfile::String,freqs::Vector,scale=[1.0;1.0;1.0])

    # Um dos dados de entrada para a otimização é a fração de volume
    # como γ = 1 é material sólido, vf está dizendo quanto do projeto
    # final terá de material sódido (1-vf será a quantidade de ar)
    vf = 0.5

    # Evita chamar um .geo
    occursin(".geo",meshfile) && error("Chamar com .msh..")
    
    # Verificamos se existem frequências sendo informadas
    isempty(freqs) && error("Analise Harmonica:: freqs deve ser um vetor não vazio")

    # Evita passar as frequências como algo diferente de um Vetor de floats
    isa(freqs,Vector{Float64}) || error("freqs deve ser um vetor de floats")
    
    # Verifica se os arquivos de entrada existem
    isfile(meshfile) || throw("Otim:: arquivo de entrada $meshfile não existe")

    # Arquivo .yaml
    arquivo_yaml = meshfile[1:end-3]*"yaml"
    isfile(arquivo_yaml) || throw("Otim:: arquivo de entrada $(arquivo_yaml) não existe")

    # Le dados da malha
    nn, coord, ne, connect, materials, nodes_open, velocities, damping, nodes_probe, nodes_target = Parsemsh_Daniele(meshfile)

    # Le os dados do arquivo yaml
    raio_filtro, niter, er = Le_YAML(arquivo_yaml)

    # Agora que queremos otimizar o SPL, vamos precisar OBRIGATÓRIAMENTE de nodes_target,
    # que vai funcionar como nodes_probe aqui
    isempty(nodes_target) && error("Otim:: nodes_target deve ter ao menos um nó informado")

    # Vamos evitar coordenadas negativas 
    for i=1:3  
        minx = minimum(coord[:,i])
        if minx<0
           coord[:,i] .= coord[:,i] .- minx 
        end
    end

    # Apply scale to the coordinates
    coord[:,1] ./= scale[1]
    coord[:,2] ./= scale[2]
    coord[:,3] ./= scale[3]

    # Precisamos de um material
    isempty(materials) && error("Analise:: at least one material is necessary")
    
    # Calcula a matriz com os centróides de cada elemento da malha
    @time centroides = Centroides(ne,connect,coord)

    # Obtém os vizinhos de cada elemento da malha
    @time vizinhos, pesos = Vizinhanca(ne,centroides,raio_filtro)

    # Vamos inicializar o vetor de variáveis de projeto.
    # γ = 0 --> ar
    # γ = 1 --> sólido
    #
    # Não podemos começar com todas as posições nulas, pois 
    # isso vai fazer com que a atualização de volume seja 
    # 0*(1+er) = sempre zero.
    # Então, podemos começar com um padrão que seja fisicamente
    # adequado para o problema em questão.
    γ = rand(ne)

    # DOFs livres do problema
    livres = setdiff(collect(1:nn),nodes_open)

    # Inicializa um arquivo de pós-processamento do gmsh
    nome = "otim.pos"
    etype = connect[:,1]
    Lgmsh_export_init(nome,nn,ne,coord,etype,connect[:,3:end])
    
    # Grava a topologia inicial 
    Lgmsh_export_element_scalar(nome,γ,"Iter 0")

    # Começo do loop principal de otimização topológica

    # Volume of each element
    V = Volumes(ne,connect,coord)

    # Total volume sem a parametrização 
    volume_full = sum(V)

    # Target volume
    Vast = vf*volume_full

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
        volume_atual = sum(γ.*V)

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

        println("Iteração       ", iter)
        println("Volume atual   ", volume_atual)
        println("Volume próxima ", vol)
        println("Volume target  ", Vast)
        println()

        # Agora podemos calcular a resposta do problema - Equilíbrio

        # Faz o sweep. A matriz MP tem dimensão nn × nf, ou seja, 
        # cada coluna é o vetor P para uma frequência de excitação
        MP,K,M =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities)

        # Calcula o SPL para esta iteração 
        objetivo = Objetivo(MP,nodes_target)

        # Armazena no historico de SPL
        historico_SLP[iter] = objetivo

        # Calcula a derivada da função objetivo em relação ao vetor γ
        dΦ = Derivada(ne,nn,γ,connect,coord,K,M,livres,freqs,dfρ,dfκ,nodes_target,MP) 
  
        # Valores extremos da derivada
        max_dΦ = maximum(dΦ)

        # Se o valor máximo for positivo, transladamos todos os valores para 
        # ficarem negativos
        #if max_dΦ > 0 
        #    dΦ = dΦ .- 1.1*max_dΦ
        #end
    

        # Como podemos ter variação de sinal na derivada, devemos tomar cuidado 
        # com a lógica dos esquemas que funcionam para compliance (derivadas sempre
        # negativas)

        # ESED - Normaliza a derivada do objetivo
        #        e corrige o sinal para a definição do Índice de Sensibilidade
        SN = -dΦ ./ V
        
        # Filtro de vizinhança espacial
        ESED_F =  Filtro(ne,vizinhos,pesos,SN)

        # Mean value using the last iteration
        if iter > 1
           ESED_F_media .= (ESED_F .+ ESED_F_ANT)./2
        else
           ESED_F_media .= ESED_F
        end

        # Store the value for the next iteration
        ESED_F_ANT .= ESED_F_media 

        # Update the relative densities
        γn, niter_beso = BESO(γ, ESED_F_media, V, vol)

        # Se niter_beso for nula, então o problema stagnou
        if niter_beso==0
           println("BESO não atualizou as variáveis")
        end

        # Grava a topologia para visualização 
        Lgmsh_export_element_scalar(nome,γn,"Iter $iter")

        # Compara a variação de variáveis de projeto
        # entre as duas iterações
        println("Variação máxima de γ ", norm(γn-γ,Inf))

        # Atualiza o γ para a próxima iteração 
        γ .= γn
     
    end # iterações externas

    return historico_V, historico_SLP

    # Calcula a função objetivo SPL_w
    # objetivo = Objetivo(MP,nodes_target)

end # main_otim

#
# Programar depois para fazer a validação das derivadas por DFC
#
function Verifica_derivada(γ,nn,ne,coord,connect,fρ,fκ,freqs,livres,velocities)
    
    # Vamos validar a derivada usando diferenças finitas
    function f_(γ,nn,ne,coord,connect,fρ,fκ,freqs,livres,velocities)

        MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities) 

        # Calcula a função objetivo SPL_w
        objetivo = Objetivo(MP,nodes_target)

        return objetivo

    end
    f(γ) = f_(γ,nn,ne,coord,connect,fρ,fκ,freqs,livres,velocities)

    # Calcula a derivada por DFC
    println("Entrando em numérica")
    d_numerica = df(γ,f,1E-8)

end



#
# Realiza uma sequência de análises harmônicas em uma lista de nf 
# frequências de excitação e guarda a solução em uma matriz nn × nf 
#
function Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities)

    # Calcula as matrizes globais
    K,M = Monta_KM2(ne,coord,connect,γ,fρ,fκ)
    
    # E a de amortecimento
    # TODO adicionar amortecimento depois
    # C = Matriz_C(nn,damping,coord,connect)

    # Número de frequências
    nf = length(freqs)

    # Aloca matriz com os valores a serem monitorados
    MP = zeros(ComplexF64,nn,nf)

    # Aloca o vetor de forças 
    P = Array{Float64}(undef,nn)

    # Loop pelas frequências
    contador = 1
    for f in freqs

        # Converte a freq de Hz para rad/s
        ω = 2*pi*f

        # Monta a matriz de rigidez dinâmica
        # Kd = K[livres,livres] .+ im*ω*C[livres,livres] .- (ω^2)*M[livres,livres]
        Kd = K[livres,livres]  .- (ω^2)*M[livres,livres]

        # Monta o vetor de forças, que depende da frequência  
        Vetor_P!(0.0,velocities,coord,connect,P,ω=ω)

        # Soluciona apenas para os gls livres do problema
        MP[livres,contador] .= Kd\P[livres]
  
        # Incrementa o contador
        contador += 1

    end

    # Retorna MP, K e M
    return MP, K, M

end