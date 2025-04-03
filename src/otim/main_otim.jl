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

    # Evita chamar um .geo
    if occursin(".geo",meshfile)
        error("Chamar com .msh..")
    end

    # Le dados da malha
    nn, coord, ne, connect, materials, nodes_open, velocities, damping, nodes_probe, nodes_target = Parsemsh_Daniele(meshfile)

    # Agora que queremos otimizar o SPL, vamos precisar OBRIGATÓRIAMENTE de nodes_target,
    # que vai funcionar como nodes_probe aqui
    isempty(nodes_target) && error("Otim:: nodes_target deve ter ao menos um nó informado")

    # Número de nós em nodes_target
    nt = length(nodes_target)

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
    if isempty(materials)
        error("Analise:: at least one material is necessary")
    end

    # Vamos inicializar o vetor de variáveis de projeto toda em 0.0
    # ou seja, ar
    γ = zeros(ne) 

    # Calcula as matrizes globais
    K,M = Monta_KM2(ne,coord,connect,γ,fρ,fκ)

    # DOFs livres do problema
    livres = setdiff(collect(1:nn),nodes_open)
    
    # Inicializa um arquivo de pós-processamento do gmsh
    nome = "otim.pos"
    etype = connect[:,1]
    Lgmsh_export_init(nome,nn,ne,coord,etype,connect[:,3:end])
    
    # Aloca o vetor de forças (para Harmônico e transiente)
    P = zeros(nn)

    # E a de amortecimento
    # TODO adicionar amortecimento depois
    # C = Matriz_C(nn,damping,coord,connect)
    
    # Verificamos se existem frequências sendo informadas
    if isempty(freqs)
        error("Analise Harmonica:: freqs deve ser um vetor não vazio")
    end

    # Número de frequências
    nω = length(freqs)
        
    # Número de nós a monitorar
    np = length(nodes_probe) 

    # Aloca matriz com os valores a serem monitorados
    target = zeros(ComplexF64,nn,nω)

    # pre-aloca vetor de resposta
    U = zeros(ComplexF64, nn)

    # Loop pelas frequências
    contador = 1
    for f in freqs

        # Converte a freq para rad/s
        ω = 2*pi*f

        # Monta a matriz de rigidez dinâmica
        # Kd = K[livres,livres] .+ im*ω*C[livres,livres] .- (ω^2)*M[livres,livres]
        Kd = K[livres,livres]  .- (ω^2)*M[livres,livres]

        # Monta o vetor de forças, que depende da frequência  
        Vetor_P!(0.0,velocities,coord,connect,P,ω=ω)

        # Soluciona apenas para os gls livres do problema
        U[livres] .= Kd\P[livres]

        # Adiciona ao arquivo
        Lgmsh_export_nodal_scalar(nome,abs.(U),"Pressão em $f Hz [abs]")
        
        # Armazena os resultados na matriz de monitoramento
        target[:,contador] .= U

        # Incrementa o contador
        contador += 1

    end

    # Calcula a função objetivo SPL_w
    objetivo = Objetivo(target,nodes_target)

    # Calcula a derivada da função objetivo em relação ao vetor γ
    dΦ = Derivada(ne,nn,γ,connect,K,M,livres,freqs,dfρ,dfκ,nodes_target,target) 

    return target, objetivo, dΦ

end
