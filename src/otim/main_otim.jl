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
    nn, coord, ne, connect, materials, nodes_open, velocities, damping, nodes_probe = Parsemsh_Daniele(meshfile)

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
    # C = Matriz_C(nn,damping,coord,connect)

    
    # Verificamos se existem frequências sendo informadas
    if isempty(freqs)
        error("Analise Harmonica:: freqs deve ser um vetor não vazio")
    end

    # Número de frequências
    nω = length(freqs)
        
    # Número de nós a monitorar
    np = length(nodes_probe) 

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
    for f in freqs

        # Converte a freq para rad/s
        ω = 2*pi*f

        # Monta a matriz de rigidez dinâmica
        # Kd = K[livres,livres] .+ im*ω*C[livres,livres] .- (ω^2)*M[livres,livres]
        Kd = K[livres,livres]  .- (ω^2)*M[livres,livres]

        # Monta o vetor de forças, que depende da frequência  
        Vetor_P!(0.0,nn,materials,velocities,coord,connect,P,ω=ω)

        # Soluciona 
        U[livres] .= Kd\P[livres]

        # Adiciona ao arquivo
        Lgmsh_export_nodal_scalar(nome,abs.(U),"Pressão em $f Hz [abs]")
        
        # Armazena os resultados na matriz de monitoramento
        monitor[:,contador] .= U[nodes_probe]

        # Incrementa o contador
        contador += 1

    end

    return monitor

end
