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
    occursin(".geo",meshfile) || error("Chamar com .msh..")
    
    # Verificamos se existem frequências sendo informadas
    isempty(freqs) && error("Analise Harmonica:: freqs deve ser um vetor não vazio")

    # Evita passar as frequências como algo diferente de um Vetor de floats
    isa(freqs,Vector{Float64}) || error("freqs deve ser um vetor de floats")
    
    # Le dados da malha
    nn, coord, ne, connect, materials, nodes_open, velocities, damping, nodes_probe, nodes_target = Parsemsh_Daniele(meshfile)

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
    
    # Vamos inicializar o vetor de variáveis de projeto toda em 0.0
    # ou seja, ar
    γ = 0.1*rand(ne) #zeros(ne) 

    # DOFs livres do problema
    livres = setdiff(collect(1:nn),nodes_open)

    # Inicializa um arquivo de pós-processamento do gmsh
    nome = "otim.pos"
    etype = connect[:,1]
    Lgmsh_export_init(nome,nn,ne,coord,etype,connect[:,3:end])
    
    # Faz o sweep. A matriz MP tem dimensão nn × nf, ou seja, 
    # cada coluna é o vetor P para uma frequência de excitação
    MP,K,M =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities) 

    # Calcula a função objetivo SPL_w
    objetivo = Objetivo(MP,nodes_target)

    # Calcula a derivada da função objetivo em relação ao vetor γ
    dΦ = Derivada(ne,nn,γ,connect,coord,K,M,livres,freqs,dfρ,dfκ,nodes_target,MP) 

    # Vamos validar a derivada usando diferenças finitas
    function f_(γ,nn,ne,coord,connect,fρ,fκ,freqs,livres,velocities)

        MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities) 

        # Calcula a função objetivo SPL_w
        objetivo = Objetivo(MP,nodes_target)

        return objetivo

    end
    f(γ) = f_(γ,nn,ne,coord,connect,fρ,fκ,freqs,livres,velocities)

    # Calcula a derivada por DFC
    d_numerica = df(γ,f,1E-8)
    

    return target, objetivo, dΦ, d_numerica

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