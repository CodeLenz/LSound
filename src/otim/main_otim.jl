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
    γ = 0.1*rand(ne) #zeros(ne) 

    # DOFs livres do problema
    livres = setdiff(collect(1:nn),nodes_open)

    # Verificamos se existem frequências sendo informadas
    if isempty(freqs)
        error("Analise Harmonica:: freqs deve ser um vetor não vazio")
    end

    # Inicializa um arquivo de pós-processamento do gmsh
    nome = "otim.pos"
    etype = connect[:,1]
    Lgmsh_export_init(nome,nn,ne,coord,etype,connect[:,3:end])
    
    # Faz o sweep
    target,K,M =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities) 

    # Gera a visualização 
    # Lgmsh_export_nodal_scalar(nome,abs.(U),"Pressão em $f Hz [abs]") 

    # Calcula a função objetivo SPL_w
    objetivo = Objetivo(target,nodes_target)

    # Calcula a derivada da função objetivo em relação ao vetor γ
    dΦ = Derivada(ne,nn,γ,connect,coord,K,M,livres,freqs,dfρ,dfκ,nodes_target,target) 

    # Vamos validar a derivada usando diferenças finitas
    function f_(γ,nn,ne,coord,connect,fρ,fκ,freqs,livres,velocities)

        target,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities) 

        # Calcula a função objetivo SPL_w
        objetivo = Objetivo(target,nodes_target)

        return objetivo

    end
    f(γ) = f_(γ,nn,ne,coord,connect,fρ,fκ,freqs,livres,velocities)
    d_numerica = df(γ,f)
    

    return target, objetivo, dΦ, d_numerica

end


function Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities)

    # Calcula as matrizes globais
    K,M = Monta_KM2(ne,coord,connect,γ,fρ,fκ)
    
    # E a de amortecimento
    # TODO adicionar amortecimento depois
    # C = Matriz_C(nn,damping,coord,connect)

    # Número de frequências
    nω = length(freqs)

    # Aloca matriz com os valores a serem monitorados
    target = zeros(ComplexF64,nn,nω)

    # Aloca o vetor de forças 
    P = zeros(nn)

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
        target[livres,contador] .= Kd\P[livres]
  
        # Incrementa o contador
        contador += 1

    end

    # Retorna target, K e M
    return target, K, M

end