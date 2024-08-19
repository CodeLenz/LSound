#
# Versão lendo direto do .msh (gmsh)
#
function Transiente(meshfile::String,metodo=:Bathe)

    # Evita chamar um .geo
    if occursin(".geo",meshfile)
        error("Chamar com .msh..")
    end

    # Verifica se o método é válido
    metodo in [:Bathe, :Newmark] || error("Métodos disponíveis são :Bathe e :Newmark")

    # Le dados da malha
    nn, coord, ne, connect, materials, nodes_open, velocities = Parsemsh_Daniele(meshfile)

    # Aplica condição inicial de pressão
    # funcao(x,y) = Gauss(x,y,0.5,0.5,0.2)
    # funcao(x,y) = Zero(x,y)
    funcao(x,y) = Degrau(x,y,0.5,0.0005,50*0.0011,0.0011)
    U0 = Applica_U0(nn,coord,funcao)

    # Calcula as matrizes globais
    @time K,M = Monta_KM(nn,ne,coord,connect,materials)

    # DOFs livres do problema
    livres = setdiff(collect(1:nn),nodes_open)
  
    # Inicializa um arquivo de pós-processamento do gmsh
    nome = "transiente.pos"
    etype = connect[:,1]
    Lgmsh_export_init(nome,nn,ne,coord,etype,connect[:,3:end])

    # Faz a jogadinha para chamar o Newmark
    P = zeros(nn)
    F(t) = Vetor_P!(t,nn,materials,velocities,coord,connect,P)

    # Chama o Newmark
    C = zeros(nn,nn)
    Δt = 3E-6
    Tf = 0.1

    if metodo===:Newmark
      println("Usando o método Newmark")
      tempos, MP = Newmark(M, C, K, F, livres, Δt, Tf,U0=U0)
    elseif metodo===:Bathe
      println("Usando o método Bathe b1 b2")
      tempos, MP = B1B2Bathe2d(M, C, K, F, livres, Δt, Tf,U0=U0)
    end

    # Inicializa um vetor com todos os gls
    vv = zeros(nn)

    # Exporta os tempos
    for  i=1:size(MP,2)

        # Copia o deslocamento para as posições livres
        vv .= MP[:,i]

        # Adiciona ao arquivo
        Lgmsh_export_nodal_scalar(nome,vv,"Pressao")

    end

    return MP, livres, nodes_open
  
end

