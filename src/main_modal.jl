

#
# Versão lendo direto do .msh (gmsh)
#
function Modal(meshfile::String, nev=4)

    # Le dados da malha
    nn, coord, ne, connect, materials, nodes_open, velocities, damping = Parsemsh_Daniele(meshfile)

    # Calcula as matrizes globais
    K,M = Monta_KM(nn,ne,coord,connect,materials)

    # DOFs livres do problema
    livres = setdiff(collect(1:nn),nodes_open)

    # Cria as vistas para as posições das matrizes
    KL = K[livres,livres]
    ML = M[livres,livres]

    # Vamos ter que apelar para um solver modal mais porrada
    flag, lamb, X = Solve_Eigen_(KL, ML, nev)

    # Numero efetivo de modos calculados
    nemodos = length(lamb)

    # Avisa se temos a situação em que nem todos os modos
    # retornados são utilizados (por exemplo, uma freq neg)
    # numero efetivo de modos
    if nemodos<nev
        println("Número efetivo de modos ($nemodos) é menor do que o solicitado")
    end

    # Frequências em Hz
    freq =  sqrt.(lamb)./(2*pi)
  
    # Inicializa um arquivo de pós-processamento do gmsh
    nome = "modal.pos"
    etype = connect[:,1]
    Lgmsh_export_init(nome,nn,ne,coord,etype,connect[:,3:end])

    # Exporta os modos 1 até nev
    for  i=1:nemodos

        # Inicializa um vetor com todos os gls
        vv = zeros(nn)

        # Copia o modo para as posições livres
        vv[livres] = X[:,i]

        # Adiciona ao arquivo
        Lgmsh_export_nodal_scalar(nome,vv,"Modo $i, freq $(freq[i])")

    end


end

