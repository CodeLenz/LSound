 #
  # Devolve os vetores nos, X e Y com os nós 
  # e as coordenadas do elemento
  #
  function Nos_Coordenadas(ele,etype,coord,connect)
  
    # Mapeia o número de nós pelo etype
    # 2 - tri - 3nós
    # 3 - quad - 4nós
    nnos = Lgmsh_nodemap()[etype]

    # Descobre os nós do elemento
    nos = connect[ele,3:2+nnos]

    # Aloca vetores de coordenadas
    X = Vector{Int64}(undef,nnos)
    Y = Vector{Int64}(undef,nnos)

    # Descobre as coordenadas x e y de cada nó
    # do elemento
    for i in LinearIndices(nos)

        # Nó
        no = nos[i]

        # Coordenadas deste nó
        x,y = coord[no,:]

        # Grava em X e Y
        X[i] = x
        Y[i] = y

    end

    # Retorns nós e coordenadas
    return nos, X, Y
end
