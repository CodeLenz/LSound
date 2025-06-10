#
#  Devolve um vetor com os n nós do elemento
#  e uma matriz com as coordenadas [X], no formato
#  x1 y1 z1
#  x2 y2 z2
#   ..... 
#  xn yn zn 
#
function Nos_Coordenadas(ele,etype,coord,connect)
  
  # Mapeia o número de nós pelo etype
  nnos = Lgmsh_nodemap()[etype]

  # Descobre os nós do elemento
  nos = connect[ele,3:2+nnos]

  # Aloca as coordenadas
  X = Array{Float64}(undef,nnos,3)
    
  # Descobre as coordenadas de cada nó do elemento
  for i in LinearIndices(nos)

    # Nó
    no = nos[i]

    # Coordenadas deste nó
    x,y,z = coord[no,:]

    # Grava 
    X[i,1] = x
    X[i,2] = y
    X[i,3] = z

  end

  # Retorna nós e coordenadas
  return nos, X

end