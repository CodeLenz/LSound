#
# Devolve a matriz [N] para um  ponto r,s
#
function Matriz_N(r,s)

    N1 = (1/4)*(1-r)*(1-s)
    N2 = (1/4)*(1+r)*(1-s)
    N3 = (1/4)*(1+r)*(1+s)
    N4 = (1/4)*(1-r)*(1+s)
  
    return [N1 N2 N3 N4]
  
  end
  
  #
  # Devolve as derivadas das funções de interpolação
  # em um ponto r,s
  #
  function dNrs(r,s)
      # Deriva das funções de interpolação em relação
      # a r e s
      #               N1     N2     N3       N4
      dNr = (1/4)*[-(1-s);  (1-s); (1+s); -(1+s)]
      dNs = (1/4)*[-(1-r); -(1+r); (1+r);  (1-r)]
      
      return dNr, dNs
  end
   
  #
  # Calcula a matriz Jacobiana do elemento
  #
  function Jacobiana(r,s,X::Vector,Y::Vector)
  
      # Derivadas das funções de interpolação
      # em relação a    r e s
      dNr, dNs = dNrs(r,s)
  
      # Inicializa a matriz J
      J = zeros(2,2)
  
      # Loop pelos somatórios
      for i=1:4
          J[1,1] += dNr[i]*X[i]
          J[1,2] += dNr[i]*Y[i]
          J[2,1] += dNs[i]*X[i]
          J[2,2] += dNs[i]*Y[i]
      end
   
      # Devolve a matriz Jacobiana para o elemento
      # no ponto r,s
      return J
  
  end
  
  #
  # Monta a matriz B de um elemento na posiçao r,s
  #
  function Matriz_B(r,s,X::Vector,Y::Vector)
  
      # Derivadas das funções de interpolação
      # em relação a    r e s
      dNr, dNs = dNrs(r,s)
  
      # Calcula a matriz Jacobiana no ponto r,s
      J = Jacobiana(r,s,X,Y)
  
      # Inicializa a matriz B
      B = zeros(2,4)
  
      # Inverte a J
      iJ = inv(J)
  
      # Loop pelas colunas de B
      for i=1:4
  
          # Corrige as derivadas de rs para xy
          dNxy = iJ*[dNr[i];dNs[i]]
   
          # Posiciona na coluna
          B[1,i] = dNxy[1]
          B[2,i] = dNxy[2]
  
      end
  
      # Devolve B e o dJ
      return B, det(J)
  
  end
  
  #
  # Devolve os vetores nos, X e Y com os nós 
  # e as coordenadas do elemento
  #
  function Nos_Coordenadas(ele,coord,connect)
  
      # Descobre os nós do elemento
      nos = connect[ele,3:end]
  
      # Aloca vetores de coordenadas
      X = zeros(4)
      Y = zeros(4)
  
      # Descobre as coordenadas x e y de cada nó
      # do elemento
      for i=1:4
  
          # Nó
          no = nos[i]
  
          # Coordenadas deste nó
          x,y = coord[no,:]
  
          # Grava em X e Y
          X[i] = x
          Y[i] = y
  
      end
  
      return nos, X, Y
  end
  
  #
  # Calcula as matrizes Ke e Me para um elemento 
  #
  function KMe(ele,c,coord,connect)
  
      # Aloca as matrizes
      Ke = zeros(4,4)
      Me = zeros(4,4)
  
      # Descobre nos, X e Y para este elemento
      nos, X, Y = Nos_Coordenadas(ele,coord,connect) 
  
      # Integração por quadratura de Gauss-Legendre
      pg = (1/sqrt(3))*[-1;1]
      wg = ones(2)
  
      for i=1:2
          # Ponto e peso nesta dimensão
          r = pg[i]
          wr = wg[i]
  
          for j=1:2
              # Ponto e peso nesta dimensão
              s = pg[j]
              ws = wg[j]
  
              # Calcula DJ e B 
              B, dJ = Matriz_B(r,s,X,Y)
  
              # Calcula N(r,s)
              N = Matriz_N(r,s) 
  
              # Somatórios
              Me = Me + N'*N*dJ
              Ke = Ke + B'*B*dJ
  
          end
      end
  
      return Ke, (1/c^2)*Me, nos
  
  end
  
  

  #
#    Edges, normals and tangents
#
#               n
#          t_ 3 |
#       (4)-------(3)
#    n _ |         | | t
#      4 |         | 2 
#      | |         | _ n 
#      t(1)-------(2)
#           | 1 -t
#           n
#
#
#
function Map_edge_bi4(edge,ζ,X,Y)

    # Basic test
    edge in 1:4 || throw("Map_edge_bi4::Invalid edge")

    # N at each node of the edge
    Na = (1/2)*(1-ζ)
    Nb = (1/2)*(1+ζ)

    # N's
    N1 = 0.0
    N2 = 0.0
    N3 = 0.0
    N4 = 0.0

    # Test for each case
    if edge==1 

        # bottom
        # s = -1
        # r free
        mult = 1.0
        dN1r = -1/2
        dN2r =  1/2
        t1 = X[1]*dN1r + X[2]*dN2r
        t2 = Y[1]*dN1r + Y[2]*dN2r

        # Functions
        N1 = Na
        N2 = Nb
        
    elseif edge==2

        # right
        # r = 1
        # s free
        mult = 1.0
        dN2s = -1/2
        dN3s =  1/2
        t1 = X[2]*dN2s + X[3]*dN3s
        t2 = Y[2]*dN2s + Y[3]*dN3s

        # Functions
        N2 = Na
        N3 = Nb
            
    elseif edge==3

        # top
        # s = 1
        # r free
        mult = -1.0
        dN3r =  1/2
        dN4r = -1/2
        t1 = X[4]*dN4r + X[3]*dN3r
        t2 = Y[4]*dN4r + Y[3]*dN3r
    
        # Functions
        N3 = Nb
        N4 = Na

    else

        # left
        # r = -1
        # s free
        mult = -1.0
        dN4s =  1/2
        dN1s = -1/2
        t1 = X[4]*dN4s + X[1]*dN1s
        t2 = Y[4]*dN4s + Y[1]*dN1s

        # Functions
        N1 = Na
        N4 = Nb
            
    end

    # Determinant of the Jacobian
    dJ = sqrt(t1^2 + t2^2)

    # Normalize and compute t and n
    t = mult*[t1;t2]./dJ
    n = [t[2];-t[1]]

    # Matriz N
    N = [N1 N2 N3 N4]

    # Return n, t, dJ and N
    return n, t, dJ, N

end

#
# Force vector for a bi4 element 
# local (normal) surface load.
#
function Edge_load_local(edge,qn,X,Y)

    # As we assume cte load
    # and the element is linear
    # we can use one Gauss Point

    # Compute Mappings
    n, t, dJ, N = Map_edge_bi4(edge,0.0,X,Y)

    # Sums
    F   = (N')*(qn*dJ)*2.0
  
    # Return F
    return F

end