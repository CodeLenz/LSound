#
# Devolve a matriz [N] para um  ponto r,s,t
#
function Matriz_N_hex8(r,s,t)

    N1 = (1/8)*(1-r)*(1-s)*(1-t)
    N2 = (1/8)*(1+r)*(1-s)*(1-t)
    N3 = (1/8)*(1+r)*(1+s)*(1-t)
    N4 = (1/8)*(1-r)*(1+s)*(1-t)
    N5 = (1/8)*(1-r)*(1-s)*(1+t)
    N6 = (1/8)*(1+r)*(1-s)*(1+t)
    N7 = (1/8)*(1+r)*(1+s)*(1+t)
    N8 = (1/8)*(1-r)*(1+s)*(1+t)
  

    return @SMatrix [N1 N2 N3 N4 N5 N6 N7 N8]
  
  end
  
  #
  # Devolve as derivadas das funções de interpolação
  # em um ponto r,s,t
  #
  function dNrs_hex8(r,s,t)
      # Derivada das funções de interpolação em relação
      # a r, s e t
      #                    
      dNr = SMatrix{8,1}([ -(((s-1)*t-s+1)/8) ; ((s-1)*t-s+1)/8; -(((s+1)*t-s-1)/8); ((s+1)*t-s-1)/8 ;
                            ((s-1)*t+s-1)/8;  -(((s-1)*t+s-1)/8); ((s+1)*t+s+1)/8 ; -(((s+1)*t+s+1)/8)  ])
      dNs = SMatrix{8,1}([ -(((r-1)*t-r+1)/8) ; ((r+1)*t-r-1)/8 ; -(((r+1)*t-r-1)/8) ;
                            ((r-1)*t-r+1)/8 ; ((r-1)*t+r-1)/8 ; -(((r+1)*t+r+1)/8) ;
                             ((r+1)*t+r+1)/8 ; -(((r-1)*t+r-1)/8) ])
      dNt = SMatrix{8,1}([ -(((r-1)*s-r+1)/8) ; ((r+1)*s-r-1)/8 ; -(((r+1)*s+r+1)/8) ; 
                             ((r-1)*s+r-1)/8; ((r-1)*s-r+1)/8; -(((r+1)*s-r-1)/8); 
                             ((r+1)*s+r+1)/8 ; -(((r-1)*s+r-1)/8) ])
      
      return dNr,  dNs, dNt
  end
   
  #
  # Calcula a matriz Jacobiana do elemento
  #
  function Jacobiana_hex8(r,s,t,X::Array)
  
      # Derivadas das funções de interpolação
      # em relação a    r e s
      dNr, dNs, dNt = dNrs_hex8(r,s,t)
  
      # Inicializa a matriz J
      J = @MMatrix zeros(3,3)
  
      # Loop pelos somatórios
      for i=1:8
          J[1,1] += dNr[i]*X[i,1]
          J[1,2] += dNr[i]*X[i,2]
          J[1,3] += dNr[i]*X[i,3]

          J[2,1] += dNs[i]*X[i,1]
          J[2,2] += dNs[i]*X[i,2]
          J[2,3] += dNs[i]*X[i,3]
          
          J[3,1] += dNt[i]*X[i,1]
          J[3,2] += dNt[i]*X[i,2]
          J[3,3] += dNt[i]*X[i,3]
          
      end
   
      # Devolve a matriz Jacobiana para o elemento
      # no ponto r,s,t
      return J
  
  end
  
  #
  # Monta a matriz B de um elemento na posiçao r,s
  #
  function Matriz_B_hex8(r,s,t,X::Array)
  
      # Derivadas das funções de interpolação
      # em relação a r,s,t
      dNr, dNs, dNt = dNrs_hex8(r,s,t)
  
      # Calcula a matriz Jacobiana no ponto r,s,t
      J = Jacobiana_hex8(r,s,t,X)
  
      # Inicializa a matriz B
      B = @MMatrix zeros(3,8)
  
      # Inverte a J
      iJ = inv(J)
  
      # Loop pelas colunas de B
      for i=1:8
  
          # Corrige as derivadas de rs para xy
          dNxy = iJ*[dNr[i];dNs[i];dNt[i]]
   
          # Posiciona na coluna
          B[1,i] = dNxy[1]
          B[2,i] = dNxy[2]
          B[3,i] = dNxy[3]
  
      end
  
      # Devolve B e o dJ
      return B, det(J)
  
  end
  
 
  #
  # Calcula as matrizes Ke e Me para um elemento 
  #
  function KMe_hex8(ele,c,X)
  
      # Aloca as matrizes
      Ke = @MMatrix zeros(8,8)
      Me = @MMatrix zeros(8,8)
  
      # Integração por quadratura de Gauss-Legendre
      pg = (1/sqrt(3))*[-1;1]
      wg = ones(2)
  
      @inbounds for i=1:2
          # Ponto e peso nesta dimensão
          r = pg[i]
          wr = wg[i]
  
          @inbounds for j=1:2
              # Ponto e peso nesta dimensão
              s = pg[j]
              ws = wg[j]

              @inbounds for k=1:2
                # Ponto e peso nesta dimensão
                t = pg[k]
                wt = wg[k]

                # Calcula DJ e B 
                B, dJ = Matriz_B_hex8(r,s,t,X)
    
                # Calcula N(r,s,t)
                N = Matriz_N_hex8(r,s,t) 
    
                # Somatórios
                Me = Me + N'*N*dJ
                Ke = Ke + B'*B*dJ

              end  #k 
          end #j
      end #i
  
      return Ke, (1/c^2)*Me
  
  end
  

#=

  ############## MODIFICAR DAQUI PARA BAIXO #####################

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
function Edge_load_local_bi4(edge,qn,X,Y)

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

#
# Damping matrix Ce
#
function Damping_local_bi4(edge,damp,X,Y)

    # As we assume cte damp
    # and the element is linear
    # we can use one Gauss Point

    # Compute Mappings (edge's center)
    n, t, dJ, N = Map_edge_bi4(edge,0.0,X,Y)

    # Matrix
    C   = (N'*N)*damp*(dJ*2.0)
  
    # Return C
    return C

end

=#