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
  function KMe_hex8(c,X)
  
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
  

#
# Faces
# 1) 1 2 3 4 ;
# 2) 5 6 7 8 ;
# 3) 1 2 6 5 ;
# 4) 2 3 7 6 ;
# 5) 4 3 7 8 ;
# 6) 1 4 8 5 
#
#
function Map_face_hex8(face,ζ,η,X)

    # Basic test
    edge in 1:6 || throw("Map_face_hex8::Invalid face")

    # N at each node of the face
    Na = (1/4)*(1-ζ)*(1-η)
    Nb = (1/4)*(1+ζ)*(1-η)
    Nc = (1/4)*(1+ζ)*(1+η)
    Nd = (1/4)*(1-ζ)*(1+η)
    
    # N's
    N1 = 0.0
    N2 = 0.0
    N3 = 0.0
    N4 = 0.0
    N5 = 0.0
    N6 = 0.0
    N7 = 0.0
    N8 = 0.0

    # Test for each case
    if face==1 

        # Jacobian matrix 
        J_ = Jacobiana_hex8(ζ,η,-1,X)

        # We use just the block 1:2, 1:2
        J = J_[1:2,1:2]
        
        # Map the N's
        N1 = Na
        N2 = Nb
        N3 = Nc
        N4 = Nd

    elseif face==2

        # Jacobian matrix 
        J_ = Jacobiana_hex8(ζ,η,1,X)

        # We use just the block 1:2, 1:2
        J = J_[1:2,1:2]
        
        # Map the N's
        N5 = Na
        N6 = Nb
        N7 = Nc
        N8 = Nd
                  
    elseif face==3

        # Jacobian matrix 
        J_ = Jacobiana_hex8(ζ,-1,η,X)

        # We use just positions 1 3
        p = [1;3]
        J = J_[p,p]
        
        # Map the N's
        N1 = Na
        N2 = Nb
        N6 = Nc
        N5 = Nd
      
    elseif face==4

        # Jacobian matrix 
        J_ = Jacobiana_hex8(1,ζ,η,X)

        # We use just positions 2:3
        J = J_[2:3,2:3]
        
        # Map the N's
        N2 = Na
        N3 = Nb
        N7 = Nc
        N6 = Nd
    
    elseif face==5

         # Jacobian matrix 
         J_ = Jacobiana_hex8(ζ,1,η,X)

         # We use just positions 1 and 3
         p = [1;3]
         J = J_[p,p]
         
         # Map the N's
         N4 = Na
         N3 = Nb
         N7 = Nc
         N8 = Nd

    else 

         # Jacobian matrix 
         J_ = Jacobiana_hex8(-1,ζ,η,X)

         # We use just positions 2:3
         J = J_[2:3,2:3]
         
         # Map the N's
         N1 = Na
         N4 = Nb
         N8 = Nc
         N5 = Nd


    end

    # Determinant is
    dJ = det(J)

    # Matriz N
    N = [N1 N2 N3 N4 N5 N6 N7 N8]

    # Return N and dJ
    return N, dJ

end


#
# Force vector for a hex8 element 
#
function Face_load_local_hex8(face,qn,X)

    # As we assume cte load
    # and the element is linear
    # we can use one Gauss Point

    # Compute Mappings
    N, dJ = Map_face_hex8(face,0.0,0.0,X)

    # Sums
    F   = (N')*(qn*dJ)*2.0
  
    # Return F
    return F

end

#
# Damping matrix Ce
#
function Damping_local_hex8(edge,damp,X)

    # As we assume cte damp
    # and the element is linear
    # we can use one Gauss Point

    # Compute Mappings (edge's center)
    dJ, N = Map_face_hex8(face,0.0,0.0,X)

    # Matrix
    C   = (N'*N)*damp*(dJ*2.0)
  
    # Return C
    return C

end