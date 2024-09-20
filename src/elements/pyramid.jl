#
# Ainda falta testar!

# Devolve a matriz [N] para um  ponto r,s,t
#
function Matriz_N_pyr5(r,s,t)

    N1 = - (1/8)*(((r-1)*s-r+1)*t+(1-r)*s+r-1)
    N2 =   (1/8)*(((r+1)*s-r-1)*t+(-r-1)*s+r+1)
    N3 = - (1/8)*(((r+1)*s+r+1)*t+(-r-1)*s-r-1)
    N4 =   (1/8)*(((r-1)*s+r-1)*t+(1-r)*s-r+1)
    N5 =   (1/2)*(t+1)

    return @SMatrix [N1 N2 N3 N4 N5]
  
  end
  
  #
  # Devolve as derivadas das funções de interpolação
  # em um ponto r,s,t
  #
  function dNrs_pyr5(r,s,t)
      # Derivada das funções de interpolação em relação
      # a r, s e t
      #                    
      dNr = SMatrix{5,1}([ -(((s-1)*(t-1))/8) ; ((s-1)*(t-1))/8 ; -(((r+1)*(t-1))/8) ;
                            ((s+1)*(t-1))/8   ;  0 ])

      dNs = SMatrix{5,1}([ -(((r-1)*(t-1))/8) ; ((r+1)*(t-1))/8 ; -(((r+1)*(t-1))/8) ; 
                            ((r-1)*(t-1))/8  ;  0 ])

      dNt = SMatrix{5,1}([ -(((r-1)*(s-1))/8) ; ((r+1)*(s-1))/8 ; -(((r+1)*(s-1))/8) ; 
                            ((r-1)*(s+1))/8  ; 1/2 ])
      
      return dNr,  dNs, dNt
  end
   
  #
  # Calcula a matriz Jacobiana do elemento
  #
  function Jacobiana_pyr5(r,s,t,X::Array)
  
      # Derivadas das funções de interpolação
      # em relação a    r e s
      dNr, dNs, dNt = dNrs_pyr5(r,s,t)
  
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
  function Matriz_B_pyr5(r,s,t,X::Array)
  
      # Derivadas das funções de interpolação
      # em relação a r,s,t
      dNr, dNs, dNt = dNrs_pyr5(r,s,t)
  
      # Calcula a matriz Jacobiana no ponto r,s,t
      J = Jacobiana_pyr5(r,s,t,X)
  
      # Inicializa a matriz B
      B = @MMatrix zeros(3,8)
  
      # Inverte a J
      iJ = inv(J)
  
      # Loop pelas colunas de B
      for i=1:5
  
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
  function KMe_pyr5(c,X)
  
      # Aloca as matrizes
      Ke = @MMatrix zeros(5,5)
      Me = @MMatrix zeros(5,5)
  
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
                B, dJ = Matriz_B_pyr5(r,s,t,X)
    
                # Calcula N(r,s,t)
                N = Matriz_N_pyr5(r,s,t) 
    
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
# 1) 1 2 3 4 ; <-- Esta face permanece
# 2) 2 3 5 ;   <-- Conferir estas faces 
# 3) 3 4 5 ;
# 4) 4 1 5 ;
# 5) 1 2 5 ;
# 
#
# ######################################################
#
# ######################################################
function Map_face_pyr5(face,ζ,η,X)

    # Basic test
    face in 1:5 || throw("Map_face_hex8::Invalid face")

    # Test for each case 
    if face==1 

        v12 = @SVector [X[2,1] - X[1,1] ;  X[2,2] - X[1,2]; X[2,3] - X[1,3]]
        v14 = @SVector [X[4,1] - X[1,1] ;  X[4,2] - X[1,2]; X[4,3] - X[1,3]]
        A1 = 0.5*norm(cross(v12,v14))

        v23 = @SVector [X[3,1] - X[2,1] ;  X[3,2] - X[2,2]; X[3,3] - X[2,3]]
        v24 = @SVector [X[4,1] - X[2,1] ;  X[4,2] - X[2,2]; X[4,3] - X[2,3]]
        A2 = 0.5*norm(cross(v23,v24))

        # Determinante do Jacobiano para essa face
        dJ = (A1+A2)/4
        
        # N
        N = Matriz_N_pyr5(ζ,η,-1) 

    elseif face==2

        v23 = @SVector [X[3,1] - X[2,1] ;  X[3,2] - X[2,2]; X[3,3] - X[2,3]]
        v25 = @SVector [X[5,1] - X[2,1] ;  X[5,2] - X[2,2]; X[5,3] - X[2,3]]
        A1 = 0.5*norm(cross(v23,v25))

        # Determinante do Jacobiano para essa face <-- Conferir 
        dJ = A1/2
        
        # N
        N = Matriz_N_pyr5(ζ,η,1) 
                  
    elseif face==3

        v34 = @SVector [X[4,1] - X[3,1] ;  X[4,2] - X[3,2]; X[4,3] - X[3,3]]
        v35 = @SVector [X[5,1] - X[3,1] ;  X[5,2] - X[3,2]; X[5,3] - X[3,3]]
        A1 = 0.5*norm(cross(v34,v35))

        # Determinante do Jacobiano para essa face
        dJ = A1/2
        
        # N
        N = Matriz_N_pyr5(ζ,-1,η) 

      
    elseif face==4

        v41 = @SVector [X[1,1] - X[4,1] ;  X[1,2] - X[4,2]; X[1,3] - X[4,3]]
        v45 = @SVector [X[5,1] - X[4,1] ;  X[5,2] - X[4,2]; X[5,3] - X[4,3]]
        A1 = 0.5*norm(cross(v41,v45))

        # Determinante do Jacobiano para essa face
        dJ = A1/2

        # N
        N = Matriz_N_pyr5(1,ζ,η) 
            
    elseif face==5

        v12 = @SVector [X[2,1] - X[1,1] ;  X[2,2] - X[1,2]; X[2,3] - X[1,3]]
        v15 = @SVector [X[5,1] - X[1,1] ;  X[5,2] - X[1,2]; X[5,3] - X[1,3]]
        A1 = 0.5*norm(cross(v12,v15))
        # Determinante do Jacobiano para essa face
        dJ = A1/2 
         
         # N
         N = Matriz_N_pyr5(ζ,1,η) 

    end

    # Return N and dJ
    return N, dJ

end


#
# Force vector for a hex8 element 
#
function Face_load_local_pyr5(face,qn,X)

    # As we assume cte load
    # and the element is linear
    # we can use one Gauss Point

    # Compute Mappings
    N, dJ = Map_face_pyr5(face,0.0,0.0,X)

    # Sums
    F   = (N')*(qn*dJ)*2.0
  
    # Return F
    return F

end

#
# Damping matrix Ce
#
function Damping_local_pyr5(face,damp,X)

    # As we assume cte damp
    # and the element is linear
    # we can use one Gauss Point

    # Compute Mappings (face's center)
    N, dJ = Map_face_pyr5(face,0.0,0.0,X)

    # Matrix
    C   = (N'*N)*damp*(dJ*2.0)
  
    # Return C
    return C

end