#
# Devolve a matriz [N] para um  ponto r,s,t
# (matriz com as funções de interpolação para este elemento)
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

# ===================================================================================
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

# ===================================================================================
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
 
# ===================================================================================
# Monta a matriz B de um elemento na posição r,s
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

# ===================================================================================
# Calcula as matrizes Ke e Me para um elemento 
#
function KMe_hex8(iρ,iκ,X)
      
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
  
    return  iρ*Ke, iκ*Me
  
end

# ===================================================================================
# Faces
# 1) 1 2 3 4 ;
# 2) 5 6 7 8 ;
# 3) 1 2 6 5 ;
# 4) 2 3 7 6 ;
# 5) 4 3 7 8 ;
# 6) 1 4 8 5 ;
# ===================================================================================
function Map_face_hex8(face,ζ,η,X)

    # Basic test
    face in 1:6 || throw("Map_face_hex8::Invalid face")

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
        N = Matriz_N_hex8(ζ,η,-1) 
        quad
    elseif face==2

        v56 = @SVector [X[6,1] - X[5,1] ;  X[6,2] - X[5,2]; X[6,3] - X[5,3]]
        v58 = @SVector [X[8,1] - X[5,1] ;  X[8,2] - X[5,2]; X[8,3] - X[5,3]]
        A1 = 0.5*norm(cross(v56,v58))

        v67 = @SVector [X[7,1] - X[6,1] ;  X[7,2] - X[6,2]; X[7,3] - X[6,3]]
        v68 = @SVector [X[8,1] - X[6,1] ;  X[8,2] - X[6,2]; X[8,3] - X[6,3]]
        A2 = 0.5*norm(cross(v67,v68))

        # Determinante do Jacobiano para essa face
        dJ = (A1+A2)/4
        
        # N
        N = Matriz_N_hex8(ζ,η,1) 
                  
    elseif face==3

        v12 = @SVector [X[2,1] - X[1,1] ;  X[2,2] - X[1,2]; X[2,3] - X[1,3]]
        v15 = @SVector [X[5,1] - X[1,1] ;  X[5,2] - X[1,2]; X[5,3] - X[1,3]]
        A1 = 0.5*norm(cross(v12,v15))

        v26 = @SVector [X[6,1] - X[2,1] ;  X[6,2] - X[2,2]; X[6,3] - X[2,3]]
        v25 = @SVector [X[5,1] - X[2,1] ;  X[5,2] - X[2,2]; X[5,3] - X[2,3]]
        A2 = 0.5*norm(cross(v26,v25))

        # Determinante do Jacobiano para essa face
        dJ = (A1+A2)/4
        
        # N
        N = Matriz_N_hex8(ζ,-1,η) 

    elseif face==4

        v23 = @SVector [X[3,1] - X[2,1] ;  X[3,2] - X[2,2]; X[3,3] - X[2,3]]
        v26 = @SVector [X[6,1] - X[2,1] ;  X[6,2] - X[2,2]; X[6,3] - X[2,3]]
        A1 = 0.5*norm(cross(v23,v26))

        v37 = @SVector [X[7,1] - X[3,1] ;  X[7,2] - X[3,2]; X[7,3] - X[3,3]]
        v36 = @SVector [X[6,1] - X[3,1] ;  X[6,2] - X[3,2]; X[6,3] - X[3,3]]
        A2 = 0.5*norm(cross(v37,v36))

        # Determinante do Jacobiano para essa face
        dJ = (A1+A2)/4

        # N
        N = Matriz_N_hex8(1,ζ,η) 
            
    elseif face==5

        v43 = @SVector [X[3,1] - X[4,1] ;  X[3,2] - X[4,2]; X[3,3] - X[4,3]]
        v48 = @SVector [X[8,1] - X[4,1] ;  X[8,2] - X[4,2]; X[8,3] - X[4,3]]
        A1 = 0.5*norm(cross(v43,v48))

        v37 = @SVector [X[7,1] - X[3,1] ;  X[7,2] - X[3,2]; X[7,3] - X[3,3]]
        v38 = @SVector [X[8,1] - X[3,1] ;  X[8,2] - X[3,2]; X[8,3] - X[3,3]]
        A2 = 0.5*norm(cross(v37,v38))

        # Determinante do Jacobiano para essa face
        dJ = (A1+A2)/4
         
        # N
        N = Matriz_N_hex8(ζ,1,η) 
        
    else 

        v14 = @SVector [X[4,1] - X[1,1] ;  X[4,2] - X[1,2]; X[4,3] - X[1,3]]
        v15 = @SVector [X[5,1] - X[1,1] ;  X[5,2] - X[1,2]; X[5,3] - X[1,3]]
        A1 = 0.5*norm(cross(v14,v15))

        v48 = @SVector [X[8,1] - X[4,1] ;  X[8,2] - X[4,2]; X[8,3] - X[4,3]]
        v45 = @SVector [X[5,1] - X[4,1] ;  X[5,2] - X[4,2]; X[5,3] - X[4,3]]
        A2 = 0.5*norm(cross(v48,v45))

        # Determinante do Jacobiano para essa face
        dJ = (A1+A2)/4
         
        # N
        N = Matriz_N_hex8(-1,ζ,η) 

    end

    # Return N and dJ
    return N, dJ

end

# ===================================================================================
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

# ===================================================================================
# Damping matrix Ce
#
function Damping_local_hex8(face,damp,X)

    # As we assume cte damp
    # and the element is linear
    # we can use one Gauss Point

    # Compute Mappings (face's center)
    N, dJ = Map_face_hex8(face,0.0,0.0,X)

    # Matrix
    C   = (N'*N)*damp*(dJ*2.0)
  
    # Return C
    return C

end

# ===================================================================================
# Calcula a volume do elemento
#
function Volume_hex8(X::Matrix)

    # Inicializa o somatório do volume
    V = 0.0 

    # Integração por quadratura de Gauss-Legendre
    pg = (1/sqrt(3))*[-1;1]
    
    for i=1:2
        # Ponto nesta dimensão
        r = pg[i]
        
        for j=1:2
            # Ponto nesta dimensão
            s = pg[j]
            
            for k=1:2
                #Ponto nesta dimensão
                t = pg[k]

                # Calcula a matriz Jacobiana no ponto r,s
                J = Jacobiana_hex8(r,s,t,X)

                # Adiciona o determinante do Jacobiano 
                V = V + det(J)

            end # k

        end # j

    end #i

    # Retorna o volume
    return V

end