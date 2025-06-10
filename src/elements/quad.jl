#
# Devolve a matriz [N] para um  ponto r,s
# (matriz com as funções de interpolação para este elemento)
#
function Matriz_N_bi4(r,s)

    N1 = (1/4)*(1-r)*(1-s)
    N2 = (1/4)*(1+r)*(1-s)
    N3 = (1/4)*(1+r)*(1+s)
    N4 = (1/4)*(1-r)*(1+s)
  
    return @SMatrix [N1 N2 N3 N4]
  
end
  
# ===================================================================================
# Devolve as derivadas das funções de interpolação
# em um ponto r,s
#
function dNrs_bi4(r,s)
    # Deriva das funções de interpolação em relação
    # a 'r' e 's'
    #                           N1     N2     N3       N4
    dNr = SMatrix{4,1}((1/4)*[-(1-s);  (1-s); (1+s); -(1+s)])
    dNs = SMatrix{4,1}((1/4)*[-(1-r); -(1+r); (1+r);  (1-r)])
      
    return dNr,  dNs
end

# ===================================================================================
# Calcula a matriz Jacobiana do elemento
#
function Jacobiana_bi4(r,s,X::Array)
  
    # Derivadas das funções de interpolação
    # em relação a    r e s
    dNr, dNs = dNrs_bi4(r,s)
  
    # Inicializa a matriz J
    J = @MMatrix zeros(2,2)
  
    # Loop pelos somatórios
    for i=1:4
        J[1,1] += dNr[i]*X[i,1]
        J[1,2] += dNr[i]*X[i,2]

        J[2,1] += dNs[i]*X[i,1]
        J[2,2] += dNs[i]*X[i,2]
    end
   
    # Devolve a matriz Jacobiana para o elemento
    # no ponto r,s
    return J
  
end

# ===================================================================================
# Monta a matriz B de um elemento na posiçao r,s
#
function Matriz_B_bi4(r,s,X::Array)
  
    # Derivadas das funções de interpolação
    # em relação a    'r' e 's'
    dNr, dNs = dNrs_bi4(r,s)
  
    # Calcula a matriz Jacobiana no ponto r,s
    J = Jacobiana_bi4(r,s,X)
  
    # Inicializa a matriz B
    B = @MMatrix zeros(2,4)
  
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
  
# ===================================================================================
# Calcula as matrizes Ke e Me para um elemento 
#
function KMe_bi4(iρ,iκ,X)

    # Aloca as matrizes
    Ke = @MMatrix zeros(4,4)
    Me = @MMatrix zeros(4,4)
  
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
  
            # Calcula DJ e B 
            B, dJ = Matriz_B_bi4(r,s,X)
  
            # Calcula N(r,s)
            N = Matriz_N_bi4(r,s) 
  
            # Somatórios
            Me = Me + N'*N*dJ
            Ke = Ke + B'*B*dJ
  
        end
    end
  
    return iρ*Ke, iκ*Me
  
end

# ===================================================================================
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
# ===================================================================================
function Map_edge_bi4(edge,ζ,X)

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
        t1 = X[1,1]*dN1r + X[2,1]*dN2r
        t2 = X[1,2]*dN1r + X[2,2]*dN2r

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
        t1 = X[2,1]*dN2s + X[3,1]*dN3s
        t2 = X[2,2]*dN2s + X[3,2]*dN3s

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
        t1 = X[4,1]*dN4r + X[3,1]*dN3r
        t2 = X[4,2]*dN4r + X[3,2]*dN3r
    
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
        t1 = X[4,1]*dN4s + X[1,1]*dN1s
        t2 = X[4,2]*dN4s + X[1,2]*dN1s

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

# ===================================================================================
# Force vector for a bi4 element 
# local (normal) surface load.
#
function Edge_load_local_bi4(edge,qn,X)

    # As we assume cte load
    # and the element is linear
    # we can use one Gauss Point

    # Compute Mappings
    n, t, dJ, N = Map_edge_bi4(edge,0.0,X)

    # Sums
    F   = (N')*(qn*dJ)*2.0
  
    # Return F
    return F

end

# ===================================================================================
# Damping matrix Ce
#
function Damping_local_bi4(edge,damp,X)

    # As we assume cte damp
    # and the element is linear
    # we can use one Gauss Point

    # Compute Mappings (edge's center)
    n, t, dJ, N = Map_edge_bi4(edge,0.0,X)

    # Matrix
    C   = (N'*N)*damp*(dJ*2.0)
  
    # Return C
    return C

end

# ===================================================================================
# Calcula a área do elemento
#
function Area_bi4(X::Matrix)

    # Inicializa o somatório da área
    A = 0.0 

    # Integração por quadratura de Gauss-Legendre
    pg = (1/sqrt(3))*[-1;1]
    
    for i=1:2
        # Ponto nesta dimensão
        r = pg[i]
        
        for j=1:2
            # Ponto nesta dimensão
            s = pg[j]
        
            # Calcula a matriz Jacobiana no ponto r,s
            J = Jacobiana_bi4(r,s,X)

            # Adiciona o determinante do Jacobiano para o ponto
            A = A + det(J)

        end # j

    end #i

    # Retorna a área
    return A

end