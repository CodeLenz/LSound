 
#
# Calcula as matrizes Ke e Me para um elemento 
#
function KMe_tri3(iρ,iκ,X)

  # Mapeamento para facilitar a notação
  x1,x2,x3 = X[:,1]
  y1,y2,y3 = X[:,2]

  # Termo em comum para a rigidez
  comum = ((2*x2-2*x1)*y3+(2*x1-2*x3)*y2+(2*x3-2*x2)*y1)

  # Termos da rigidez
  k11 = (y3^2-2*y2*y3+y2^2+x3^2-2*x2*x3+x2^2)/comum
  k12 = -(y3^2+(-y2-y1)*y3+y1*y2+x3^2+(-x2-x1)*x3+x1*x2)/comum
  k13 = ((y2-y1)*y3-y2^2+y1*y2+(x2-x1)*x3-x2^2+x1*x2)/comum
  k22 = (y3^2-2*y1*y3+y1^2+x3^2-2*x1*x3+x1^2)/comum
  k23 = -((y2-y1)*y3-y1*y2+y1^2+(x2-x1)*x3-x1*x2+x1^2)/comum
  k33 = (y2^2-2*y1*y2+y1^2+x2^2-2*x1*x2+x1^2)/comum
  Ke =  iρ * @SMatrix [k11 k12 k13 ; k12 k22 k23 ; k13 k23 k33]

  # Termos da massa
  m  = iκ *((x2-x1)*y3+(x1-x3)*y2+(x3-x2)*y1)/12
  Me = @SMatrix [m m/2 m/2 ; m/2 m m/2 ; m/2 m/2 m]
    
  return Ke, Me

end
  
# ===================================================================================
# Force vector for a bi3 element 
# local (normal) surface load.
#
function Edge_load_local_tri3(edge,qn,X::Matrix)

  # As we assume cte load
  # and the element is linear
  # Basic test
  edge in 1:3 || throw("Map_load_local_tri3::Invalid edge")

  F1 = 1.0
  F2 = 1.0
  F3 = 1.0

  if edge==1

    dx = X[2,1]-X[1,1]
    dy = X[2,2]-X[1,2]
    F3 = 0.0

  elseif edge==2

    # (2)->(3)
    dx = X[3,1]-X[2,1]
    dy = X[3,2]-X[2,2]
    F1 = 0.0
            
  else

    # (3)->(1)
    dx = X[3,1]-X[1,1]
    dy = X[3,2]-X[1,2]
    F2 = 0.0
    
  end

  # Comprimento da aresta
  L  = sqrt(dx^2 + dy^2)

  # Resposta
  F   = (L/2)*qn*[F1;F2;F3]
  
  # Return F
  return F

end

# ===================================================================================
# Damping matrix Ce
#
function Damping_local_tri3(edge,damp,X)
  
  # Initialize the damping matrix
  C = zeros(3,3)

  # face 1 (12)
  if edge==1
     cte1 = damp*sqrt(X[2,2]^2-2*X[1,2]*X[2,2]+X[1,2]^2+X[2,1]^2-2*X[1,1]*X[2,1]+X[1,1]^2)
     C[1,1] = cte1/3
     C[1,2] = cte1/6
     C[2,1] = C[1,2]
     C[2,2] = cte1/3
  
  elseif edge==2

  # face 2 (23)
     cte2 = damp*sqrt(X[3,2]^2-2*X[2,2]*X[3,2]+X[2,2]^2+X[3,1]^2-2*X[2,1]*X[3,1]+X[2,1]^2)
     C[2,2] = cte2/3
     C[2,3] = cte2/6
     C[3,2] = C[2,3]
     C[3,3] = cte2/3
  
  else
 
    # face 3 (13)
    cte3 = damp*sqrt(X[3,2]^2-2*X[1,2]*X[3,2]+X[1,2]^2+X[3,1]^2-2*X[1,1]*X[3,1]+X[1,1]^2)
    C[1,1] = cte3/3
    C[1,3] = cte3/6
    C[3,1] = C[1,3]
    C[3,3] = cte3/3

  end
  
  # Return C
  return C

end

# ===================================================================================
# Calcula a área do elemento
#
function Area_tri3(X::Matrix)

    # Monta a matriz para o cálculo do determinante
    MA = @SMatrix    [1 X[1,1] X[1,2] ;
                      1 X[2,1] X[2,2] ;
                      1 X[3,1] X[3,2] ]

    # Retorna o determinante
    0.5*det(MA)

end