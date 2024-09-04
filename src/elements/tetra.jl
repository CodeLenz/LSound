#
# Devolve a matriz [N] para um  ponto r,s,t
#
  function Matriz_N_tet4(r,s,t)

    N1 = -((((r-1)*s-r+1)*t+(1-r)*s+r-1)/8)
    N2 = (r*s*t)/8+(s*t)/8-(r*t)/8-t/8-(r*s)/8-s/8+r/8+1/8
    N3 = -((s*t)/4)-t/4+s/4+1/4
    N4 = t/2+1/2

    return @SMatrix [N1 N2 N3 N4]
  
  end
  
  #
  # Volume do tetraedro
  #
  function Vol_tet4(X::Array)

      # Coordinates
      x1,x2,x3,x4 = X[:,1] 
      y1,y2,y3,y4 = X[:,2]
      z1,z2,z3,z4 = X[:,3]

      # Static matrix
      V = @SMatrix [1 x1 y1 z1 ; 
                    1 x2 y2 z2 ; 
                    1 x3 y3 z3 ;
                    1 x4 y4 z4]

      return det(V)/6

  end


  #
  # Determinante da matriz Jacobiana
  # 

  #
  # Calcula a matriz Jacobiana do elemento
  #
  function Jacobiana_tetra4(r,s,t,X::Array)
  
       # Coordinates
       x1,x2,x3,x4 = X[:,1] 
       y1,y2,y3,y4 = X[:,2]
       z1,z2,z3,z4 = X[:,3]

      # Posições da matriz
      J11 = (((s-1)*t-s+1)*x2+((1-s)*t+s-1)*x1)/8
      J12 = (((s-1)*t-s+1)*y2+((1-s)*t+s-1)*y1)/8
      J13 = (((s-1)*t-s+1)*z2+((1-s)*t+s-1)*z1)/8

      J21 = -(((2*t-2)*x3+((-r-1)*t+r+1)*x2+((r-1)*t-r+1)*x1)/8)
      J22 = -(((2*t-2)*y3+((-r-1)*t+r+1)*y2+((r-1)*t-r+1)*y1)/8)
      J23 = -(((2*t-2)*z3+((-r-1)*t+r+1)*z2+((r-1)*t-r+1)*z1)/8)

      J31 = (4*x4+(-(2*s)-2)*x3+((r+1)*s-r-1)*x2+((1-r)*s+r-1)*x1)/8
      J32 = (4*y4+(-(2*s)-2)*y3+((r+1)*s-r-1)*y2+((1-r)*s+r-1)*y1)/8
      J33 = (4*z4+(-(2*s)-2)*z3+((r+1)*s-r-1)*z2+((1-r)*s+r-1)*z1)/8
   
      # Devolve a matriz Jacobiana para o elemento no ponto r,s,t
      return @SMatrix [J11 J12 J13 ; J21 J22 J23 ; J31 J32 J33]
  
  end
  
  #
  # Monta a matriz B de um elemento tetraedrico linear
  #
  function Matriz_B_tet4(X::Array)
  
    # Coordinates
    x1,x2,x3,x4 = X[:,1] 
    y1,y2,y3,y4 = X[:,2]
    z1,z2,z3,z4 = X[:,3]

    # Common term
    cte = (((x2-x1)*y3+(x1-x3)*y2+(x3-x2)*y1)*z4+((x1-x2)*y4+(x4-x1)*y2+(x2-x4)*y1)*z3+
           ((x3-x1)*y4+(x1-x4)*y3+(x4-x3)*y1)*z2+((x2-x3)*y4+(x4-x2)*y3+(x3-x4)*y2)*z1)

    B11 = -((y3-y2)*z4+(y2-y4)*z3+(y4-y3)*z2)/cte
    B12 = ((y3-y1)*z4+(y1-y4)*z3+(y4-y3)*z1)/cte
    B13 = -(((y2-y1)*z4+(y1-y4)*z2+(y4-y2)*z1))/cte       
    B14 = ((y2-y1)*z3+(y1-y3)*z2+(y3-y2)*z1)/cte

    B21 = ((x3-x2)*z4+(x2-x4)*z3+(x4-x3)*z2)/cte
    B22 = -((x3-x1)*z4+(x1-x4)*z3+(x4-x3)*z1)/cte
    B23 = ((x2-x1)*z4+(x1-x4)*z2+(x4-x2)*z1)/cte
    B24 = -((x2-x1)*z3+(x1-x3)*z2+(x3-x2)*z1)/cte

    B31 = -((x3-x2)*y4+(x2-x4)*y3+(x4-x3)*y2)/cte
    B32 = ((x3-x1)*y4+(x1-x4)*y3+(x4-x3)*y1)/cte
    B33 = -((x2-x1)*y4+(x1-x4)*y2+(x4-x2)*y1)/cte
    B34 = ((x2-x1)*y3+(x1-x3)*y2+(x3-x2)*y1)/cte

    return @SMatrix [B11 B12 B13 B14 ; B21 B22 B23 B24 ; B31 B32 B33 B34] 

  end
  
 
  #
  # Calcula as matrizes Ke e Me para um elemento 
  #
  function KMe_tet4(c,X)

      # Coordinates
      x1,x2,x3,x4 = X[:,1] 
      y1,y2,y3,y4 = X[:,2]
      z1,z2,z3,z4 = X[:,3]

      # A matriz de rigidez pode ser calculada de maneira bem rápida,
      # pois B é cte (não depende de r,s,t)
  
      # Volume do tetraedro
      V = Vol_tet4(X)

      # Calcula a matriz B (cte)
      B = Matriz_B_tet4(X)

      # Rigidez
      Ke = transpose(B)*B*V

      # Calculo da matriz da massa
      cte = (((x2-x1)*y3+(x1-x3)*y2+(x3-x2)*y1)*z4+((x1-x2)*y4+(x4-x1)*y2+
             (x2-x4)*y1)*z3+((x3-x1)*y4+(x1-x4)*y3+(x4-x3)*y1)*z2+
             ((x2-x3)*y4+(x4-x2)*y3+(x3-x4)*y2)*z1) / c^2

      cte1 = cte / 60
      cte2 = cte1 / 2       
     
      Me = @SMatrix [cte1 cte2 cte2 cte2 ;
                     cte2 cte1 cte2 cte2 ; 
                     cte2 cte2 cte1 cte2 ; 
                     cte2 cte2 cte2 cte1 ]

      return Ke, Me
  
  end
  
