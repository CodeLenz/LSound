#
#
function Gauss(x,y,cx=0.5,cy=0.5,b=0.2)
    exp(-(x-cx)^2/b^2)*exp(-(y-cy)^2/b^2)
end

function Zero(x,y)
    0.0
end


# Degrau com centro em (xc,yx)
# base b e altura h
function Degrau(x,y,xc,yc,b,h)

    saida = 0.0
    if (xc-b<=x<=xc+b) && (yc-h<=y<=yc+h)
        saida = 1.0
    end

    return saida

end


#
function Applica_U0(nn,coord,funcao::Function)

    # Inicializa o U0
    U0 = zeros(nn)

    # Loop pelos nos
    for no in axes(coord,1)

        # Coordenadas do nó
        x,y = coord[no,:]

        # Calcula a amplitude na coordenada
        U0[no] = funcao(x,y)
        
    end #no

    return U0

end