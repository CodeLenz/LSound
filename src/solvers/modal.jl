function Modal(K,M,livres,nev)

    # Cria as vistas para as posições das matrizes
    KL = K[livres,livres]
    ML = M[livres,livres]

    # Vamos ter que apelar para um solver modal mais porrada
    flag, lamb, X = Solve_Eigen_(KL, ML, nev)

    # Numero efetivo de modos calculados
    nemodos = length(lamb)

    # Avisa se temos a situação em que nem todos os modos
    # retornados são utilizados (por exemplo, uma freq neg)
    # numero efetivo de modos
    if nemodos<nev
        println("Número efetivo de modos ($nemodos) é menor do que o solicitado")
    end

    # Frequências em Hz
    freq =  sqrt.(lamb)./(2*pi)

    # Retorna frequências e modos
    return freq, X  

end

