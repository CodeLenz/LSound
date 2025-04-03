#
# Calcula o SPL (Sound Pressure Level) para uma dada distribuição de pressao p

#
function SPL(::Vector,p0)

    # Número de pontos 
    nn = length(p)

    # Calcula a soma das pressões em nodes_target ao quadrado 
    P2 = sum(abs.(p).^2)

    # Média (pelo número de pontos em nodes_target)
    P2avg = P2 / nn

    # Como estamos calculando as integrais e as médias por nós, ao invés de 
    # efetivamente considerar comprimentos e áreas, podemos simplificar o P2ref
    P2ref = p0^2

    # Calcula a medida logaritmica
    return 10*log10(P2avg/P2ref)

end

# Média dos SPL em cada uma das frequências consideradas
#
# p0 é a pressão de referência (20μ Pa)
#
function Objetivo(target::Matrix, nodes_target::Vector p0=20E-6)

    # Incializa a soma
    soma = 0.0

    # Para cada linha em target, calcula o SPL
    # Cada coluna em target é p[nodes_target] para uma frequência específica
    for coluna in eachcol(target)
        soma = soma + SPL(coluna[nodes_target],p0)
    end

    # Retorna o valor médio (média pelo número de frequências)
    return soma / size(target,2)


end

