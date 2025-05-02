#
# Calcula o SPL (Sound Pressure Level) para uma dada distribuição de pressao p
# equivalente a uma frequência ω_n. O vetor P foi obtido via Sweep

#
function SPLn(P::Vector,p0)

    # Número de pontos 
    nt = length(P)

    # Calcula a soma das pressões em nodes_target ao quadrado 
    P2 = sum((abs.(P)).^2)

    # Média (pelo número de pontos em nodes_target)
    P2avg = P2 / nt

    # Como estamos calculando as integrais e as médias por nós, ao invés de 
    # efetivamente considerar comprimentos e áreas, podemos simplificar o P2ref
    P2ref = p0^2

    # Calcula a medida logaritmica
    return 10*log10(P2avg/P2ref)

end

# Média dos SPL em cada uma das frequências consideradas
#
# MP é uma matriz nn × nf com pressões em diferentes frequências
#
# p0 é a pressão de referência (20μ Pa)
#
function Objetivo(MP::Matrix, nodes_target::Vector, p0=20E-6)

    # Incializa a soma
    soma = 0.0

    # Para cada coluna em MP, calcula o SPL
    # Cada coluna em MP é P[nodes_target] para uma frequência específica
    for P in eachcol(MP)
        soma = soma + SPLn(P[nodes_target],p0)
    end

    # Número de frequências é o número de colunas em target
    Nf = size(target,2)

    # Retorna o valor médio (média pelo número de frequências)
    return soma / Nf


end

