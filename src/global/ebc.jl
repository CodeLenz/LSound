#
#
# Calcula o vetor de "forças" devido à aplicação das condições de contorno 
# de pressão 
# 
# A matriz K deve ser a matriz sem a aplicação da máscara livres 
#
#
# TODO fase deveria entrar como a parte complexa ?
#
#
function P_pressure(nnos::Int64, K::AbstractArray{T}, pressures::Vector, freq=0.0) where T

    # Aloca o vetor de saída
    F = zeros(T,nnos)

    # O caso mais simples é quando pressures for vazio
    if isempty(pressures)
        return F
    end

    # Verifica a questão da frequência
    flag_freq = true
    if freq==0.0
        flag_freq = false
    end
 
    # Loop pelas entradas de pressures
    for p in pressures
 
        # Recupera a frequência
        f = p["freq"]
         
        # Testa o flag e a frequência 
        if !flag_freq || isapprox(f,freq)
            
           # Recupera o valor da pressão 
           press = p["value"] 
 
           # Recupera os nós 
           nodes = p["nodes"]
 
           # Loop pelos nós, calculando a contribuição 
           # da pressão no vetor de forças 
           for node in nodes
              
              F .-= K[:,node]*press

           end # node
 
        end #if 
 
    end # p 

    # Retorna o vetor de forças devido à pressão imposta
    return F

end

#
# Essential boundary conditions (pressures)
#
# Esta rotina simplemente copia os valores impostos 
# para as posições (nós) em P, modificando P no processo
#
# Como podemos aplicar pressões para diferentes 
# frequências, podemos indicar um valor específico 
# se for do interesse. Caso freq=0, consideramos 
# todas as frequências 
#
# Deve ser utilizada após a solução do sistema 
#
# K[livres,livres] \ P[livres]
#
# com 
#
# Mask_ebc!(P,pressures)
#
#
# TODO fase deveria entrar como a parte complexa ?
#
#
function Mask_ebc!(P::Vector{T},pressures,freq=0.0) where T

    # Se pressures for vazio, então saimos 
    if isempty(pressures) 
       return nothing
    end

    # Verifica a questão da frequência
    flag_freq = true
    if freq==0.0
        flag_freq = false
    end

    # Loop pelas entradas de pressures
    for p in pressures

        # Recupera a frequência
        f = p["freq"]
 
        # Testa o flag e a frequência 
        if !flag_freq || isapprox(f,freq)
           
           # Recupera o valor da pressão 
           press = p["value"] 

           # Recupera os nós 
           nodes = p["nodes"]

           # Aplica o valor 
           P[nodes] .= press

        end #if 

    end # p 

    # Não retornamos nada nesta rotina
    return nothing

end
