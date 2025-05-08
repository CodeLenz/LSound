

#
# Realiza uma sequência de análises harmônicas em uma lista de nf 
# frequências de excitação e guarda a solução em uma matriz nn × nf 
#
function Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities)

    # Calcula as matrizes globais
    K,M = Monta_KM2(ne,coord,connect,γ,fρ,fκ)
    
    # E a de amortecimento
    # TODO adicionar amortecimento depois
    # C = Matriz_C(nn,damping,coord,connect)

    # Número de frequências
    nf = length(freqs)

    # Aloca matriz com os valores a serem monitorados
    MP = zeros(ComplexF64,nn,nf)

    # Aloca o vetor de forças 
    P = Array{Float64}(undef,nn)

    # Loop pelas frequências
    contador = 1
    for f in freqs

        # Converte a freq de Hz para rad/s
        ω = 2*pi*f

        # Monta a matriz de rigidez dinâmica
        # Kd = K[livres,livres] .+ im*ω*C[livres,livres] .- (ω^2)*M[livres,livres]
        Kd = K[livres,livres]  .- (ω^2)*M[livres,livres]

        # Monta o vetor de forças, que depende da frequência  
        Vetor_P!(0.0,velocities,coord,connect,P,ω=ω)

        # Soluciona apenas para os gls livres do problema
        MP[livres,contador] .= Kd\P[livres]
  
        # Incrementa o contador
        contador += 1

    end

    # Retorna MP, K e M
    return MP, K, M

end