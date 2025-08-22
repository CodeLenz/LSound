#
# Realiza uma sequência de análises harmônicas em uma lista de nf 
# frequências de excitação e guarda a solução em uma matriz nn × nf 
#
function Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities,pressures::Vector)

    # Calcula as matrizes globais
    K,M = Monta_KM_param(ne,coord,connect,γ,fρ,fκ)
    
    # E a de amortecimento
    # TODO adicionar amortecimento depois
    # C = Matriz_C(nn,damping,coord,connect)

    # Número de frequências
    nf = length(freqs)

    # Aloca o vetor de pressão (U)
    U = zeros(ComplexF64,nn)

    # Aloca matriz com os valores a serem monitorados
    MP = zeros(ComplexF64,nn,nf)

    # Aloca o vetor de forças 
    P = Array{ComplexF64}(undef,nn)

    # Loop pelas frequências
    contador = 1
    for f in freqs

        # Converte a freq de Hz para rad/s
        ω = 2*pi*f

        # Monta a matriz de rigidez dinâmica
        # Kd = K[livres,livres] .+ im*ω*C[livres,livres] .- (ω^2)*M[livres,livres]
        Kd = K  .- (ω^2)*M

        # Monta o vetor de forças, que depende da frequência  
        Vetor_P!(0.0,velocities,coord,connect,P,ω=ω)

        # Monta o vetor de forças devido às pressões impostas 
        # 
        # TODO VER A QUESTÃO DA FREQUÊNCIA
        #
        F_P =  P_pressure(nn, Kd, pressures)

        # A "força" total será a soma das duas parcelas
        P += F_P
            
        # Soluciona 
        U[livres] .= Kd[livres,livres]\P[livres]

        # Precisamos registrar os valores das pressões impostas em U
        # TODO VER A QUESTÃO DA FREQUÊNCIA
        Mask_ebc!(U,pressures)

        # Grava na matriz MP
        MP[:,contador] .= U
        
        # Incrementa o contador
        contador += 1

    end

    # Retorna MP, K e M
    return MP, K, M

end