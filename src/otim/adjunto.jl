#
# Monta o vetor adjunto para uma frequência específica
#
function F_adj(nodes_target::Vector{T1},P::Vector{T2}) where {T1,T2}

    # Inicializa o vetor de carregamento adjunto
    F = similar(P) 
    
    # Zera o vetor
    fill!(F,zero(T2))

    # Somatório nas posições de nodes_target
    for p in nodes_target

        F[p] += 2*conj(P[p])

    end

    # Retorna F adjunto, sem os termos constantes que multiplicam
    return F 

end

#
# Calcula a derivada da matriz de rigidez dinâmica
#
function Derivada_KM(et,γe,dfρ::Function,dfκ::Function,X::Array)

    # Calcular as derivadas da inversa de ρ e da inversa
    # de κ em relação à γe
    diρ = dfρ(γe)
    diκ = dfκ(γe)

    # Monta as matrizes dos elementos
    # usando as derivadas das propriedades em 
    # relação à γe
     if et==3
        Ke, Me = KMe_bi4(diρ,diκ,X)
     elseif et==2
        Ke, Me = KMe_tri3(diρ,diκ,X)
     elseif et==4
         Ke, Me = KMe_tet4(diρ,diκ,X)   
     elseif et==5
        Ke, Me = KMe_hex8(diρ,diκ,X)
     elseif et==7
         Ke, Me = KMe_pyr5(diρ,diκ,X) 
      else
         error("Derivada_KM::Elemento não definido")
     end

     # Devolve as derivadas
     return Ke, Me

end

#
# Calcula a derivada da função objetivo
# Média simples do SPL em cada frequência
#
function Derivada(ne,nn,γ::Vector{T0},connect::Matrix{T1},coord::Matrix,
                  K::AbstractMatrix{T0},M::AbstractMatrix{T0},
                  livres::Vector{T1},freqs::Vector{T0}, dfρ::Function, dfκ::Function,
                  nodes_target::Vector{T1},target::Array{T2},p0=20E-6) where {T0,T1,T2}

    # Define o vetor de derivadas
    d = zeros(ne)

    # Como estamos assumindo que todos os nós target tem o mesmo p0
    P02 = p0^2 

    # Calcula a constante -10/(ln(10)*nt*P02)
    cte = -10/(log(10)*length(nodes_target))

    # Número de frequências
    Nf = length(freqs)

    # Define λ fora do loop, para reaproveitar
    λ = zeros(T2,nn)

    # Loop pelas frequências
    coluna = 1
    for f in freqs
        
        # Converte a freq para rad/s
        ω = 2*pi*f

        # Recupera as pressões para essa frequência (coluna de target)
        P = target[:,coluna]

        # Monta a matriz de rigidez dinâmica
        Kd = K[livres,livres]  .- (ω^2)*M[livres,livres]

        # Monta o vetor adjunto para essa frequência
        F = F_adj(nodes_target,P)

        # Escalona F pela cte e pelo Nf
        F .= F*cte/Nf
        
        # Soluciona o problema adjunto, obtendo λ^n
        λ[livres] .= Kd[livres,livres]\F[livres]

        # Loop pelos elementos
        for ele = 1:ne

            # Tipo de elemento
            etype = connect[ele,1]

            # Localizações 
            nos, X = Nos_Coordenadas(ele,etype,coord,connect)
  
            # Pressão nos nós do elementos
            pe = P[nos]

            # Vetor adjunto nos nós do elemento
            λe = λ[nos]

            # Variável de projeto do elemento
            γe = γ[ele]

            # Calcula a derivada da rigidez dinâmica do elemento
            dKe, dMe = Derivada_KM(etype,γe,dfρ,dfκ,X)

            # Derivada da matriz dinâmica do elemento
            dKde = dKe - dMe*ω^2  

            # Calcula a derivada e sobrepõe na posição do elemento
            d[ele] += 2*real(λe'*dKde*pe)

        end

        # Atualiza a coluna em target
        coluna += 1

    end

    # Retorna a derivada
    return d
         
end