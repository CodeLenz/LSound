#
# Monta o vetor adjunto para uma frequência específica
#
# nodes_target -> vector of Target nodes   (vetor com os 'nós monitorados' )
# O vetor P foi obtido via Sweep.
#
function F_adj(nodes_target::Vector{T1},P::Vector{T2}) where {T1,T2}

    # Inicializa o vetor de carregamento adjunto
    F = similar(P) 
    
    # Zera o vetor
    fill!(F,zero(T2))

    # Somatório nas posições de nodes_target
    for p in nodes_target

        F[p] += -conj(P[p])

    end

    # Retorna F adjunto, sem os termos constantes que multiplicam
    return F 

end

# ===================================================================================
# Calcula a derivada da matriz de rigidez dinâmica
#
# et  ->  tipo de elemento
# γe  ->  variável de projeto do elemento
# dfρ ->  função que parametriza a derivada da densidade 
# dfκ ->  função que parametriza a derivada do módulo de compressibilidade
# X   ->  matriz com as coordenadas do elemento 
#
function Derivada_KM(et,γe,dfρ::Function,dfκ::Function,X::Array)

    # Calcular as derivadas da inversa de ρ e da inversa
    # de κ em relação à γe
    diρ = dfρ(γe)
    diκ = dfκ(γe)

    # Monta as matrizes dos elementos
    # usando as derivadas das propriedades 
    # em relação à γe
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

# ===================================================================================
# Calcula a derivada da função objetivo
#
# Média simples do SPL em cada frequência
#
function Derivada(ne,nn,γ::Vector{T0},connect::Matrix{T1},coord::Matrix{T0},
                  K::AbstractMatrix{T0},M::AbstractMatrix{T0},
                  livres::Vector{T1},freqs::Vector{T0},
                  pressures::Vector, dfρ::Function, dfκ::Function,
                  nodes_target::Vector{T1},MP::Matrix{T2},
                  elements_design::Vector,A::Vector,p0=20E-6) where {T0,T1,T2}


    #
    # TESTE 
    #       
    #γ[γ.==0]    .= 1E-3   
    
    # Define o vetor de derivadas
    d = zeros(ne)

    # Número de frequências
    Nf = length(freqs)

    # Número de nós para monitorar a pressão
    nt = length(nodes_target)

    # Calcula a constante 10/(ln(10)*nt)
    cte = 10/(log(10)*nt)
 
    # Define λ fora do loop, para reaproveitar
    λn = zeros(T2,nn)

    # Aloca antes do loop
    P = MP[:,1]

    # Aloca Fn 
    Fn = similar(P)

    # Aloca antes do loop
    Kd = similar(K[livres,livres])

    # Loop pelas frequências
    coluna = 1
    for f in freqs
        
        # Converte a freq para rad/s
        ωn = 2*pi*f

        # Sensibilidade para essa frequência 
        An = A[coluna]

        # Recupera as pressões para essa frequência (coluna de target)
        P .= MP[:,coluna]

        # Monta a matriz de rigidez dinâmica
        Kd .= K[livres,livres]  .- (ωn^2)*M[livres,livres]

        # Monta o vetor adjunto para essa frequência
        Fn .= F_adj(nodes_target,P)

        # Calcula o Pn2
        P2 = sum((abs.(P[nodes_target])).^2)

        # Média (pelo número de pontos em nodes_target)
        P2avg = P2 / nt

        # Escalona F pela cte e pelo Nf
        Fn .= An*Fn*cte/(P2avg)
        
        # Soluciona o problema adjunto, obtendo λ^n
        λn[livres] .= Kd\Fn[livres]

        # Loop pelos elementos
        for ele in elements_design

            # Tipo de elemento
            etype = connect[ele,1]

            # Localizações 
            nos, X = Nos_Coordenadas(ele,etype,coord,connect)
  
            # Pressão nos nós do elementos
            pe = P[nos]

            # Vetor adjunto nos nós do elemento
            λe = λn[nos]

            # Variável de projeto do elemento
            γe = γ[ele]

            # Calcula a derivada da rigidez dinâmica do elemento
            dKe, dMe = Derivada_KM(etype,γe,dfρ,dfκ,X)

            # Derivada da matriz dinâmica do elemento
            dKde = dKe - dMe*ωn^2  

            # Derivada de Fn [dFn/dγm]
            # dFp =  Derivada_forca_pressao(nos,pressures,dKde)
            # -2*real(transpose(λe)*dFp)

            # Calcula a derivada e sobrepõe na posição do elemento
            d[ele] += 2*real(transpose(λe)*dKde*pe) 

        end # Elemento

        # Atualiza a coluna em target
        coluna += 1

    end # Frequência

    #
    # TESTE
    #
    # γ[γ.==0]    .= 0
    


    # Retorna a derivada, lembrando de dividir pelo número de 
    # frequências
    return d./Nf
         
end

# ===================================================================================
# Derivada de Fn [dFn/dγm]
#
# "Forças" devido a pressões impostas para um elemento com 
#  vetor de nós = nos
#
# dKde já é a derivada da matriz Kd do elemento em relação a sua 
# variável de projeto γ_m
#
function Derivada_forca_pressao(nos::Vector,pressures::Vector,dKde::AbstractMatrix)

     # Aloca o vetor de saída
     # com a dimensão do número de nós do elemento
     dF = zeros(ComplexF64,length(nos))

     # Loop pelas entradas de pressures
     for p in pressures
       
        # Recupera o valor da pressão 
        png = p["value"] 
  
        # Recupera os nós 
        nodes_p = p["nodes"]
  
        # Loop pelos nós, calculando a contribuição 
        # da pressão no vetor de forças 
        for node in nodes_p

            # Precisamos verificar se o nó no qual a pressão está sendo aplicada
            # existe no elemento e, caso verdadeiro, utilizar a coluna da matriz LOCAL
            # correspondente
            pos = findfirst(x->x==node,nos)

            # Se existir uma posição correspondente, podemos
            # calcular a derivada
            if !isnothing(pos)
                dF .-= dKde[:,pos]*png
            end 
 
        end # node
  
    end # p 
 
    # Retorna o vetor da derivada das forças devido à pressão imposta
    return dF

end