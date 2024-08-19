# The Bathe time integration method revisited for prescribing desired
# numerical dissipation
# Mohammad Mahdi Malakiyeh a , Saeed Shojaee a , Klaus-Jürgen Bathe b,⇑
#
# F(t) é uma função de t que 
# devolve um vetor 
#       M = Matriz de Massa 
#       K = Matriz de Rigidez
#       Δt = Intervalo de Tempo
#       TF = Tempo Final de simulação

function B1B2Bathe2d(M, C, K, F::Function,  livres, Δt, Tf;  U0=Float64[], V0=Float64[],
    γ  = 1/2, 
    δ  = 1/4,
    β1 = 0.39,
    β2 = 0.78)

    # Número de gls livres
    nl = length(livres) 

    # Dimensão do sistema (diferente do n da rotina principal)
    # pois aqui são todos os gls (livres e presos)
    n = size(M, 1)

    # Condições iniciais
    if isempty(U0)
        U0 = zeros(nl)
    else
        U0 = U0[livres] 
    end

    if isempty(V0)
        V0 = zeros(nl)
    else
        V0 = V0[livres]
    end

    # Constantes do método (aqui vou usar a nomenclatura do artigo
    # para não dar confusão nas equações)

    # Algumas constantes úteis
    c1 = γ*Δt
    c2 = c1*c1
    b1 = β2*(1-γ)*Δt
    b2 = b1*b1

    # Usadas para calcular o R2
    comum1 = Δt*(β2*(1-γ))^2
    d1 = (γ*(1-β1)+β2*(1-γ))/comum1
    d2 = (γ*β1 + (1-β2)*(1-γ))/comum1
    d3 = (γ*(1-β1)) / (β2*(1-γ))
    d4 = (γ*β1 + (1-β2)*(1-γ))/(β2*(1-γ))
    d5 = 1/(β2*(1-γ))

    # Para a aceleração final
    credo = c1*β1 + (1-γ)*Δt*(1-β2)

    
    # Copia só os gls livres
    KL = K[livres,livres] 
    CL = C[livres,livres] 
    ML = M[livres,livres] 

    # Calcula a aceleração no tempo 0
    A0 = zeros(nl)
    A0 .= ML \ (F(0.0)[livres] .- CL*V0  .- KL*U0)

    # Gera um lista de tempos discretos
    tempos = range(start=0.0, stop=Tf+Δt, step=Δt)

    # Matriz para guardar os deslocamentos em 
    # cada passo de tempo. Cada coluna é um 
    # passo de tempo
    MU = zeros(n, length(tempos))

    # Já guarda a condição inicial U0 aqui
    MU[livres,1] .= U0

    # Calcula a Matriz de coeficientes 1
    # Equação 6 do artigo
    K1 = (4/c2)*ML .+ (2/c1)*CL .+ KL

    # Pre-fatora a matriz de coeficientes. Somente 
    # posições livres
    LK1 = lu(K1)

    # Calcula a Matriz de coeficientes 2
    # Equação 16 do artigo
    K2 = (1/b2)*ML .+ (1/b1)*CL .+ KL

    # Pre-fatora a matriz de coeficientes. Somente 
    # posições livres
    LK2 = lu(K2)

    # Pré-aloca os vetores que serão utilizados na solução
    R1 = similar(U0)
    U1 = similar(U0)
    V1 = similar(U0)
    A1 = similar(U0)
    R2 = similar(U0)
    Rt = similar(U0)
    Ut = similar(U0)
    Vt = similar(V0)
    At = similar(A0)

    # Vamos ao loop principal
    coluna = 2 
    @showprogress "Calculando..." for i in 1:length(tempos)-1

        ############### PRIMEIRA ETAPA DO MÉTODO ########### 
        # Monta o vetor R1 de forças no tempo intermediário
        tγ = tempos[i] + c1

        # Eq. 7 do artigo
        R1 .= F(tγ)[livres] .+ ML*((4/c2)*U0 .+ (4/c1)*V0 .+ A0) .+ 
              CL*((2/c1)*U0 .+ V0)

        # Obtem os deslocamentos no tempo intermediário
        U1 .= LK1\R1

        # Agora vamos precisar calcular a velocidade e a 
        # aceleração no tempo intermediário. Essas equações
        # não estáo explícitas no artigo, mas podemos isolar
        # a velocidade usando a Eq. 2
        V1 .= (2/c1)*(U1.-U0) .- V0

        # Com Vc1, podemos calcular a aceleração no tempo
        # intermediário, usando a Eq. 1
        A1 .= (2/c1)*(V1 .- V0) .- A0

        ############### SEGUNDA ETAPA DO MÉTODO ########### 

        # Equação 17
        # usando as cts que definimos no começo da rotina
        R2 .= F(tempos[i+1])[livres] .+ 
            ML*((1/b2)*U0 .+ d1*V0 .+ d2*V1 .+ d3*A0 .+ d4*A1) .+ 
            CL*(d5*U0 .+ d3*V0 .+ d4*V1)   

        # Agora podemos obter Ut
        Ut .= LK2\R2

        # Usando a Eq. 13 podemos calcular as velocidades 
        Vt .= (1/b1)*(Ut.-U0) .- d3*V0 .- d4*V1

        # E, usando a Eq. 10, as acelerações
        At .= (1/b1)*(Vt.-V0 .-(c1*(1-β1))*A0 - credo*A1)

        # Guarda na coluna
        MU[livres, coluna] .= Ut

        # Incrementa a coluna
        coluna += 1

        # Adianta as medidas para o próximo tempo
        U0 .= Ut
        V0 .= Vt
        A0 .= At

    end

    return collect(tempos), MU

end
