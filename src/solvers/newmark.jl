
#
# F(t) é uma função de t que 
# devolve um vetor 
#       M = Matriz de Massa 
#       K = Matriz de Rigidez
#       Δt = Intervalo de Tempo
#       TF = Tempo Final de simulação

function Newmark(M, C, K, F::Function, livres, Δt, Tf;  U0=Float64[], V0=Float64[], γ = 1/2, β = 1/4)

    # Vamos trabalhar somente com os gls livres
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

    # Constantes do método 
    

    # Alias para as coordenadas livres
    ML =  M[livres,livres]
    KL =  K[livres,livres]
    CL =  C[livres,livres]

    # Cache para o vetor de forças
    forca =  F(0.0)[livres]

    # Calcula a aceleração no tempo 0
    A0 = zeros(nl)
    A0 .= ML \ (forca .- CL*V0  .- KL*U0)

    # Gera um lista de tempos discretos
    tempos = range(start=0.0, stop=Tf+Δt, step=Δt)

    # Matriz para guardar os deslocamentos em 
    # cada passo de tempo. Cada coluna é um 
    # passo de tempo
    MU = zeros(n, length(tempos))

    # Já guarda a condição inicial U0 aqui
    MU[livres,1] .= U0

    # Calcula a Matriz de coeficientes
    MA = ML .+ (KL .* β .* Δt^2) .+ (CL .* γ .* Δt)

    # Inicializa o solver linear
    prob = LinearProblem(MA, U0)
    linsolve = init(prob)
    
    # Pre-aloca Ut, Vt, A, U e V
    Ut = similar(U0)
    U  = similar(U0)
    Vt = similar(V0)
    V  = similar(V0)
    A  = similar(A0)

    # Vamos ao loop principal
    coluna = 2 
    @showprogress "Newmark... " for i in 1:length(tempos)-1

        # Calcula o preditor de deslocamento (Eq. a)
        @. Ut = U0 + Δt*V0 + (Δt^2/2) * (1 - 2*β)*A0
        
        # Calcula o preditor da velocidade (Eq. b)
        @. Vt = V0 + (1 - γ)*Δt*A0

        # Força equivalente 
        linsolve.b =  (F(tempos[i+1])[livres] .- KL*Ut .- CL*Vt)

        #Calcula a aceleração para frente
        sol = solve!(linsolve)
        A .= sol.u
        
        # Calcula o deslocamento em t+Δt
        @. U = Ut + A*β*Δt^2
         
        # Guarda na coluna
        MU[livres, coluna] .= U

        # Incrementa a coluna
        coluna += 1

        # Calcula a velocidade em t+Δt
        @. V = Vt + A*γ* Δt
        
        # Adianta as medidas para o próximo tempo
        U0 .= U
        V0 .= V
        A0 .= A

    end

    return collect(tempos), MU

end
