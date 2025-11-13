

#
# n => número de variáveis de projeto
#
function LP(c, A, b, γ)

   #
   # Cria o modelo vazio e associa a um otimizador
   # O Alpine utiliza mais de um pacote de otimização, dependendo 
   # da etapa que ele está realizando. O mip_solver é bem crítico para 
   # a primeira etapa, em que é contínua. O Gurobi é uma boa opção, mas 
   # precisamos instalar a licença no computador. Uma outra opção é utilizar
   # o Ipopt ou  HigHS.
   # O Cbc é o otimizador para a etapa discreta (branch and bound), que é 
   # realizada depois da etapa contínua.
   ipopt  = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
   gurobi = optimizer_with_attributes(Gurobi.Optimizer, "output_flag" => false)
   highs  = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
   cbc    = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)
    
   # minlp_solver é o solver para o problema binário (0/1)
   # nlp_solver é o solver para o problema contínuo do LP
   # mip_solver é o utilizado na primeira etapa 
   model = Model(
      optimizer_with_attributes(
         Alpine.Optimizer,
         "minlp_solver" => highs, #<- para binario. Aqui também dá para usar o cbc
         #"nlp_solver" => ipopt,  #<- para contínuo
         "mip_solver" => gurobi,
         #"mip_solver" => highs,
      ),
   )

   # Vamos descobrir o número de variáveis de projeto e de restrições 
   nc = length(c)
   mb = length(b)
   ma,na = size(A)

   # E testar a consistência dos dados
   nc==na || error("LP: dimensões inconsistentes")
   mb==ma || error("LP: dimensões inconsistentes")

   # Cria o vetor com as restrições laterais para cada variável
   Δxi = zeros(Int,na)
   Δxs = zeros(Int,na)

   for i in LinearIndices(Δxi)
       Δxi[i] = -round(Int,γ[i])
       Δxs[i] =  round(Int,1-γ[i])
   end

   # Cria um vetor de variáveis de projeto
   @variable(model, Δxi[i] <= Δx[i=1:na] <= Δxs[i], Int)

   # Monta todas as restrições ao mesmo tempo 
   @constraint(model, A * Δx .<= b)

   # Monta a função objetivo
   @objective(model, Min, c' * Δx)

   # Resolve o problema 
   #optimize!(model)
   redirect_stdout((()->optimize!(model)),open("nul", "w"))

   # Valor do objetivo
   # objective_value(model)

   # Vetor de variáveis de projeto no ponto de ótimo
   Δxopt = value(Δx)

   # Vamos testar as restrições
   @show [A*Δxopt b]
   @show termination_status(model)

   # @show c
   # @show A
   # @show xopt

   # Retorna o vetor de variáveis de projeto no ponto de ótimo
   return Δxopt

end


