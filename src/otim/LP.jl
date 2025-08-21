

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
   # o Ipopt u o HigHS.
   # O Cbc é o otimizador para a etapa discreta (branch and bound), que é 
   # realizada depois da etapa contínua.
   ipopt  = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
   #gurobi = optimizer_with_attributes(Gurobi.Optimizer, "output_flag" => false)
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
         #"mip_solver" => gurobi,
         "mip_solver" => highs,
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
   xi = zeros(Int,na)
   xs = zeros(Int,na)

   for i in LinearIndices(xi)
       xi[i] = -round(Int,γ[i])
       xs[i] =  round(Int,1-γ[i])
   end

   # Cria um vetor de variáveis de projeto
   @variable(model, xi[i] <= x[i=1:na] <= xs[i], Int)

   # Monta todas as restrições ao mesmo tempo 
   @constraint(model, A * x .<= b)

   # Monta a função objetivo
   @objective(model, Min, c' * x)

   # Resolve o problema 
   optimize!(model)

   # Valor do objetivo
   objective_value(model)

   # Vetor de variáveis de projeto no ponto de ótimo
   xopt = value(x)

   # Vamos testar as restrições
   @show [A*xopt b]

   # @show c
   # @show A
   # @show xopt

   # Retorna o vetor de variáveis de projeto no ponto de ótimo
   return xopt

end


