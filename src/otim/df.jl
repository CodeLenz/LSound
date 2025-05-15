function df(γ0::Vector,f::Function, δ=1E-6)

    # Iniciliza o vetor gradiente
    d = similar(γ0)

    # Para cada componente do vetor, perturba e calcula a diferença
    for i in eachindex(γ0)

        # Valor atual nesta posição
        γe = γ0[i]
        
        # Perturba para frente
        γ0[i] = γe + δ
        if γ0[i]>=1.0  
           error("variável não pode passar de 1.") 
        end
        ff    = f(γ0)

        # Perturba para trás
        γ0[i] = γe - δ
        if γ0[i]<=0  
            error("variável não pode ser menor do que zero") 
         end
         
        ft    = f(γ0)

        # Aproxima a derivada
        d[i] = (ff-ft)/(2*δ)

        # Desfaz a perturbação
        γ0[i] = γe

    end

    # Retorna a derivada
    return d

end