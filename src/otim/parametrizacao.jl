
# Parametrização do inverso de ρ
function fρ(γe,ψ=2.0)

    # Teste de consistência
    0<γe<=1 || error("fρ:: γe inválido")

    # ρ do ar
    ρa = 1.21
    ρr = ρa*1E7

    # Calcula a parametrização
    return (1/ρa) + (γe^ψ)*(1/ρr - 1/ρa)

end 

# Derivada da parametrização do inverso de ρ
function dfρ(γe,ψ=2.0)

    # Teste de consistência
    0<γe<=1 || error("dfρ:: γe inválido")

    # ρ do ar
    ρa = 1.21
    ρr = ρa*1E7

    # Calcula a parametrização
    return ψ*(γe^(ψ-1))*(1/ρr - 1/ρa)

end 


# Parametrização do inverso de κ
function fκ(γe,ψ=2.0)

    # Teste de consistência
    0<γe<=1 || error("fκ:: γe inválido")

    # κ do ar
    κa = 1.42E5
    κr = κa*1E9

    # Calcula a parametrização
    return (1/κa) + (γe^ψ)*(1/κr - 1/κa)

end 

# Derivada da parametrização do inverso de κ
function dfκ(γe,ψ=2.0)

    # Teste de consistência
    0<γe<=1 || error("dfκ:: γe inválido")

    # κ do ar
    κa = 1.42E5
    κr = κa*1E9

    # Calcula a parametrização
    return ψ*(γe^(ψ-1))*(1/κr - 1/κa)

end 

