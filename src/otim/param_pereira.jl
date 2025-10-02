
#=
Topology optimization of acoustic systems with a multiconstrained BESO approach,
Finite Elements in Analysis and Design,
Volume 201,
2022,
103701,
ISSN 0168-874X,
https://doi.org/10.1016/j.finel.2021.103701.
(https://www.sciencedirect.com/science/article/pii/S0168874X21001736)
=#

# Parametrização do inverso de ρ
function fρ(γe,ψ=2.0;ρ_ar=1.21,ρ2=1.21*1E7)

    # Teste de consistência
    0<=γe<=1 || error("fρ:: γe inválido")

    # ρ do ar
    ρa = ρ_ar

    # Fase sólida
    ρr = ρ2
    
    # Calcula a parametrização
    return (1/ρa) + (γe^ψ)*(1/ρr - 1/ρa)

end 

# Derivada da parametrização do inverso de ρ
function dfρ(γe,ψ=2.0;ρ_ar=1.21,ρ2=1.21*1E7)

    # Teste de consistência
    0<=γe<=1 || error("dfρ:: γe inválido")

    # ρ do ar
    ρa = ρ_ar

    # Fase sólida
    ρr = ρ2
    
    # Calcula a parametrização
    return ψ*(γe^(ψ-1))*(1/ρr - 1/ρa)

end 


# Parametrização do inverso de κ
function fκ(γe,ψ=2.0;κ_ar=1.42E5,κ2=1.42E5*1E9)

    # Teste de consistência
    0<=γe<=1 || error("fκ:: γe inválido")

    # κ do ar
    κa = κ_ar

    # fase sólida 
    κr = κ2

    # Calcula a parametrização
    return (1/κa) + (γe^ψ)*(1/κr - 1/κa)

end 

# Derivada da parametrização do inverso de κ
function dfκ(γe,ψ=2.0;κ_ar=1.42E5,κ2=1.42E5*1E9)

    # Teste de consistência
    0<=γe<=1 || error("dfκ:: γe inválido")

    # κ do ar
    κa = κ_ar

    # fase sólida 
    κr = κ2

    # Calcula a parametrização
    return ψ*(γe^(ψ-1))*(1/κr - 1/κa)

end 

