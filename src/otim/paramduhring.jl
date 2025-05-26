#
# Parametrização de acordo ao artigo da Maria Duhring   (tentativa)
#
# "2.2. Design variables and material interpolation"
#
# γe = ϵ   <-- para conferir com a nomenclatura do artigo⋆
#

# Parametrização do inverso de ρ
function fρ(γe)     

    # Teste de consistência
    0<γe<=1 || error("fρ:: γe inválido")

    # ρ do ar
    ρ1 = 1.204 

    # ρ do material 
    # no caso alumínio    < -- poderia colocar como entrada da função? 
    ρ2 = 2643.0

    # Calcula a parametrização    [(ρ2/ρ1)^(-1) = (ρ1/ρ2)]
    return 1 + γe *((ρ1/ρ2) - 1)

end 

# Derivada da parametrização do inverso de ρ
function dfρ(γe)

    # Teste de consistência
    0<γe<=1 || error("dfρ:: γe inválido")

    # ρ do ar
    ρ1 = 1.204 

    # ρ do material 
    # no caso alumínio
    ρ2 = 2643.0

    # Calcula derivada da parametrização de ρ
    return ((ρ1/ρ2) - 1)

end 


# Parametrização do inverso de κ
function fκ(γe)

    # Teste de consistência
    0<γe<=1 || error("fκ:: γe inválido")

    # κ do ar
    κ1 = 141.921E3

    # κ do material 
    # no caso alumínio
    κ2 = 6.87E10

    # Calcula a parametrização    [(κ2/κ1)^(-1) = (κ1/κ2)]
    return 1 + γe *((κ1/κ2) - 1)

end 

# Derivada da parametrização do inverso de κ
function dfκ(γe)

    # Teste de consistência
    0<γe<=1 || error("dfκ:: γe inválido")

     # κ do ar
    κ1 = 141.921E3

    # κ do material 
    # no caso alumínio
    κ2 = 6.87E10

    # Calcula a parametrização    [(κ2/κ1)^(-1) = (κ1/κ2)]
    return ((κ1/κ2) - 1)

end 


