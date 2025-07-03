
#=
Maria B. Dühring, Jakob S. Jensen, Ole Sigmund,
Acoustic design by topology optimization,
Journal of Sound and Vibration,
Volume 317, Issues 3–5,
2008,
Pages 557-575,
ISSN 0022-460X,
https://doi.org/10.1016/j.jsv.2008.03.042.
(https://www.sciencedirect.com/science/article/pii/S0022460X08002812)
=#


# Parametrização do inverso de ρ
#
# γe: variável de projeto
# ρ_ar: densidade do ar
# ρ2  : densidade da fase sólida
#
function fρ_duhring(γe;ρ_ar=1.204,ρ2=2643.0)     

    # Teste de consistência
    0<γe<1 || error("fρ:: γe inválido")

    # ρ do ar
    ρ1 = ρ_ar

    
    # Calcula a parametrização    [(ρ2/ρ1)^(-1) = (ρ1/ρ2)]
    return 1 + γe *((ρ1/ρ2) - 1)

end 

# Derivada da parametrização do inverso de ρ
function dfρ_duhring(γe;ρ_ar=1.204,ρ2=2643.0)

    # Teste de consistência
    0<γe<1 || error("dfρ:: γe inválido")

    # ρ do ar
    ρ1 = ρ_ar 
    
    # Calcula derivada da parametrização de ρ
    return ((ρ1/ρ2) - 1)

end 


# Parametrização do inverso de κ
function fκ_duhring(γe ; κ_ar=141.921E3,κ2=6.87E10)

    # Teste de consistência
    0<γe<1 || error("fκ:: γe inválido")

    # κ do ar
    κ1 = κ_ar

    # Calcula a parametrização    [(κ2/κ1)^(-1) = (κ1/κ2)]
    return 1 + γe *((κ1/κ2) - 1)

end 

# Derivada da parametrização do inverso de κ
function dfκ_duhring(γe; κ_ar=141.921E3,κ2=6.87E10)

    # Teste de consistência
    0<γe<1 || error("dfκ:: γe inválido")

     # κ do ar
    κ1 = κ_ar

    # Calcula a parametrização    [(κ2/κ1)^(-1) = (κ1/κ2)]
    return ((κ1/κ2) - 1)

end 


