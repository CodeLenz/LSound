
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
# ρ_ar : densidade do ar
# ρ_s  : densidade da fase sólida
#
function fρ(γe;ρ_ar=1.204,ρ_s=2643.0)     

    # Teste de consistência
    0<=γe<=1 || error("fρ:: γe inválido")

    # Parametrização 
    return (1/ρ_ar) + γe*(1/ρ_s-1/ρ_ar)


end 

# Derivada da parametrização do inverso de ρ
function dfρ(γe;ρ_ar=1.204,ρ_s=2643.0)

    # Teste de consistência
    0<=γe<=1 || error("dfρ:: γe inválido")
    
    # Calcula derivada da parametrização de ρ
    return (1/ρ_s-1/ρ_ar)

end 


# Parametrização do inverso de κ
function fκ(γe ; κ_ar=141.921E3,κ_s=6.87E10)

    # Teste de consistência
    0<=γe<=1 || error("fκ:: γe inválido")
    
    # Calcula a parametrização 
    return (1/κ_ar) + γe*(1/κ_s-1/κ_ar)

end 

# Derivada da parametrização do inverso de κ
function dfκ(γe; κ_ar=141.921E3,κ_s=6.87E10)

    # Teste de consistência
    0<=γe<=1 || error("dfκ:: γe inválido")

    # Calcula a parametrização 
    return (1/κ_s-1/κ_ar)

end 


