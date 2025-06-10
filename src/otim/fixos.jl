# ===================================================================================
# γ                 -> vetor de variáveis de projeto
#                      vector of design variables
# 
# elements_fixed    -> vetor de elementos que não são de projeto (estão 'fixos')
#                      vector of non-design elements (they are 'fixed')
#
# values_fixed      -> vetor com os valores para os elementos fixos
#                      vector with values ​​for the fixed elements
# ===================================================================================
#
# Garante que os valores de γ nas posições fixas tenham os 
# valores prescritos
#
function Fix_γ!(γ::Vector,elements_fixed::Vector, values_fixed::Vector)

    γ[elements_fixed] .= values_fixed

end

#
# Zera as posições de derivada com elementos fixos
#
function Fix_D!(D::Vector{T},elements_fixed::Vector) where T

    D[elements_fixed] .= zero(T)

end
