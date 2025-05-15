#
# Garante que os valores de γ nas posições fixas tenham os 
# valores prescritos
#
function Fix_γ!(γ::Vector,elements_fixed::Vector, values_fixed::Vector)

    γ[elements_fixed] .= values_fixed

end

#
# Zeras as posições de derivada com elementos fixos
#
function Fix_D!(D::Vector{T},elements_fixed::Vector) where T

    D[elements_fixed] .= zero(T)

end
