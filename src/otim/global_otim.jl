

# Parametrização do inverso de ρ
function fρ(γe,ψ=2.0)

    # ρ do ar
    ρa = 1.21
    ρr = ρa*1E7

    # Calcula a parametrização
    return (1/ρa) + (γe^ψ)*(1/ρr - 1/ρa)

end 

# Parametrização do inverso de κ
function fκ(γe,ψ=2.0)

    # κ do ar
    κa = 1.42E5
    κr = κa*1E9

    # Calcula a parametrização
    return (1/κa) + (γe^ψ)*(1/κr - 1/κa)

end 


#
# Monta_KM, versão para otimização
#
# γ é um vetor de variáveis de projeto [0,1]
# fρ é uma função que parametriza ρ
# fκ é uma função que parametriza κ
function Monta_KM2(ne,coord,connect,γ::Vector,fρ::Function,fκ::Function)  
    
    # Aloca vetores para a montagem eficiente 
    # das matrizes esparsas
    I = Int64[]
    J = Int64[]
    VK = Float64[]
    VM = Float64[]
 
    # Loop pelos elementos
    for ele=1:ne

        # Variável de projeto do elemento
        γe = γ[ele]

        # Calcula as inversas 
        iρ = fρ(γe)
        iκ = fκ(γe)

        # Tipo de elemento
        et = connect[ele,1]

        # Descobre nos, X e Y para este elemento
        nos, X = Nos_Coordenadas(ele,et,coord,connect) 

        # Monta as matrizes dos elementos
        if et==3
           Ke, Me = KMe_bi4(iρ,iκ,X)
        elseif et==2
           Ke, Me = KMe_tri3(iρ,iκ,X)
        elseif et==4
            Ke, Me = KMe_tet4(iρ,iκ,X)    
        elseif et==5
           Ke, Me = KMe_hex8(iρ,iκ,X) 
        elseif et==7
            Ke, Me = KMe_pyr5(iρ,iκ,X)    
         else
            error("Elemento não definido")
        end
 
        # Sobreposição das locais nas globais
        for i in LinearIndices(nos)
            ni = nos[i]
            for j in LinearIndices(nos)
                push!(I,ni)
                push!(J,nos[j])
                push!(VK,Ke[i,j])
                push!(VM,Me[i,j])
            end
        end

    end

    # Retorna as matrizes globais
    return sparse(I,J,VK), sparse(I,J,VM)

end

