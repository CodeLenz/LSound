
#
# Programar depois para fazer a validação das derivadas por DFC
#
function Verifica_derivada(γ,nn,ne,coord,connect,fρ,fκ,freqs,livres,velocities)
    
    # Vamos validar a derivada usando diferenças finitas
    function f_(γ,nn,ne,coord,connect,fρ,fκ,freqs,livres,velocities)

        MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities) 

        # Calcula a função objetivo SPL_w
        objetivo = Objetivo(MP,nodes_target)

        return objetivo

    end
    f(γ) = f_(γ,nn,ne,coord,connect,fρ,fκ,freqs,livres,velocities)

    # Calcula a derivada por DFC
    println("Entrando em numérica")
    d_numerica = df(γ,f,1E-8)

end

