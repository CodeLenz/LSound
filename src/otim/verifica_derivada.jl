#
# Validação das derivadas por DFC
#
function Verifica_derivada(γ,nn,ne,coord,connect,fρ,fκ,freqs,livres,velocities,pressures,nodes_target,elements_design)
    
    # Vamos validar a derivada usando diferenças finitas
    function f_(γ,nn,ne,coord,connect,fρ,fκ,freqs,livres,velocities,pressures,nodes_target)

        MP,_ =  Sweep(nn,ne,coord,connect,γ,fρ,fκ,freqs,livres,velocities,pressures) 

        # Calcula a função objetivo SPL_w
        objetivo = Objetivo(MP,nodes_target)

        return objetivo

    end
    f(γ) = f_(γ,nn,ne,coord,connect,fρ,fκ,freqs,livres,velocities,pressures,nodes_target)

    # Calcula a derivada por DFC
    d_numerica = df(γ,f,elements_design,1E-6)

end