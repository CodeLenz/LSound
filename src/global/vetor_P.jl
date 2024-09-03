 # Rotina que monta o vetor de "forças" para um determinando tempo t
 function Vetor_P!(t,nn,materials,velocities,coord,connect,P)

    # Zera todo o vetor 
    fill!(P,0.0)

    # Loop pelo vetor velocities. Cada linha deste vetor é um dicionário
    for dvn in velocities

      # Recover data from Dictionary
      valor = dvn["value"]
      freq  = 2*pi*dvn["freq"]  
      fase  = 2*pi*dvn["phase"]
      elements   = dvn["elements"] 

      # Derivative of vn w.r.t time
      qn = -valor*freq*cos(freq*t + fase)

      # Agora precisamos fazer um loop sobre os elementos
      # e suas arestas
      for i in axes(elements,1)

        # Element and edge
        ele  = elements[i,1]
        edge = elements[i,2]

        # Element type
        et = connect[ele,1]

        # Material for this element
        mat = connect[ele,2]

        # Find material density
        ρ = materials[mat,1]

        # Find nodes and coordinates
        nos,X = Nos_Coordenadas(ele,et,coord,connect)

        # Local vector 
        if et==3
            Pn = Edge_load_local_bi4(edge,-ρ*qn,X)
        elseif et==2
            Pn = Edge_load_local_tri3(edge,-ρ*qn,X)
        else
          error("Elemento não definido")
        end

        # Add to the global vector
        P[nos] .+=  Pn

      end #i
 
    end # dict

    # Mesmo sendo uma função !, vamos retornar
    # o vetor, pois depois vamos fazer um f(t)
    # com essa rotina
    return P

 end