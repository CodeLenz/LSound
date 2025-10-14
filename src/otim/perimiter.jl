#
# Perimiter of the design
#
function Perimiter(γ, neighedge, elements_design)

   # Initialize the permiter 
   P = 0.0

   # Loop over the design elements
   for ele in elements_design

       # Neighbor
       vizinhos = neighedge[ele]

       # Densidade do elemento 
       γe = γ[ele]

       # Loop over neighedge
       for viz in vizinhos
           P += 0.5*(γe - γ[viz])^2
       end #viz

   end #ele

   # return the perimiter 
   return P

end

#
# Delta de Kroenecker
#
function DK(a::Int,b::Int)
    ifelse(a==b,1.0,0.0)
end



#
# Gradiente do perímetro
#
function dPerimiter(ne, γ, neighedge, elements_design)

   # Vetor de derivadas
   dP = zeros(ne)

   # Loop over the positions of dP
   for m in elements_design

        # Loop over the design elements
        for ele in elements_design

            # Neighbor
            vizinhos = neighedge[ele]

            # Densidade do elemento 
            γe = γ[ele]

            # Pré-computa DK(ele,m)
            p1 = DK(ele,m)
            
            # Loop over neighedge
            for viz in vizinhos

                # Pej
                Pej = (γe - γ[viz])
                
                # Derivada
                dP[m] += Pej*(p1 - DK(viz,m))
                

            end #viz

        end #ele

    end # m

   # return the sensitivity
   return dP

end