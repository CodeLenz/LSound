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

            P += sqrt((γe - γ[viz])^2)

       end #viz

   end #ele

   # return the perimiter 
   return P

end

#
# Delta de kroenecker
#
function DK(a::Int,b::Int)
    ifelse(a==b,1.0,0.0)
end



#
# Gradiente do perímetro
#
function dPerimiter(ne, γ, neighedge, elements_design, δ=1E-6)

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

            # Loop over neighedge
            for viz in vizinhos

                # Pej
                Pej = sqrt((γe - γ[viz])^2 + δ^2)
                
                # Derivada
                dP[m] += (1/Pej)*(γe - γ[viz])*(DK(ele,m) - DK(viz,m))

            end #viz

        end #ele

    end # m

   # return the sensitivity
   return dP

end