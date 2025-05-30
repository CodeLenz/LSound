#
# Read a .msh file and process data to a specific format
#
# Daniele (Acoustics)
#

#
# Physical groups que esta rotina consegue processar
#
#
# Material,nome,id,dens,c,Z [ surfaces (and volumes) ]
#
# Open [ lines and/or nodes and/or surfaces]
# 
# Vn,value,freq,phase (in degrees) [ lines and/or nodes and/or surfaces]
#
# Pressure,value,freq,phase (in degrees) [ lines and/or nodes and/or surfaces]
#
# Yn, value [ lines and/or nodes and/or surfaces]
#
# Probe [ lines and/or nodes and/or surfaces]
#
#
#
# For optimization, there are some specific phisical groups
#
# Target [ lines and/or nodes and/or surfaces]
# 
# Fixed, value [surfaces or volumes]
#
#
# Elementos que estão implementados
#
#          2D
# 2 -> triangular (linear)
# 3 -> quadrangular (linear)
#
#          3D
# 4 -> Tetrahedra (linear)
# 5 -> hexaedra (linear)
# 7 -> pyramid (linear)
#
#
function Parsemsh_Daniele(meshfile::String,verbose=false)
 
    # Primeiro precisamos definir se a malha é 2D ou 3D
    elist = Lgmsh_import_etypes(meshfile)

    # Se tivermos elementos do 4/5/7, então é 3D. Do contrário,
    # é 2D. Observe que ter 2/3 não é uma indicação direta de 
    # que a malha é 2D, pois o gmsh também gera esses elementos
    # para malhas 3D.
    dimensao = 2
    et = [2,3]
    if (4 in elist) || (5 in elist) || (7 in elist)
        dimensao = 3
        et = [4,5,7]
    end

    if verbose
        println("Solucionando um problema de dimensão $dimensao")
    end

    # Maximum number of nodes in the elements of the mesh
    nmax = maximum(Lgmsh_nodemap()[et])

    # Read mesh
    nn, coord, ne, etypes, connect, etags = Readmesh(meshfile,et)

    # Le todos os grupos físicos do arquivo 
    pgroups, pgnames = Lgmsh_import_physical_groups(meshfile)

    # Vector with Dicts of materials
    materials = Dict{String,Union{Float64,Int64,Vector{Int64}}}[]

    # Local dict inside the loop
    localD_m = Dict{String,Union{Float64,Int64,Vector{Int64}}}()

    # Vector with Dicts of normal velocities on nodes
    velocities = Dict{String,Union{Float64,Matrix{Int64}}}[]

    # Vector with Dicts of pressure on nodes
    pressures = Dict{String,Union{Float64,Vector{Int64}}}[]

    # Local dict inside the loop
    localD_vn = Dict{String,Union{Float64,Matrix{Int64}}}()

    # Local dict inside the loop
    localD_pressure = Dict{String,Union{Float64,Vector{Int64}}}()

    # Vector with Dicts of damping in faces
    damping = Dict{String,Union{Float64,Matrix{Int64}}}[]

    # Local dict inside the loop
    localD_damp = Dict{String,Union{Float64,Matrix{Int64}}}()

    # Vector of OPEN nodes
    nodes_open = Int64[]

    # Vector of PRESS nodes
    nodes_pressure = Int64[]

    # Vector of Probe nodes
    nodes_probe = Int64[]

    # Vector of Target nodes
    nodes_target = Int64[]

    # Vector with fixed elements
    elements_fixed = Int64[]

    # Vector with the values for the 
    # fixed elements
    values_fixed = Float64[]

    # Maximum id in materials
    max_id = 0

    # Loop over groups
    for g in LinearIndices(pgnames)

      # Name
      name = pgnames[g]

      # Split the string by ","
      st = split(name,",")

      # Check if Material
      if occursin("Material",st[1])

            # Clean dictionary to store local data
            empty!(localD_m)

            # We expect the following data
            # name, id, dens, c, Z
            id   = parse(Int64,st[3])
            dens = parse(Float64,st[4])
            c    = parse(Float64,st[5])
            Z    = parse(Float64,st[6])

            # Populate local dict
            localD_m["id"]   = id
            localD_m["dens"] = dens
            localD_m["c"]    = c
            localD_m["Z"]    = Z
          
            # Store maximum id
            max_id = max(max_id,id)

            # Now we must find wich elements are associated to this group
            elems_domain = Lgmsh.Readelementsgroup(meshfile,name,etags)

            # And store in the dict in Forward order (smaller to the largest)
            localD_m["elements"] = sort(elems_domain)

            # Copy the dict to the vector of materials
            push!(materials,copy(localD_m))

      # Check if Fixed
      elseif occursin("Fixed",st[1])

            # Valor a ser fixado
            value   = parse(Float64,st[2])
                
            # Now we must find wich elements are associated to this group
            elems_domain = Lgmsh.Readelementsgroup(meshfile,name,etags)

            # Adiciona os elementos ao vetor de elementos fixos
            elements_fixed = vcat(elements_fixed,sort(elems_domain))

            # E grava os valores fixos na sequência
            fixos = value*ones(length(elems_domain))

            # E concatena com o vetor de valores dos elementos fixos
            values_fixed = vcat(values_fixed,fixos)

      elseif  occursin("Open",st[1])

            # Find nodes 
            nodes = Lgmsh.Readnodesgroup(meshfile,name)

            # Append
            nodes_open = vcat(nodes_open,nodes)

      elseif  occursin("Probe",st[1])

            # Find nodes 
            nodes = Lgmsh.Readnodesgroup(meshfile,name)

            # Append
            nodes_probe = vcat(nodes_probe,nodes)

       elseif  occursin("Target",st[1])

            # Find nodes 
            nodes = Lgmsh.Readnodesgroup(meshfile,name)

            # Append
            nodes_target = vcat(nodes_target,nodes)

      elseif  occursin("Vn",st[1])

            # Clean dictionary to store local data
            empty!(localD_vn)

            # Normal velocity, frequency and phase 
            localD_vn["value"] = parse(Float64,st[2])
            localD_vn["freq"]  = parse(Float64,st[3])
            localD_vn["phase"] = parse(Float64,st[4])

            # Find nodes 
            nodes_vn = Lgmsh.Readnodesgroup(meshfile,name)

            # If 2D  - Find element and edges
            # else   - Find element faces
            # Vamos continuar chamando de edges, mesmo em 3D
            eleedges = Int64[]
            edges = Int64[]
            for tt in et
                if dimensao==2
                   eleedges_,edges_ = FindElementsEdges(tt,ne,etypes,connect,nodes_vn)
                else
                    eleedges_,edges_ = FindElementsFaces(tt,ne,etypes,connect,nodes_vn)
                end
                if !isempty(eleedges_)
                    push!(eleedges,eleedges_...)
                    push!(edges,edges_...)
                end
            end

            # Append
            localD_vn["elements"] = [eleedges edges]

            # Copy the dict to the vector of velocities
            push!(velocities,copy(localD_vn))

        elseif  occursin("Pressure",st[1])

            # Clean dictionary to store local data
            empty!(localD_pressure)

            # Pressure, frequency and phase 
            localD_pressure["value"] = parse(Float64,st[2])
            localD_pressure["freq"]  = parse(Float64,st[3])
            localD_pressure["phase"] = parse(Float64,st[4])

            # Find nodes 
            nodes_pressure_local = Lgmsh.Readnodesgroup(meshfile,name)

            # Append
            localD_pressure["nodes"] = nodes_pressure_local

            # Copy the dict to the vector of pressures
            push!(pressures,copy(localD_pressure))

            # Append nodes_pressure_local no vetor com TODOS os nodes_pressure
            nodes_pressure = vcat(nodes_pressure,nodes_pressure_local)

      elseif  occursin("Yn",st[1])

            # Clean dictionary to store local data
            empty!(localD_damp)

            # Valor
            localD_damp["value"] = parse(Float64,st[2])
            
            # Find nodes 
            nodes_damp = Lgmsh.Readnodesgroup(meshfile,name)

            # Find element and edges
            eleedges = Int64[]
            edges = Int64[]
            for tt in et
                if dimensao==2
                   eleedges_,edges_ = FindElementsEdges(tt,ne,etypes,connect,nodes_damp)
                else
                   eleedges_,edges_ = FindElementsFaces(tt,ne,etypes,connect,nodes_damp)
                end
                push!(eleedges,eleedges_...)
                push!(edges,edges_...)
            end
            
            # Append
            localD_damp["elements"] = [eleedges edges]

            # Copy the dict to the vector of dampings
            push!(damping,copy(localD_damp))


      end #if

    end 

    # We can now process data to build connectivities with material and
    # element types 
    connect2 = zeros(Int64,ne,nmax+2)

    # Some data are already processed
    connect2[:,1] .= etypes
    connect2[:,3:end] .= connect

    # Materials as a matrix
    materials2 = zeros(max_id,3)
   
    # loop over vector of material dicts
    for mat in materials

        # id
        id = mat["id"]

        # elements
        elements = mat["elements"]
        
        # Copy
        connect2[elements,2] .= id

        v = [mat["dens"] mat["c"] mat["Z"]]

        # fill line of materials2
        materials2[id,:] = v
         
    end

    # Testing...
    #=
    Mesh(nn, coord, ne, connect2, materials2, unique!(nodes_open), velocities, unique!(nodes_pressure), 
    pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed)
    =#

    # Return processed data
    return nn, coord, ne, connect2, materials2, unique!(nodes_open), velocities, unique!(nodes_pressure), 
           pressures, damping, nodes_probe, nodes_target, elements_fixed, values_fixed

end