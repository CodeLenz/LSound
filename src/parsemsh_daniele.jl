#
# Read a .msh file and process data to a specific format
#
# Daniele (Acoustics)
#

using Lgmsh

function Parsemsh_Daniele(meshfile::String)

    # Element types to  read
    et = [3]

    # Maximum number of nodes in the elements of the mesh
    nmax = maximum(Lgmsh_nodemap()[et])

    # Read mesh
    nn, coord, ne, etypes, connect, etags = Readmesh(meshfile,et)

    #
    # Expected Physical groups
    #
    # Material,nome,id,dens,c,Z [ surfaces (and volumes) ]
    #
    # Open [ lines and/or nodes]
    # 
    # Vn,value,freq,phase (in degrees)
    #
    #
    pgroups, pgnames = Lgmsh_import_physical_groups(meshfile)

    @show pgroups, pgnames

    # Vector with Dicts of materials
    materials = Dict{String,Union{Float64,Int64,Vector{Int64}}}[]

    # Local dict inside the loop
    localD_m = Dict{String,Union{Float64,Int64,Vector{Int64}}}()

    # Vector with Dicts of normal velocities on nodes
    velocities = Dict{String,Union{Float64,Matrix{Int64}}}[]

    # Local dict inside the loop
    localD_vn = Dict{String,Union{Float64,Matrix{Int64}}}()

    # Vector of OPEN nodes
    nodes_open = Int64[]

    # Loop over groups
    for g=1:length(pgnames)

      # Name
      name = pgnames[g]

      # Split the string by ","
      st = split(name,",")

      # Check if Material
      if occursin("Material",st[1])

            println("reading ", st)

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
          
            # Now we must find wich elements are associated to this group
            elems_domain = Lgmsh.Readelementsgroup(meshfile,name,etags)

            # And store in the dict in Forward order (smaller to the largest)
            localD_m["elements"] = sort(elems_domain)

            # Copy the dict to the vector of materials
            push!(materials,copy(localD_m))

      elseif  occursin("Open",st[1])

            # Find nodes 
            nodes = Lgmsh.Readnodesgroup(meshfile,name)

            # Append
            nodes_open = vcat(nodes_open,nodes)

      elseif  occursin("Vn",st[1])

            # Clean dictionary to store local data
            empty!(localD_vn)

            # Pressure, frequency and phase 
            localD_vn["value"] = parse(Float64,st[2])
            localD_vn["freq"]  = parse(Float64,st[3])
            localD_vn["phase"] = parse(Float64,st[4])

            # Find nodes 
            nodes_vn = Lgmsh.Readnodesgroup(meshfile,name)

            # Find element and edges
            eleedges,edges = FindElementsEdges(3,ne,etypes,connect,nodes_vn)

            # Append
            localD_vn["elements"] = [eleedges edges]

            # Copy the dict to the vector of velocities
            push!(velocities,copy(localD_vn))

      end #if

    end 

    # We can now process data to build connectivities with material and
    # element types 
    connect2 = zeros(Int64,ne,nmax+2)

    # Some data are already processed
    connect2[:,1] .= etypes
    connect2[:,3:end] .= connect

    # Materials as a matrix
    materials2 = zeros(length(materials),3)
   
    # loop over vector of material dicts
    for mat in materials

        # id
        id = mat["id"]

        # elements
        elements = mat["elements"]
        
        # Copy
        connect2[elements,2] .= id

        v = [mat["dens"] mat["c"] mat["Z"]]

        @show id, v

        # fill line of materials2
        materials2[id,:] = v
         
    end

    # Return processed data
    return nn, coord, ne, connect2, materials2, nodes_open, velocities

end