struct Mesh

    # Number of nodes
    nn::Int64

    # Coordinates 
    coord::Matrix{Float64}

    # Number of elements
    ne::Int64

    # Connectivities
    connect::Matrix{Int64}

    # Vector of Dictionaries with material data
    materials::Vector{Dict{String,Union{Float64,Int64,Vector{Int64}}}}

    # Vector with open nodes (zero pressure)
    nodes_open::Vector{Int64}

    # Vector with Dicts of normal velocities on nodes
    velocities = Vector{Dict{String,Union{Float64,Matrix{Int64}}}}

    # Vector of nodes with prescribed pressures 
    nodes_pressure::Vector{Int64}

    # Vector with Dicts of pressure on nodes
    pressures::Vector{Dict{String,Union{Float64,Vector{Int64}}}}

    # Damping in faces of elements
    damping::Vector{Dict{String,Union{Float64,Matrix{Int64}}}}

    # Vector of nodes to probe 
    nodes_probe::Vector{Int64}

    # Vector of target nodes (optimization)
    nodes_target::Vector{Int64}

    # Vector of fixed elements (optimization)
    elements_fixed::Vector{Int64}

    # Values tobe fixed (in element_fixed)
    values_fixed::Vector{Float64}

end