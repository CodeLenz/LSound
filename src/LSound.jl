module LSound

    # Imports
    using LinearAlgebra
    using SparseArrays
    using StaticArrays
    using Statistics
    using LinearSolve
    using LinearMaps
    using ArnoldiMethod
    using ProgressMeter
    using DelimitedFiles
    using Lgmsh
    using YAML
    using Gmsh


    # Main structure for mesh 
    include("mesh.jl")

    include("parse/parsemsh_daniele.jl")
    include("global/auxiliar.jl")
    include("global/ebc.jl")
    include("elements/quad.jl")
    include("elements/tri.jl")
    include("elements/hexa.jl")
    include("elements/tetra.jl")
    include("elements/pyramid.jl")
    include("global/global.jl")
    include("solvers/autovalores.jl")
    include("solvers/newmark.jl")
    include("global/vetor_P.jl")
    include("global/U0.jl")
    include("solvers/bathe.jl")
    include("solvers/modal.jl")
    include("main.jl")

   

    export Analise

end
