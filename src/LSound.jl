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

    include("parse/parsemsh_daniele.jl")
    include("global/auxiliar.jl")
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

    # Otimização
    include("otim/le_yaml.jl")
    include("otim/parametrizacao.jl")
    include("otim/global_otim.jl")
    include("otim/main_otim.jl")
    include("otim/objetivo.jl")
    include("otim/sensibilidade.jl")
    include("otim/df.jl")
    include("otim/volumes.jl")
    include("otim/centroides.jl")
    include("otim/vizinhanca.jl")
    include("otim/filtro_espacial.jl")
    include("otim/BESO.jl")

    export Analise

end
