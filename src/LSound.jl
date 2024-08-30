module LSound

    # Imports
    using LinearAlgebra
    using SparseArrays
    using StaticArrays
    using LinearMaps
    using ArnoldiMethod
    using ProgressMeter
    using DelimitedFiles
    using Lgmsh

    include("parse/parsemsh_daniele.jl")
    include("global/aux.jl")
    include("elements/quad.jl")
    include("elements/tri.jl")
    include("global/global.jl")
    include("solvers/autovalores.jl")
    include("solvers/newmark.jl")
    include("global/vetor_P.jl")
    include("global/matriz_C.jl")
    include("global/U0.jl")
    include("solvers/bathe.jl")
    include("solvers/modal.jl")
    include("main.jl")

    export Analise

end
