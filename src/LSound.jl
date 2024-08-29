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

    include("parsemsh_daniele.jl")
    include("elemento.jl")
    include("global.jl")
    include("autovalores.jl")
    include("newmark.jl")
    include("vetor_P.jl")
    include("matriz_C.jl")
    include("U0.jl")
    include("bathe.jl")

    include("modal.jl")
    include("main.jl")

    export Analise

end
