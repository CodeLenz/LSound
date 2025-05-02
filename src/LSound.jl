module LSound

    # Imports
    using LinearAlgebra
    using SparseArrays
    using StaticArrays
    using LinearSolve
    using LinearMaps
    using ArnoldiMethod
    using ProgressMeter
    using DelimitedFiles
    using Lgmsh

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
    include("otim/parametrizacao.jl")
    include("otim/global_otim.jl")
    include("otim/main_otim.jl")
    include("otim/objetivo.jl")
    include("otim/sensibilidade.jl")
    include("otim/df.jl")

    export Analise

end
