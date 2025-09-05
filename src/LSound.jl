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

   using JuMP
   using Alpine
   using Gurobi
   using Ipopt
   using HiGHS
   using Cbc

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

    # Otimização
    include("otim/le_yaml.jl")
    include("otim/LP.jl")
    include("otim/param_duhring.jl")
    include("otim/param_pereira.jl")
    include("otim/fixos.jl")
    include("otim/global_otim.jl")
    include("otim/objetivo.jl")
    include("otim/sensibilidade.jl")
    include("otim/df.jl")
    include("otim/volumes.jl")
    include("otim/centroides.jl")
    include("otim/vizinhanca.jl")
    include("otim/filtro_espacial.jl")
    include("otim/sweep.jl")
    include("otim/verifica_derivada.jl")
    include("otim/processa_FRF.jl")
    include("otim/perimiter.jl")
    include("otim/main_ISLP.jl")

    export Analise, Otim_ISLP, Processa_FRF

end
