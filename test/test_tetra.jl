@testset "elements/tetra.jl" begin

    meshfile = "msh/tube_tetra.msh"
  
    f,X = Analise(meshfile,:Modal,nev=4)
  
    @test isapprox(f[1],171.07,atol=1E-1) 
  
  end
  
  