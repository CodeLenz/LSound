@testset "elements/tri.jl" begin

    meshfile = "msh/tube_tri.msh"
  
    f,X = Analise(meshfile,:Modal,nev=4)
  
    @test isapprox(f[1],170,atol=1E-1) 
  
  end
  
  