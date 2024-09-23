@testset "elements/hexa.jl" begin

    meshfile = "msh/tube_hexa.msh"
  
    f,X = Analise(meshfile,:Modal,nev=4)
  
    @test isapprox(f[2],170.7,atol=1E-0) 
  
  end
  
  