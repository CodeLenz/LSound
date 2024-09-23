@testset "elements/quad.jl" begin

  meshfile = "msh/tube_quad.msh"

  f,X = Analise(meshfile,:Modal,nev=4)

  @test isapprox(f[1],170,atol=1E-1) 

end

