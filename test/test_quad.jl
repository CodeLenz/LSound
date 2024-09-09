@testset "elements/quad.jl" begin


    # Testa N 
    @test  all(LSound.Matriz_N_bi4(-1,-1).==[1.0 0.0 0.0 0.0])
    @test  all(LSound.Matriz_N_bi4( 1,-1).==[0.0 1.0 0.0 0.0])
    @test  all(LSound.Matriz_N_bi4( 1, 1).==[0.0 0.0 1.0 0.0])
    @test  all(LSound.Matriz_N_bi4(-1, 1).==[0.0 0.0 0.0 1.0])
    @test  all(LSound.Matriz_N_bi4( 0, 0).==[0.25 0.25 0.25 0.25])

    # Test dN
    dNr,dNs = LSound.dNrs_bi4(-1,-1)
    @test all(dNr .== [-0.5;0.5;0;0])
    @test all(dNs .== [-0.5;0;0;0.5])

    dNr,dNs = LSound.dNrs_bi4(1,-1)
    @test all(dNr .== [-0.5;0.5;0;0])
    @test all(dNs .== [ 0;-0.5;0.5;0])

    dNr,dNs = LSound.dNrs_bi4(1,1)
    @test all(dNr .== [0;0;0.5;-0.5])
    @test all(dNs .== [ 0;-0.5;0.5;0])
    
    dNr,dNs = LSound.dNrs_bi4(-1,1)
    @test all(dNr .== [0;0;0.5;-0.5])
    @test all(dNs .== [ -0.5;0;0;0.5])

    dNr,dNs = LSound.dNrs_bi4(0,0)
    @test all(dNr .== [-0.25; 0.25; 0.25; -0.25])
    @test all(dNs .== [-0.25; -0.25; 0.25; 0.25])


end

