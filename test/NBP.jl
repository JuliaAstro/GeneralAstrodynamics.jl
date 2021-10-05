"""
Tests for NBP dynamics.
"""
module NBPTests

    using AstrodynamicalModels, ModelingToolkit, Test

    @testset "NBP Model Constructors" begin
        NBP(10; stm=false, structural_simplify=false) 
        NBP(2; stm=true, structural_simplify=true) 
        @test true
    end

    @testset "NBP Function Constructors" begin
        NBPFunction(10; stm=false, structural_simplify=false) 
        NBPFunction(2; stm=true, structural_simplify=true) 
        @test true
    end

end