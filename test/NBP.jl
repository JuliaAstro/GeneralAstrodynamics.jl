"""
Tests for NBP dynamics.
"""
module NBPTests

using AstrodynamicalModels, ModelingToolkit, Test

@testset "NBP Model Constructors" begin
    NBP(10; stm=false)
    NBP(2; stm=true)
    @test true
end

@testset "NBP Function Constructors" begin
    NBPFunction(10; stm=false)
    NBPFunction(2; stm=true)
    @test true
end

end
