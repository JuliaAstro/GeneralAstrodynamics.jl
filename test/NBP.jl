"""
Tests for NBP dynamics.
"""
module NBPTests

using AstrodynamicalModels, ModelingToolkit, Test

@testset "NBP Model Constructors" begin
    NBSystem(10; stm=false)
    NBSystem(2; stm=true)
    @test true
end

@testset "NBP Function Constructors" begin
    NBFunction(10; stm=false)
    NBFunction(2; stm=true)
    @test true
end

end
