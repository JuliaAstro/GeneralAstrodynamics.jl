"""
Tests for NBP dynamics.
"""
module NBPTests

using AstrodynamicalModels, ModelingToolkit, Test

@testset "NBP Model Constructors" begin
    @test NBSystem(10; stm=false) isa ODESystem
    @test NBSystem(2; stm=true) isa ODESystem
end

@testset "NBP Function Constructors" begin
    NBFunction(10; stm=false)
    NBFunction(2; stm=true)
    @test true
end

end
