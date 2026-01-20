"""
Tests for NBP dynamics.
"""
module NBPTests

using Test
using AstrodynamicalModels:
    NBSystem,
    NBFunction
using ModelingToolkit: System

@testset "NBP Model Constructors" begin
    @test NBSystem(10; stm=false) isa System
    @test NBSystem(2; stm=true) isa System
end

@testset "NBP Function Constructors" begin
    NBFunction(10; stm=false)
    NBFunction(2; stm=true)
    @test true
end

end
