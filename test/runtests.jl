using Test
using AstrodynamicalModels

# Just a compilation test for now...
@testset "Vector Field Function Calls" begin
    
    # Do these function calls work?
    @test R2BPVectorField(randn(6), randn(1), 0) isa Vector{<:Number}
    @test CR3BPVectorField(randn(6), randn(1), 0) isa Vector{<:Number}
    @test CR3BPWithSTMVectorField(randn(42), randn(1), 0) isa Vector{<:Number}

end