using Test

@testset "Aqua tests" begin
    using Aqua: test_all
    using AstrodynamicalModels
    test_all(AstrodynamicalModels; ambiguities = false)
end
