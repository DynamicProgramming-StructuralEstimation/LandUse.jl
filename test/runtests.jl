using LandUse
using Test
using LinearAlgebra

ENV["GKSwstype"] = "100"  # hack to make GR run on gitlab CI

@testset "LandUse.jl" begin
	include("model_test.jl")
	include("country_test.jl")
end
