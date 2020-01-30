
@testset "A Country" begin
    p = LandUse.Param()
    @test_throws ArgumentError LandUse.Param(par = Dict(:K => 3, :Sk => rand(2)))
    @test_throws ArgumentError LandUse.Param(par = Dict(:Sk => [0.2,0.4]))

	C = LandUse.Country(p)
	@test length(C.R) == p.K
	@test isa(C.R[1], LandUse.Region)

end