
@testset "A Country" begin
    p = LandUse.Param()
    cp = LandUse.CParam()
    @test_throws ArgumentError LandUse.CParam(par = Dict(:K => 3, :kshare => rand(2)))
    @test_throws ArgumentError LandUse.CParam(par = Dict(:kshare => [0.2,0.4]))

	C = LandUse.Country(cp,[p;p])
	@test length(C.R) == cp.K
	@test isa(C.R[1], LandUse.Region)

	x,M,p = LandUse.run(p)  # get single region solutions

	@testset "test vs single region" begin


		@testset "density - city size" begin
		    @test LandUse.D(M[5].ϕ,p,M[5]) ≈ LandUse.D2(M[5].ϕ,p,M[5])
		    @test LandUse.D(M[5].ϕ,p,M[5]) ≈ LandUse.D2(M[5].ϕ,p,M[5])
		end

	end
	@testset "updating" begin
		xs = rand(2*2)
	    LandUse.update!(C,[p;p], vcat(M[1].Lr / M[1].Sr,M[1].r,M[1].pr,M[1].Sr, M[2].Sr ))

	    @test C.R[1].pr == M[1].pr
	    @test C.R[2].pr == M[1].pr

	    @test C.R[1].Lr ≈ M[1].Lr
	    @test C.R[2].Lr ≈ M[1].Lr atol=1e-2
	    @test C.R[1].Sr ≈ M[1].Sr
	    @test C.R[2].Sr ≈ M[1].Sr atol=1e-1

	end

end
