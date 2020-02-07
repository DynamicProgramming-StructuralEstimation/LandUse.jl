
@testset "A Country" begin
    p = LandUse.Param()
    @test_throws ArgumentError LandUse.Param(par = Dict(:K => 3, :Sk => rand(2)))
    @test_throws ArgumentError LandUse.Param(par = Dict(:Sk => [0.2,0.4]))

	C = LandUse.Country([p;p])
	@test length(C.R) == p.K
	@test isa(C.R[1], LandUse.Region)

	x,M,p = LandUse.run()  # get single region solutions

	@testset "test vs single region" begin


		@testset "density - city size" begin
		    @test LandUse.D(M[5].ϕ,p,M[5]) ≈ LandUse.D2(M[5].ϕ,p,M[5])
		    @test LandUse.D(M[5].ϕ,p,M[5]) ≈ LandUse.D2(M[5].ϕ,p,M[5])
		end

	end
	@testset "updating" begin
		xs = rand(2*2)
	    LandUse.update!(C,[p;p], vcat(M[1].ρr,M[1].wr,M[1].r, M[1].pr,xs))
	    @test C.R[1].pr == M[1].pr
	    @test C.R[1].ρr == M[1].ρr
	    @test C.R[2].pr == M[1].pr
	    @test C.R[2].ρr == M[1].ρr

	    @test C.R[1].Lr == xs[1]
	    @test C.R[2].Lr == xs[2]
	    @test C.R[1].Sr == xs[3]
	    @test C.R[2].Sr == xs[4]

	end

end

