
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
	    LandUse.update!(C,[p;p], vcat(M[1].Lr / M[1].Sr,M[1].r,M[1].pr,M[1].Sr, M[2].Sr,M[1].Lu, M[2].Lu ))
		#
		# 1. b: ratio of labor to land in region 1
		# 2. r: land rent
		# 3. pr: relative price rural good
		# 4. 4 : (K+3), Sr: amount of land use in rural production (complement of Srh)
		# 5. (K+4) : (2K+3), Lu: urban pop in each k

	    @test C.R[1].pr == M[1].pr
	    @test C.R[2].pr == M[1].pr

	    @test C.R[1].Lr ≈ M[1].Lr
	    # @test C.R[2].Lr ≈ M[2].Lr not true upon first update.
	    @test C.R[1].Sr ≈ M[1].Sr
	    @test C.R[2].Sr ≈ M[2].Sr
	    @test C.R[1].Lu ≈ M[1].Lu
	    @test C.R[2].Lu ≈ M[2].Lu

	end

	@testset "2 region country" begin
		p = LandUse.Param()
		cp = LandUse.CParam()
		x = LandUse.runk()
		for it in 1:length(p.T)
			@test LandUse.pop(x[2][it].R[1]) + LandUse.pop(x[2][it].R[2]) ≈ cp.L
			@test LandUse.area(x[2][it].R[1]) + LandUse.area(x[2][it].R[2]) ≈ cp.S
			@test LandUse.area(x[2][it].R[1]) ≈ x[2][it].R[1].ϕ + x[2][it].R[1].Sr + x[2][it].R[1].Srh
			@test LandUse.area(x[2][it].R[2]) ≈ x[2][it].R[2].ϕ + x[2][it].R[2].Sr + x[2][it].R[2].Srh
			@test LandUse.pop(x[2][it].R[1]) ≈ x[2][it].R[1].Lu + x[2][it].R[1].Lr
			@test LandUse.pop(x[2][it].R[2]) ≈ x[2][it].R[2].Lu + x[2][it].R[2].Lr

		end
	end

end
