
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

	@testset "2 equal-region country" begin
		sols,C,cpar,pp = LandUse.runk()

		for it in 1:length(C[1].T)
			@test LandUse.pop(C[it].R[1]) + LandUse.pop(C[it].R[2]) ≈ cpar.L
			@test LandUse.area(C[it].R[1]) + LandUse.area(C[it].R[2]) ≈ cpar.S
			@test LandUse.area(C[it].R[1]) ≈ C[it].R[1].ϕ + C[it].R[1].Sr + C[it].R[1].Srh
			@test LandUse.area(C[it].R[2]) ≈ C[it].R[2].ϕ + C[it].R[2].Sr + C[it].R[2].Srh
			@test LandUse.pop(C[it].R[1]) ≈ C[it].R[1].Lu + C[it].R[1].Lr
			@test LandUse.pop(C[it].R[2]) ≈ C[it].R[2].Lu + C[it].R[2].Lr

			# utility equal across countries
			@test C[it].R[1].U == C[it].R[2].U
		end
	end

	@testset "2 different-area country" begin
		cpar = Dict(:S => 1.0, :L => 1.0, :kshare => [0.3,0.7], :K => 2)
		sols,C,cpar,pp = LandUse.runk(cpar = cpar)
		for it in 1:length(C[1].T)
			@test LandUse.pop(C[it].R[1]) + LandUse.pop(C[it].R[2]) ≈ cpar.L
			@test LandUse.area(C[it].R[1]) + LandUse.area(C[it].R[2]) ≈ cpar.S
			@test LandUse.area(C[it].R[1]) ≈ C[it].R[1].ϕ + C[it].R[1].Sr + C[it].R[1].Srh
			@test LandUse.area(C[it].R[2]) ≈ C[it].R[2].ϕ + C[it].R[2].Sr + C[it].R[2].Srh
			@test LandUse.pop(C[it].R[1]) ≈ C[it].R[1].Lu + C[it].R[1].Lr
			@test LandUse.pop(C[it].R[2]) ≈ C[it].R[2].Lu + C[it].R[2].Lr

			# utility equal across countries
			@test C[it].R[1].U ≈ C[it].R[2].U
		end

	end

	@testset "2 different productivity same area" begin
		cpar = Dict(:S => 1.0, :L => 1.0, :K => 2, :θprop => [1.0,0.99])
		sols,C,cpar,pp = LandUse.runk(cpar = cpar)
		for it in 1:length(C[1].T)
			@test LandUse.pop(C[it].R[1]) + LandUse.pop(C[it].R[2]) ≈ cpar.L
			@test LandUse.area(C[it].R[1]) + LandUse.area(C[it].R[2]) ≈ cpar.S
			@test LandUse.area(C[it].R[1]) ≈ C[it].R[1].ϕ + C[it].R[1].Sr + C[it].R[1].Srh
			@test LandUse.area(C[it].R[2]) ≈ C[it].R[2].ϕ + C[it].R[2].Sr + C[it].R[2].Srh
			@test LandUse.pop(C[it].R[1]) ≈ C[it].R[1].Lu + C[it].R[1].Lr
			@test LandUse.pop(C[it].R[2]) ≈ C[it].R[2].Lu + C[it].R[2].Lr

			# utility equal across countries
			@test C[it].R[1].U ≈ C[it].R[2].U
		end
	end

	@testset "2 different prod and area" begin
		cpar = Dict(:S => 1.0, :L => 1.0, :kshare => [0.7,0.3], :K => 2, :θprop => [1.0,0.99])
		sols,C,cpar,pp = LandUse.runk(cpar = cpar)
		for it in 1:length(C[1].T)
			@test LandUse.pop(C[it].R[1]) + LandUse.pop(C[it].R[2]) ≈ cpar.L
			@test LandUse.area(C[it].R[1]) + LandUse.area(C[it].R[2]) ≈ cpar.S
			@test LandUse.area(C[it].R[1]) ≈ C[it].R[1].ϕ + C[it].R[1].Sr + C[it].R[1].Srh
			@test LandUse.area(C[it].R[2]) ≈ C[it].R[2].ϕ + C[it].R[2].Sr + C[it].R[2].Srh
			@test LandUse.pop(C[it].R[1]) ≈ C[it].R[1].Lu + C[it].R[1].Lr
			@test LandUse.pop(C[it].R[2]) ≈ C[it].R[2].Lu + C[it].R[2].Lr

			# utility equal across countries
			@test C[it].R[1].U ≈ C[it].R[2].U
		end
	end
	@testset "3 different prod and area" begin
		cpar = Dict(:S => 1.0, :L => 1.0, :kshare => [0.6,0.2,0.2], :K => 3, :θprop => [1.0,0.997,0.995])
		sols,C,cpar,pp = LandUse.runk(cpar = cpar)
		for it in 1:length(C[1].T)
			@test LandUse.pop(C[it].R[1]) + LandUse.pop(C[it].R[2]) + LandUse.pop(C[it].R[3]) ≈ cpar.L
			@test LandUse.area(C[it].R[1]) + LandUse.area(C[it].R[2]) + LandUse.area(C[it].R[3]) ≈ cpar.S
			@test LandUse.area(C[it].R[1]) ≈ C[it].R[1].ϕ + C[it].R[1].Sr + C[it].R[1].Srh
			@test LandUse.area(C[it].R[2]) ≈ C[it].R[2].ϕ + C[it].R[2].Sr + C[it].R[2].Srh
			@test LandUse.area(C[it].R[3]) ≈ C[it].R[3].ϕ + C[it].R[3].Sr + C[it].R[3].Srh
			@test LandUse.pop(C[it].R[1]) ≈ C[it].R[1].Lu + C[it].R[1].Lr
			@test LandUse.pop(C[it].R[2]) ≈ C[it].R[2].Lu + C[it].R[2].Lr
			@test LandUse.pop(C[it].R[3]) ≈ C[it].R[3].Lu + C[it].R[3].Lr

			# utility equal across countries
			@test C[it].R[1].U ≈ C[it].R[2].U
			@test C[it].R[1].U ≈ C[it].R[3].U
		end
	end

end
