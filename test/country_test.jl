
@testset "A Country" begin

	x,M,p = LandUse.runm()  # the baseline model
	xk,C,pk = LandUse.runk(estimateθ = false, istest = true)  # 2 identical regions in a country
	
	@test length(C[1].R) == 2
	@test isa(C[1].R[1], LandUse.Region)

	@testset "test vs single region" begin

		for it in 1:length(p.T)
			m = M[it]
			c = C[it]

			@test m.ϕ ≈ c.R[1].ϕ
			@test m.ϕ ≈ c.R[2].ϕ

			@test m.Lr ≈ c.R[1].Lr atol=1e-6
			@test m.Lr ≈ c.R[2].Lr atol=1e-6
			
			@test m.Lu ≈ c.R[1].Lu
			@test m.Lu ≈ c.R[2].Lu

			@test m.pr ≈ c.R[1].pr
			@test m.pr ≈ c.R[2].pr
			
			@test m.Sr ≈ c.R[1].Sr
			@test m.Sr ≈ c.R[2].Sr
			
			@test m.Srh ≈ c.R[1].Srh
			@test m.Srh ≈ c.R[2].Srh

			@test m.Srh + m.Sr + π * (m.ϕ)^2 ≈ p.S
			@test c.R[1].Srh + c.R[1].Sr + π * (c.R[1].ϕ)^2 ≈ c.Sk[1]
			@test c.R[1].Srh + c.R[1].Sr + π * (c.R[1].ϕ)^2 ≈ p.S
			@test c.R[2].Srh + c.R[2].Sr + π * (c.R[2].ϕ)^2 ≈ c.Sk[2]
			@test c.R[2].Srh + c.R[2].Sr + π * (c.R[2].ϕ)^2 ≈ p.S
		end
	end
	@testset "20 region country runs" begin
		K = 20
		x,C1,p1 = LandUse.k(K)
		@test length(C1[1].R) == K
	end
	# @testset "updating" begin
	# 	xs = rand(2*2)
	#     LandUse.update!(C,[p;p], vcat(M[1].Lr / M[1].Sr,M[1].r,M[1].pr,M[1].Sr, M[2].Sr,M[1].Lu, M[2].Lu ))
	# 	#
	# 	# 1. b: ratio of labor to land in region 1
	# 	# 2. r: land rent
	# 	# 3. pr: relative price rural good
	# 	# 4. 4 : (K+3), Sr: amount of land use in rural production (complement of Srh)
	# 	# 5. (K+4) : (2K+3), Lu: urban pop in each k

	#     @test C.R[1].pr == M[1].pr
	#     @test C.R[2].pr == M[1].pr

	#     @test C.R[1].Lr ≈ M[1].Lr
	#     # @test C.R[2].Lr ≈ M[2].Lr not true upon first update.
	#     @test C.R[1].Sr ≈ M[1].Sr
	#     @test C.R[2].Sr ≈ M[2].Sr
	#     @test C.R[1].Lu ≈ M[1].Lu
	#     @test C.R[2].Lu ≈ M[2].Lu

	# end

	# @testset "2 equal-region country" begin
	# 	sols,C,cpar,pp = LandUse.runk()

	# 	for it in 1:length(C[1].T)
	# 		@test LandUse.pop(C[it].R[1]) + LandUse.pop(C[it].R[2]) ≈ cpar.L atol=0.0001
	# 		@test LandUse.area(C[it].R[1]) + LandUse.area(C[it].R[2]) ≈ cpar.S atol=0.0001
	# 		@test LandUse.area(C[it].R[1]) ≈ C[it].R[1].ϕ + C[it].R[1].Sr + C[it].R[1].Srh atol=0.0001
	# 		@test LandUse.area(C[it].R[2]) ≈ C[it].R[2].ϕ + C[it].R[2].Sr + C[it].R[2].Srh atol=0.0001
	# 		@test LandUse.pop(C[it].R[1]) ≈ C[it].R[1].Lu + C[it].R[1].Lr atol=0.0001
	# 		@test LandUse.pop(C[it].R[2]) ≈ C[it].R[2].Lu + C[it].R[2].Lr atol=0.0001

	# 		# utility equal across countries
	# 		@test C[it].R[1].U == C[it].R[2].U
	# 	end
	# end

	# @testset "2 different-area country" begin
	# 	cpar = Dict(:S => 1.0, :L => 1.0, :kshare => [0.3,0.7], :K => 2)
	# 	sols,C,cpar,pp = LandUse.runk(cpar = cpar)
	# 	for it in 1:length(C[1].T)
	# 		@test LandUse.pop(C[it].R[1]) + LandUse.pop(C[it].R[2]) ≈ cpar.L atol=0.00001
	# 		@test LandUse.area(C[it].R[1]) + LandUse.area(C[it].R[2]) ≈ cpar.S atol=0.00001
	# 		@test LandUse.area(C[it].R[1]) ≈ C[it].R[1].ϕ + C[it].R[1].Sr + C[it].R[1].Srh atol=0.00001
	# 		@test LandUse.area(C[it].R[2]) ≈ C[it].R[2].ϕ + C[it].R[2].Sr + C[it].R[2].Srh atol=0.00001
	# 		@test LandUse.pop(C[it].R[1]) ≈ C[it].R[1].Lu + C[it].R[1].Lr atol=0.00001
	# 		@test LandUse.pop(C[it].R[2]) ≈ C[it].R[2].Lu + C[it].R[2].Lr atol=0.00001

	# 		# utility equal across countries
	# 		@test C[it].R[1].U ≈ C[it].R[2].U
	# 	end

	# end

	# @testset "2 different productivity same area" begin
	# 	cpar = Dict(:S => 1.0, :L => 1.0, :K => 2, :θprop => [1.0,0.99])
	# 	sols,C,cpar,pp = LandUse.runk(cpar = cpar)
	# 	for it in 1:length(C[1].T)
	# 		@test LandUse.pop(C[it].R[1]) + LandUse.pop(C[it].R[2]) ≈ cpar.L atol=0.00001
	# 		@test LandUse.area(C[it].R[1]) + LandUse.area(C[it].R[2]) ≈ cpar.S atol=0.00001
	# 		@test LandUse.area(C[it].R[1]) ≈ C[it].R[1].ϕ + C[it].R[1].Sr + C[it].R[1].Srh atol=0.00001
	# 		@test LandUse.area(C[it].R[2]) ≈ C[it].R[2].ϕ + C[it].R[2].Sr + C[it].R[2].Srh atol=0.00001
	# 		@test LandUse.pop(C[it].R[1]) ≈ C[it].R[1].Lu + C[it].R[1].Lr atol=0.00001
	# 		@test LandUse.pop(C[it].R[2]) ≈ C[it].R[2].Lu + C[it].R[2].Lr atol=0.00001

	# 		# utility equal across countries
	# 		@test C[it].R[1].U ≈ C[it].R[2].U
	# 	end
	# end

	# @testset "2 different prod and area" begin
	# 	cpar = Dict(:S => 1.0, :L => 1.0, :kshare => [0.7,0.3], :K => 2, :θg => [1.0,0.99])
	# 	sols,C,cpar,pp = LandUse.runk(cpar = cpar)
	# 	for it in 1:length(C[1].T)
	# 		@test LandUse.pop(C[it].R[1]) + LandUse.pop(C[it].R[2]) ≈ cpar.L atol=0.0001
	# 		@test LandUse.area(C[it].R[1]) + LandUse.area(C[it].R[2]) ≈ cpar.S atol=0.0001
	# 		@test LandUse.area(C[it].R[1]) ≈ C[it].R[1].ϕ + C[it].R[1].Sr + C[it].R[1].Srh atol=0.0001
	# 		@test LandUse.area(C[it].R[2]) ≈ C[it].R[2].ϕ + C[it].R[2].Sr + C[it].R[2].Srh atol=0.0001
	# 		@test LandUse.pop(C[it].R[1]) ≈ C[it].R[1].Lu + C[it].R[1].Lr atol=0.0001
	# 		@test LandUse.pop(C[it].R[2]) ≈ C[it].R[2].Lu + C[it].R[2].Lr atol=0.0001

	# 		# utility equal across countries
	# 		@test C[it].R[1].U ≈ C[it].R[2].U
	# 	end
	# end
	# @testset "3 different prod and area" begin
	# 	cpar = Dict(:S => 1.0, :L => 1.0, :kshare => [0.6,0.2,0.2], :K => 3, :θg => [1.0,1.01,1.02])
	# 	sols,C,cpar,pp = LandUse.runk(cpar = cpar)
	# 	for it in 1:length(C[1].T)
	# 		@test LandUse.pop(C[it].R[1]) + LandUse.pop(C[it].R[2]) + LandUse.pop(C[it].R[3]) ≈ cpar.L atol=0.00001
	# 		@test LandUse.area(C[it].R[1]) + LandUse.area(C[it].R[2]) + LandUse.area(C[it].R[3]) ≈ cpar.S atol=0.00001
	# 		@test LandUse.area(C[it].R[1]) ≈ C[it].R[1].ϕ + C[it].R[1].Sr + C[it].R[1].Srh atol=0.00001
	# 		@test LandUse.area(C[it].R[2]) ≈ C[it].R[2].ϕ + C[it].R[2].Sr + C[it].R[2].Srh atol=0.00001
	# 		@test LandUse.area(C[it].R[3]) ≈ C[it].R[3].ϕ + C[it].R[3].Sr + C[it].R[3].Srh atol=0.00001
	# 		@test LandUse.pop(C[it].R[1]) ≈ C[it].R[1].Lu + C[it].R[1].Lr atol=0.00001
	# 		@test LandUse.pop(C[it].R[2]) ≈ C[it].R[2].Lu + C[it].R[2].Lr atol=0.00001
	# 		@test LandUse.pop(C[it].R[3]) ≈ C[it].R[3].Lu + C[it].R[3].Lr atol=0.00001

	# 		# utility equal across countries
	# 		@test C[it].R[1].U ≈ C[it].R[2].U
	# 		@test C[it].R[1].U ≈ C[it].R[3].U
	# 	end
	# end

end
