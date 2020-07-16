@testset "Single Region" begin
	@testset "components in flat model checks" begin

		p = LandUse.Param()   # flat epsilon slope, so closed form solutions apply
	    s = LandUse.get_starts(p)

		# use m0 values as starting values
		m = LandUse.Region(p)
		LandUse.update!(m,p,s[1])

		l = rand()
		if l >= m.ϕ
			@test LandUse.cost(l,m.ϕ,p) == p.c0 + p.c1 *m.ϕ + p.c2 * m.ϕ^2
		else
			@test LandUse.cost(l,m.ϕ,p) == p.c0 + p.c1 *l + p.c2 * l^2
		end

		# commuting cost
		@test LandUse.τ(0.4,0.5,p) > 0
		@test LandUse.τ(0.5,0.5,p) > 0
		@test LandUse.τ(0.0,m.ϕ,p) == 0


		# test wage function
		@test LandUse.w(m.Lu,0.0,m.ϕ,p) > LandUse.w(m.Lu,m.ϕ,m.ϕ,p)
		@test LandUse.w(m.Lu,0.0,m.ϕ,p) > LandUse.w(m.Lu,1.0,m.ϕ,p)
		@test LandUse.w(m.Lu,m.ϕ-eps(),m.ϕ,p) > LandUse.wu(m.Lu,m.ϕ,m.ϕ,p)
		@test LandUse.w(m.Lu,1.0,m.ϕ,p) == LandUse.w(m.Lu,m.ϕ,m.ϕ,p)
		@test LandUse.w(m.Lu,0.0,m.ϕ,p) == LandUse.wu(m.Lu,0.0,m.ϕ,p)
		@test LandUse.w(m.Lu,0.0,m.ϕ,p) == LandUse.wu0(m.Lu,p)
		@test LandUse.w(m.Lu,1.0,m.ϕ,p) < LandUse.wu0(m.Lu,p)

		# second version of wage function
		@test LandUse.w(m.Lu,0.0,m.ϕ,p) == LandUse.w(0.0,m.ϕ,p)
		@test LandUse.w(m.Lu,1.0,m.ϕ,p) == LandUse.w(1.0,m.ϕ,p)
		@test LandUse.w(m.Lu,m.ϕ-eps(),m.ϕ,p) == LandUse.w(m.ϕ-eps(),m.ϕ,p)

		# test optimal consumption rules
		@test LandUse.cu(0.0,p,m) == (1.0 - p.γ)*(1.0 - p.ν)*(LandUse.w(m.Lu,0.0,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) - p.sbar
		@test LandUse.cu(0.9,p,m) == (1.0 - p.γ)*(1.0 - p.ν)*(LandUse.w(m.Lu,0.9,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) - p.sbar
		@test LandUse.cr(0.0,p,m) == p.cbar + (1.0 - p.γ)*(p.ν)*(LandUse.w(m.Lu,0.0,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) ./ m.pr
		@test LandUse.cr(0.9,p,m) == p.cbar + (1.0 - p.γ)*(p.ν)*(LandUse.w(m.Lu,0.9,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) ./ m.pr

		# test excess consumption
		@test LandUse.wr(m.ϕ,p) == m.wr
		@test LandUse.xsr(p,m) ≈ m.wr + m.r + p.sbar - m.pr * p.cbar
		@test LandUse.xsu(0.0,p,m) == LandUse.w(0.0,m.ϕ,p) + m.r + p.sbar - m.pr * p.cbar
		@test LandUse.xsu(0.0,p,m) > LandUse.xsr(p,m)
		@test LandUse.xsu(m.ϕ,p,m) == LandUse.xsr(p,m)
		@test LandUse.xsu(1.0,p,m) == LandUse.xsr(p,m)

		# test house price function
		@test LandUse.q(0.0,p,m) >   LandUse.q(1.0,p,m)
		@test LandUse.q(0.0,p,m) >   m.qr
		@test LandUse.q(m.ϕ,p,m) ==  m.qr
		@test LandUse.q(1.0,p,m) ==  m.qr

		# test housing consumption
		@test LandUse.h(0.0,p,m) <   LandUse.h(m.ϕ,p,m)
		@test LandUse.h(1.0,p,m) ==  LandUse.h(m.ϕ,p,m)
		@test LandUse.h(1.0,p,m) ≈  p.γ * (m.wr + m.r + p.sbar - m.pr * p.cbar) / m.qr
		@test LandUse.Srh(p,m) ≈ m.Lr * LandUse.γ(m.ϕ,m.ϕ,p) * (m.wr + m.r + p.sbar - m.pr * p.cbar) / m.ρr

		# equation (17)
		@test LandUse.D(l,p,m) ≈ (LandUse.χ(l,m.ϕ,p) * LandUse.q(l,p,m)^(1+LandUse.ϵ(l,m.ϕ,p))) / (p.γ * (LandUse.w(m.Lu,l,m.ϕ,p) + m.r + p.sbar - m.pr * p.cbar))
		@test LandUse.D(l,p,m) ≈ LandUse.D2(l,p,m)

		@test LandUse.Yu(m,p) == p.θu * m.Lu
	end

	@testset "components in general model checks" begin

		p = LandUse.Param(par = Dict(:ϵsmax => 0.0))
		@test p.ϵs == 0.0
		p = LandUse.Param(par = Dict(:ϵsmax => 10.0))
		@test p.ϵs == 10.0

		p = LandUse.Param()

		s = LandUse.get_starts(p)

		# use m0 values as starting values
		m = LandUse.Region(p)
		LandUse.update!(m,p,s[1])

		l = rand()
		if l >= m.ϕ
			@test LandUse.cost(l,m.ϕ,p) == p.c0 + p.c1 *m.ϕ + p.c2 * m.ϕ^2
		else
			@test LandUse.cost(l,m.ϕ,p) == p.c0 + p.c1 *l + p.c2 * l^2
		end

		# commuting cost
		@test LandUse.τ(0.4,0.5,p) > 0
		@test LandUse.τ(0.5,0.5,p) > 0
		@test LandUse.τ(0.0,m.ϕ,p) == 0


		# test wage function
		@test LandUse.w(m.Lu,0.0,m.ϕ,p) > LandUse.w(m.Lu,m.ϕ,m.ϕ,p)
		@test LandUse.w(m.Lu,0.0,m.ϕ,p) > LandUse.w(m.Lu,1.0,m.ϕ,p)
		@test LandUse.w(m.Lu,m.ϕ-eps(),m.ϕ,p) > LandUse.wu(m.Lu,m.ϕ,m.ϕ,p)
		@test LandUse.w(m.Lu,1.0,m.ϕ,p) == LandUse.w(m.Lu,m.ϕ,m.ϕ,p)
		@test LandUse.w(m.Lu,0.0,m.ϕ,p) == LandUse.wu(m.Lu,0.0,m.ϕ,p)
		@test LandUse.w(m.Lu,0.0,m.ϕ,p) == LandUse.wu0(m.Lu,p)
		@test LandUse.w(m.Lu,1.0,m.ϕ,p) < LandUse.wu0(m.Lu,p)

		# test optimal consumption rules
		@test LandUse.cu(0.0,p,m) == (1.0 - p.γ)*(1.0 - p.ν)*(LandUse.w(m.Lu,0.0,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) - p.sbar
		@test LandUse.cu(0.9,p,m) == (1.0 - p.γ)*(1.0 - p.ν)*(LandUse.w(m.Lu,0.9,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) - p.sbar
		@test LandUse.cr(0.0,p,m) ≈ p.cbar + (1.0 - p.γ)*(p.ν)*(LandUse.w(m.Lu,0.0,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) ./ m.pr
		@test LandUse.cr(0.9,p,m) == p.cbar + (1.0 - p.γ)*(p.ν)*(LandUse.w(m.Lu,0.9,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) ./ m.pr


		# test excess consumption
		@test LandUse.xsr(p,m) ≈ m.wr + m.r + p.sbar - m.pr * p.cbar
		@test LandUse.xsu(0.0,p,m) == LandUse.w(0.0,m.ϕ,p) + m.r + p.sbar - m.pr * p.cbar
		@test LandUse.xsu(0.0,p,m) > LandUse.xsr(p,m)
		@test LandUse.xsu(m.ϕ,p,m) == LandUse.xsr(p,m)
		@test LandUse.xsu(1.0,p,m) == LandUse.xsr(p,m)

		# test house price function
		@test LandUse.q(0.0,p,m) >   LandUse.q(1.0,p,m)
		@test LandUse.q(0.0,p,m) >   m.qr
		@test LandUse.q(m.ϕ,p,m) ==  m.qr
		@test LandUse.q(1.0,p,m) ==  m.qr

		# test housing consumption
		@test LandUse.h(0.0,p,m) <   LandUse.h(m.ϕ,p,m)
		@test LandUse.h(1.0,p,m) ==  LandUse.h(m.ϕ,p,m)
		@test LandUse.h(1.0,p,m) ≈  p.γ * (m.wr + m.r + p.sbar - m.pr * p.cbar) / m.qr
		@test LandUse.h(0.0,p,m) ≈  p.γ * (m.wu0 + m.r + p.sbar - m.pr * p.cbar) / LandUse.q(0.0,p,m)
		@test LandUse.Srh(p,m) ≈ m.Lr * LandUse.γ(m.ϕ,m.ϕ,p) * (m.wr + m.r + p.sbar - m.pr * p.cbar) / m.ρr

		# equation (17)
		@test LandUse.D(l,p,m) ≈ (LandUse.χ(l,m.ϕ,p) * LandUse.q(l,p,m)^(1+LandUse.ϵ(l,m.ϕ,p))) / (p.γ * (LandUse.w(m.Lu,l,m.ϕ,p) + m.r + p.sbar - m.pr * p.cbar))

		@test LandUse.Yu(m,p) == p.θu * m.Lu
		# @test LandUse.D(l,p,m) ≈ LandUse.D2(l,p,m)

	end

	@testset "test flat elasticity ϵs = 0" begin

		p = LandUse.Param(par = Dict(:ϵsmax => 0.0))   # flat epsilon slope, so closed form solutions apply
		s = LandUse.get_starts(p)
		LandUse.setperiod!(p,1)  # make sure we are in year 1
		m = LandUse.Region(p)
		fm = LandUse.FModel(p)

		r = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,m),s[1],iterations = 100)
		rf = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,fm),s[1],iterations = 100)
		LandUse.update!(m,p,r.zero)
		LandUse.update!(fm,p,rf.zero)

		# we are testing results in the m object !
		# object fm is the object to test against.

		# is elasticity in fact constant?
		@test LandUse.ϵ(0.0,m.ϕ,p) == LandUse.ϵ(1.0,m.ϕ,p)
		@test LandUse.ϵ(0.0,m.ϕ,p) == p.ϵr
		@test LandUse.ϵ(m.ϕ,m.ϕ,p) == LandUse.ϵ(1.0,m.ϕ,p)
		@test LandUse.ϵ(m.ϕ,m.ϕ,p) == p.ϵr

		# consumption functinos

		x_analytic = zeros(6)
		x_general  = zeros(6)
		LandUse.Eqsys!(x_analytic,fm,p)
		LandUse.Eqsys!(x_general,m,p)

		@test norm(x_analytic .- x_general) < 1.e-8

		l = rand()
		@test LandUse.D(l,p,m) ≈ (LandUse.χ(l,m.ϕ,p) * LandUse.q(l,p,m)^(1+LandUse.ϵ(l,m.ϕ,p))) / (p.γ * (LandUse.w(m.Lu,l,m.ϕ,p) + m.r + p.sbar - m.pr * p.cbar))
		@test LandUse.D(l,p,m) ≈ LandUse.D2(l,p,m)

		# equation (23)
		@test m.r * p.L ≈ m.iq + m.ρr * (p.S - m.ϕ)

		# test equation (24)
		@test_broken m.r * p.L ≈ m.ρr + m.wu0 * m.iτ atol= p.S==1 ? 1e-3 : 2e-1


		# test per capita income formulation - equivalent in both formulations
		# pcy(m,p) = m.r + wr(m.Lu,m.ϕ,p) * m.Lr / p.L + m.iy / p.L
		@test LandUse.pcy(m,p) ≈ m.ρr / p.L + m.wu0 * (1.0 - LandUse.τ(m.ϕ,m.ϕ,p) * m.Lr / p.L) atol=0.03

		# test walras' law: also the rural goods market must clear at equilibrium:
		# equation (25)
		# @test p.ν * (1 - p.γ) * LandUse.pcy(m,p) + m.pr * p.cbar * (1.0 - p.ν * (1 - p.γ)) ==  m.pr * LandUse.Yr(m,p) / p.L
		@testset "Walras Law" begin
			@test isapprox(p.ν * (1 - p.γ) * (LandUse.pcy(m,p) + p.sbar - m.pr * p.cbar) + m.pr * p.cbar ,  m.pr * LandUse.Yr(m,p) / p.L )

			@test isapprox(p.ν * (1 - p.γ) * (LandUse.pcy(m,p) + p.sbar - m.pr * p.cbar) + m.pr * p.cbar - m.pr * LandUse.Yr(m,p) / p.L , LandUse.Rmk(m,p), atol = 1e-10)
			@test isapprox(LandUse.Rmk(m,p) , 0.0, atol = 1e-10)

			# try with different formulation of total rural demand:
			# rural cons in city + rural cons for each rural worker == total rural good production
			@test isapprox(m.icr + m.Lr * LandUse.cr(m.ϕ,p,m), LandUse.Yr(m,p))
		end

		@testset "equation (26) special case" begin

			@test isapprox( (1 - p.ν) * (1 - p.γ) * (LandUse.pcy(m,p) + p.sbar - m.pr * p.cbar) - p.sbar + p.ϵr * m.r , (1 - m.iτ) * LandUse.Yu(m,p) / m.Lu, atol = 0.05)
			# @test (1 - p.ν) * (1 - p.γ) * (LandUse.pcy(m,p) - m.pr * p.cbar) + p.ϵr * m.r == (m.Lu - m.iτ) * m.wu0 / p.L
			yy = fm.ρr / p.L + fm.wu0 * (1.0 - LandUse.τ(fm.ϕ,fm.ϕ,p) * fm.Lr / p.L)

			# the LHS seems correct (or at least consistently wrong across formulations :-)
			@test_broken (1 - p.ν) * (1 - p.γ) * (yy + p.sbar - fm.pr * p.cbar) + p.ϵr * fm.r ≈ (1 - p.ν) * (1 - p.γ) * (LandUse.pcy(m,p) + p.sbar - m.pr * p.cbar) + p.ϵr * m.r atol=0.02
		end


	end

	@testset "test full solution" begin
		p = LandUse.Param(par = Dict( :ηm => 1.1 ))

		x,M,p = LandUse.run(LandUse.Region,p)

		tol =  p.S == 1 ? 0.02 : 0.02
		# rutol = 1.0e-2
		# @warn("utility at precision $tol, rural cons $rutol")
		# for it in 1:length(p.T)
		for it in 1:length(p.T)
			if it > 8
				tol = 1e-1
				@warn "period $it on requires low tol=$tol" maxlog=1
			end
			@testset "period $it" begin
				LandUse.setperiod!(p,it)

				@test M[it].Sr ≈ p.S - M[it].ϕ - M[it].Srh
				@test M[it].Lu ≈ M[it].iDensity
				@test p.ϵr == LandUse.ϵ(1.0,M[it].ϕ,p)
				@test M[it].ρr ≈ (LandUse.χ(1.0,M[it].ϕ,p)/(1+LandUse.ϵ(1.0,M[it].ϕ,p))) * M[it].qr^(1+LandUse.ϵ(1.0,M[it].ϕ,p))
				# @test M[it].ρr ≈ (1-p.α) * M[it].pr * p.θr * (p.α * (M[it].Lr / M[it].Sr)^((p.σ-1)/p.σ) + (1-p.α))^(1/(p.σ-1))
				@test p.L  ≈ M[it].Lu + M[it].Lr

				# equation (23)
				@test M[it].r * p.L ≈ M[it].iq + M[it].ρr * (p.S - M[it].ϕ)

				# test equation (24)
				# no longer true with nonlinear commcost
				# @test isapprox(M[it].r * p.L , M[it].ρr + M[it].wu0 * M[it].iτ, atol = tol)

				# test per capita income formulation - equivalent in both formulations
				# pcy(m,p) = M[it].r + wr(M[it].Lu,M[it].ϕ,p) * M[it].Lr / p.L + M[it].iy / p.L
				@testset "per capita aggregate income" begin
					@test LandUse.τ(M[it].ϕ,M[it].ϕ,p) > 0.0
					# @test isapprox( LandUse.pcy(M[it],p) , M[it].ρr / p.L + M[it].wu0 * (1.0 - LandUse.τ(M[it].ϕ,M[it].ϕ,p) * M[it].Lr / p.L), atol = tol)
				end

				# test walras' law: also the rural goods market must clear at equilibrium:
				# equation (25)
				@testset "Walras' Law" begin
					@test p.ν * (1 - p.γ) * (LandUse.pcy(M[it],p) + p.sbar) + M[it].pr * p.cbar * (1.0 - p.ν * (1 - p.γ)) ≈  M[it].pr * LandUse.Yr(M[it],p) / p.L atol=tol
					# @test isapprox( p.ν * (1 - p.γ) * (LandUse.pcy(M[it],p) - M[it].pr * p.cbar) + M[it].pr * p.cbar , M[it].pr * LandUse.Yr(M[it],p) /	 p.L, atol = rutol)
					@test isapprox(p.ν * (1 - p.γ) * (LandUse.pcy(M[it],p) + p.sbar) + M[it].pr * p.cbar * (1.0 - p.ν * (1 - p.γ)) - M[it].pr * LandUse.Yr(M[it],p) / p.L  , 0.0, atol = tol)
					@test isapprox(LandUse.Rmk(M[it],p), 0.0,  atol = tol)

					# try with different formulation of total rural demand:
					# rural cons in city + rural cons for each rural worker == total rural good production
					@test isapprox(M[it].icr + M[it].Lr * LandUse.cr(M[it].ϕ,p,M[it]), LandUse.Yr(M[it],p), atol = tol)

				end

				# tight tol on utility
				l = rand()
				utol = 1e-9
				@test LandUse.D(l,p,M[it]) ≈ (LandUse.χ(l,M[it].ϕ,p) * LandUse.q(l,p,M[it])^(1+LandUse.ϵ(l,M[it].ϕ,p))) / (p.γ * (LandUse.w(M[it].Lu,l,M[it].ϕ,p) + M[it].r + p.sbar - M[it].pr * p.cbar))
				@test isapprox(LandUse.utility(1.0,p,M[it]), M[it].U, atol = utol)
				@test isapprox(LandUse.utility(0.1,p,M[it]), M[it].U, atol = utol)
				@test isapprox(LandUse.utility(M[it].ϕ,p,M[it]), M[it].U, atol = utol)
			end
		end
	end
end
