using LandUse
using Test
using LinearAlgebra

@testset "LandUse.jl" begin
	p = LandUse.Param()
	m0 = LandUse.CD0Model(p)
	CD0_closure(F,x) = LandUse.solve!(F,x,p,m0)
	r0 = LandUse.nlsolve(CD0_closure, [1; 0.5])
    @testset "struct change model" begin

    end
    @testset "model constructor updating" begin
		LandUse.update!(m0,p,r0.zero...)

		m = LandUse.Model(p)
		@test isnan(m.Lu )
		@test isnan(m.wu0)
		@test isnan(m.wr )
		@test isnan(m.Srh)
		@test isnan(m.xsr )

		LandUse.update!(m,m0,p)

		@test m.ρr == m0.ρr   # land price in rural sector
		@test m.ϕ  == m0.ϕ    # city size
		@test m.r  == m0.r    # land rent
		@test m.Lr == m0.Lr   # employment in rural sector
		@test m.pr == m0.pr   # relative price rural good
		@test m.Sr == m0.Sr   # amount of land used in rural production

		@test m.Lu   == p.L - m.Lr   # employment in urban sector
		@test m.wu0  == LandUse.wu0(m.Lu,p)   # wage rate urban sector at city center (distance = 0)
		@test m.wr   == LandUse.wr(m.Lu,m.ϕ,p) # wage rate rural sector
		@test m.Srh  == LandUse.Srh(p,m)
		@test m.xsr   == LandUse.xsr(p,m)
    end

    @testset "consistency checks" begin

		p = LandUse.Param(par = Dict(:ϵs => 0.0))   # flat epsilon slope, so closed form solutions apply
        s = LandUse.get_starts()

		# use m0 values as starting values
		m = LandUse.Model(p)
		LandUse.update!(m,p,s[1])

		l = rand()
		if l >= m.ϕ
			@test LandUse.cost(l,m.ϕ,p) == p.c0 + p.c1 *m.ϕ + p.c2 * m.ϕ^2
		else
			@test LandUse.cost(l,m.ϕ,p) == p.c0 + p.c1 *l + p.c2 * l^2
		end


		# test wage function
		@test LandUse.w(m.Lu,0.0,m.ϕ,p) > LandUse.w(m.Lu,m.ϕ,m.ϕ,p)
		@test LandUse.w(m.Lu,0.0,m.ϕ,p) > LandUse.w(m.Lu,1.0,m.ϕ,p)
		@test LandUse.w(m.Lu,m.ϕ-eps(),m.ϕ,p) > LandUse.wu(m.Lu,m.ϕ,m.ϕ,p)
		@test LandUse.w(m.Lu,1.0,m.ϕ,p) == LandUse.w(m.Lu,m.ϕ,m.ϕ,p)
		@test LandUse.w(m.Lu,0.0,m.ϕ,p) == LandUse.wu(m.Lu,0.0,m.ϕ,p)
		@test LandUse.w(m.Lu,0.0,m.ϕ,p) == LandUse.wu0(m.Lu,p)
		@test LandUse.w(m.Lu,1.0,m.ϕ,p) < LandUse.wu0(m.Lu,p)

		# test excess consumption
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
		@test LandUse.h(1.0,p,m) ≈  p.γ * (m.wr + m.r - m.pr * p.cbar) / m.qr

		# equation (17)
		@test LandUse.D(l,p,m) ≈ (LandUse.χ(l,m.ϕ,p) * LandUse.q(l,p,m)^(1+LandUse.ϵ(l,m.ϕ,p))) / (p.γ * (LandUse.w(m.Lu,l,m.ϕ,p) + m.r - m.pr * p.cbar))
    end

    @testset "test solved model from starts" begin
        s = LandUse.get_starts()
		p = LandUse.Param(par = Dict(:ϵs => 0.0))   # flat epsilon slope, so closed form solutions apply
		m = LandUse.Model(p)
		fm = LandUse.FModel(p)


		r = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,m),s[1],iterations = 100)
		rf = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,fm),s[1],iterations = 100)
		LandUse.update!(m,p,r.zero)
		LandUse.update!(fm,p,rf.zero)

		# test integration
		# LandUse.integrate!(m,p)

		x_analytic = zeros(6)
		x_general  = zeros(6)
		LandUse.Eqsys!(x_analytic,fm,p) 
		LandUse.Eqsys!(x_general,m,p) 

		@test norm(x_analytic .- x_general) < 1.e-6

		l = rand()
		@test LandUse.D(l,p,m) ≈ (LandUse.χ(l,m.ϕ,p) * LandUse.q(l,p,m)^(1+LandUse.ϵ(l,m.ϕ,p))) / (p.γ * (LandUse.w(m.Lu,l,m.ϕ,p) + m.r - m.pr * p.cbar))

		# test equation (24)
		@test m.r * p.L ≈ m.ρr + m.wu0 * m.iτ

		# test per capita income formulation - equivalent in both formulations
		# pcy(m,p) = m.r + wr(m.Lu,m.ϕ,p) * m.Lr / p.L + m.iy / p.L
		@test LandUse.pcy(m,p) ≈ m.ρr / p.L + m.wu0 * (1.0 - LandUse.τ(m.ϕ,m.ϕ,p) * m.Lr / p.L)

		# test walras' law: also the rural goods market must clear at equilibrium:
		# equation (25)
		@test p.ν * (1 - p.γ) * LandUse.pcy(m,p) + m.pr * p.cbar * (1.0 - p.ν * (1 - p.γ)) == m.pr * LandUse.Yr(m,p) / p.L


    end

    @testset "converged model consistency" begin
    	x,M,p = LandUse.run()

		tol = 1e-4
    	for it in 1:length(p.T)
			@test M[it].Sr ≈ 1.0 - M[it].ϕ - M[it].Srh
			@test M[it].Lu ≈ M[it].iDensity
			@test p.ϵr == LandUse.ϵ(1.0,M[it].ϕ,p)
			@test M[it].ρr ≈ (LandUse.χ(1.0,M[it].ϕ,p)/(1+LandUse.ϵ(1.0,M[it].ϕ,p))) * M[it].qr^(1+LandUse.ϵ(1.0,M[it].ϕ,p))
			# @test M[it].ρr ≈ (1-p.α) * M[it].pr * p.θr * (p.α * (M[it].Lr / M[it].Sr)^((p.σ-1)/p.σ) + (1-p.α))^(1/(p.σ-1))
			@test p.L  ≈ M[it].Lu + M[it].Lr
			l = rand()
			@test LandUse.D(l,p,M[it]) ≈ (LandUse.χ(l,M[it].ϕ,p) * LandUse.q(l,p,M[it])^(1+LandUse.ϵ(l,M[it].ϕ,p))) / (p.γ * (LandUse.w(M[it].Lu,l,M[it].ϕ,p) + M[it].r - M[it].pr * p.cbar))
			@test_broken isapprox(LandUse.utility(1.0,p,M[it]), M[it].U, atol = tol)
			@test_broken isapprox(LandUse.utility(0.1,p,M[it]), M[it].U, atol = tol)
			@test_broken isapprox(LandUse.utility(M[it].ϕ,p,M[it]), M[it].U, atol = tol)
    	end
    end
end
