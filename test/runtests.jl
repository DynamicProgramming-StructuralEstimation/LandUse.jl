using LandUse
using Test

@testset "LandUse.jl" begin
	p = Param()
	m0 = CD0Model(p)
	StructChange_closure(F,x) = StructChange!(F,x,p,m0)
	r0 = LandUse.nlsolve(StructChange_closure, [1; 0.5])
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

		p = LandUse.Param()
		m0 = LandUse.CD0Model(p)
		StructChange_closure(F,x) = LandUse.StructChange!(F,x,p,m0)
		r0 = LandUse.nlsolve(StructChange_closure, [1; 0.5])
		LandUse.update!(m0,p,r0.zero...)

		# use m0 values as starting values
		m = LandUse.Model(p)
		LandUse.update!(m,m0,p)

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
		@test LandUse.D(l,p,m) ≈ (LandUse.χ(l,m.ϕ,p) * LandUse.q(l,p,m)^(1+ϵ(l,m.ϕ,p))) / (p.γ * (LandUse.w(m.Lu,l,m.ϕ,p) + m.r - m.pr * p.cbar))

		# do integration
		# LandUse.integrate!(m,p)
		# @test m.iDensity == m.qr * 
    end
    @testset "converged model consistency" begin
		p = LandUse.Param()
		m0 = LandUse.CD0Model(p)
		StructChange_closure(F,x) = LandUse.StructChange!(F,x,p,m0)
		r0 = LandUse.nlsolve(StructChange_closure, [1; 0.5])
		LandUse.update!(m0,p,r0.zero...)

		# use m0 values as starting values
		m = LandUse.Model(p)
		LandUse.update!(m,m0,p)
		F_closure(F,x) = LandUse.solve!(F,x,p,m)
		r = LandUse.nlsolve(F_closure,[m.ρr;m.ϕ;m.r;m.Lr;m.pr;m.Sr],iterations = 1000)
		LandUse.update!(m,p,r.zero)  # final update 


		tol = 1e-4

		@test m.Sr ≈ 1.0 - m.ϕ - m.Srh
		@test m.Lu ≈ m.iDensity
		@test m.ρr ≈ (LandUse.χ(1,m.ϕ,p)/(1+ϵ(1,m.ϕ,p))) * m.qr^(1+ϵ(1,m.ϕ,p))
		@test m.ρr ≈ (1-p.α) * m.pr * p.θr * (p.α * (m.Lr / m.Sr)^((p.σ-1)/p.σ) + (1-p.α))^(1/(p.σ-1))
		@test p.L  ≈ m.Lu + m.Lr
		@test isapprox(LandUse.utility(1.0,p,m), m.U, atol = tol)
		@test isapprox(LandUse.utility(0.1,p,m), m.U, atol = tol)
		@test isapprox(LandUse.utility(m.ϕ,p,m), m.U, atol = tol)
        
    end
end
