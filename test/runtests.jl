using LandUse
using Flux
using Test
using LinearAlgebra

# ENV["GKSwstype"] = "100"  # hack to make GR run on gitlab CI

# # not running tests for proper CI setup at the moment
# # compiling all required packages takes too long.
# # instead: https://discourse.julialang.org/t/running-package-tests-really-slow/46587/5
# # 1. make sure using Revise
# # 2. includet("test/runtests.jl")
# # 3. modify test functions in this file and keep running those.


# function t1()
# 	p = LandUse.Param()   # flat epsilon slope, so closed form solutions apply
# 	# s = LandUse.get_starts(p)
# 	x0 = LandUse.startval(p)

# 	# use m0 values as starting values
# 	m = LandUse.Region(p)
# 	LandUse.update!(m,p,[x0...])
	
# 	l = rand()
	
# 	# commuting cost
# 	@test LandUse.τ(0.4,p) > 0
# 	@test LandUse.τ(0.5,p) > 0
# 	@test LandUse.τ(0.0,p) == 0.0
	
	
# 	# test wage function
# 	@test LandUse.w(m.ϕ/2,m.ϕ,p) == p.θu - LandUse.τ(m.ϕ/2,p)
# 	@test LandUse.w(m.Lu,0.0,m.ϕ,p) > LandUse.w(m.Lu,m.ϕ,m.ϕ,p)
# 	@test LandUse.w(m.Lu,0.0,m.ϕ,p) > LandUse.w(m.Lu,1.0,m.ϕ,p)
# 	@test LandUse.w(m.Lu,m.ϕ-eps(),m.ϕ,p) > LandUse.wu(m.Lu,m.ϕ,m.ϕ,p)
# 	@test LandUse.w(m.Lu,1.0,m.ϕ,p) == LandUse.w(m.Lu,m.ϕ,m.ϕ,p)
# 	@test LandUse.w(m.Lu,0.0,m.ϕ,p) == LandUse.wu(m.Lu,0.0,m.ϕ,p)
# 	@test LandUse.w(m.Lu,0.0,m.ϕ,p) == LandUse.wu0(m.Lu,p)
# 	@test LandUse.w(m.Lu,1.0,m.ϕ,p) < LandUse.wu0(m.Lu,p)
	
# 	# second version of wage function
# 	@test LandUse.w(m.Lu,0.0,m.ϕ,p) == LandUse.w(0.0,m.ϕ,p)
# 	@test LandUse.w(m.Lu,1.0,m.ϕ,p) == LandUse.w(1.0,m.ϕ,p)
# 	@test LandUse.w(m.Lu,m.ϕ-eps(),m.ϕ,p) == LandUse.w(m.ϕ-eps(),m.ϕ,p)
	
# 	# test optimal consumption rules
# 	@test LandUse.cu(0.0,p,m) == (1.0 - p.γ)*(1.0 - p.ν)*(LandUse.w(m.Lu,0.0,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) - p.sbar
# 	@test LandUse.cu(0.9,p,m) == (1.0 - p.γ)*(1.0 - p.ν)*(LandUse.w(m.Lu,0.9,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) - p.sbar
# 	@test LandUse.cr(0.0,p,m) == p.cbar + (1.0 - p.γ)*(p.ν)*(LandUse.w(m.Lu,0.0,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) ./ m.pr
# 	@test LandUse.cr(0.9,p,m) == p.cbar + (1.0 - p.γ)*(p.ν)*(LandUse.w(m.Lu,0.9,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) ./ m.pr
	
# 	# test excess consumption
# 	@test LandUse.wr(m.ϕ,p) == m.wr
# 	@test LandUse.xsr(p,m) ≈ m.wr + m.r + p.sbar - m.pr * p.cbar
# 	@test LandUse.xsu(0.0,p,m) ≈ LandUse.w(0.0,m.ϕ,p) + m.r + p.sbar - m.pr * p.cbar
# 	@test LandUse.xsu(0.0,p,m) > LandUse.xsr(p,m)
# 	@test LandUse.xsu(m.ϕ,p,m) == LandUse.xsr(p,m)
# 	@test LandUse.xsu(1.0,p,m) == LandUse.xsr(p,m)
	
# 	# test house price function
# 	@test LandUse.q(0.0,p,m) >   LandUse.q(1.0,p,m)
# 	@test LandUse.q(0.0,p,m) >   m.qr
# 	@test LandUse.q(m.ϕ,p,m) ==  m.qr
# 	@test LandUse.q(1.0,p,m) ==  m.qr
	
# 	# test housing consumption
# 	@test LandUse.h(0.0,p,m) <   LandUse.h(m.ϕ,p,m)
# 	@test LandUse.h(1.0,p,m) ==  LandUse.h(m.ϕ,p,m)
# 	@test LandUse.h(1.0,p,m) ≈  p.γ * (m.wr + m.r + p.sbar - m.pr * p.cbar) / m.qr
# 	@test LandUse.Srh(p,m) ≈ m.Lr * LandUse.γ(m.ϕ,m.ϕ,p) * (m.wr + m.r + p.sbar - m.pr * p.cbar) / m.ρr
	
# 	# equation (17)
# 	@test LandUse.D(l,p,m) ≈ (LandUse.χ(l,m.ϕ,p) * LandUse.q(l,p,m)^(1+LandUse.ϵ(l,m.ϕ,p))) / (p.γ * (LandUse.w(m.Lu,l,m.ϕ,p) + m.r + p.sbar - m.pr * p.cbar))
# 	@test LandUse.D(l,p,m) ≈ LandUse.D2(l,p,m)
	
# 	@test LandUse.Yu(m,p) == p.θu * m.Lu
# end

# function t2()
# 	p = LandUse.Param(par = Dict(:ϵsmax => 0.0))
# 	@test p.ϵs == 0.0
# 	p = LandUse.Param(par = Dict(:ϵsmax => 10.0))
# 	@test p.ϵs == 10.0

# 	p = LandUse.Param()
# 	x0 = LandUse.startval(p)

# 	# s = LandUse.get_starts(p)

# 	# use m0 values as starting values
# 	m = LandUse.Region(p)
# 	LandUse.update!(m,p,[x0...])

# 	l = rand()

# 	# commuting cost
# 	@test LandUse.τ(0.4,p) > 0
# 	@test LandUse.τ(0.5,p) > 0
# 	@test LandUse.τ(0.0,p) == 0


# 	# test wage function
# 	@test LandUse.w(m.Lu,0.0,m.ϕ,p) > LandUse.w(m.Lu,m.ϕ,m.ϕ,p)
# 	@test LandUse.w(m.Lu,0.0,m.ϕ,p) > LandUse.w(m.Lu,1.0,m.ϕ,p)
# 	@test LandUse.w(m.Lu,m.ϕ-eps(),m.ϕ,p) > LandUse.wu(m.Lu,m.ϕ,m.ϕ,p)
# 	@test LandUse.w(m.Lu,1.0,m.ϕ,p) == LandUse.w(m.Lu,m.ϕ,m.ϕ,p)
# 	@test LandUse.w(m.Lu,0.0,m.ϕ,p) == LandUse.wu(m.Lu,0.0,m.ϕ,p)
# 	@test LandUse.w(m.Lu,0.0,m.ϕ,p) == LandUse.wu0(m.Lu,p)
# 	@test LandUse.w(m.Lu,1.0,m.ϕ,p) < LandUse.wu0(m.Lu,p)

# 	# test optimal consumption rules
# 	@test LandUse.cu(0.0,p,m) == (1.0 - p.γ)*(1.0 - p.ν)*(LandUse.w(m.Lu,0.0,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) - p.sbar
# 	@test LandUse.cu(0.9,p,m) == (1.0 - p.γ)*(1.0 - p.ν)*(LandUse.w(m.Lu,0.9,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) - p.sbar
# 	@test LandUse.cr(0.0,p,m) ≈ p.cbar + (1.0 - p.γ)*(p.ν)*(LandUse.w(m.Lu,0.0,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) ./ m.pr
# 	@test LandUse.cr(0.9,p,m) == p.cbar + (1.0 - p.γ)*(p.ν)*(LandUse.w(m.Lu,0.9,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) ./ m.pr


# 	# test excess consumption
# 	@test LandUse.xsr(p,m) ≈ m.wr + m.r + p.sbar - m.pr * p.cbar
# 	@test LandUse.xsu(0.0,p,m) ≈ LandUse.w(0.0,m.ϕ,p) + m.r + p.sbar - m.pr * p.cbar
# 	@test LandUse.xsu(0.0,p,m) > LandUse.xsr(p,m)
# 	@test LandUse.xsu(m.ϕ,p,m) == LandUse.xsr(p,m)
# 	@test LandUse.xsu(1.0,p,m) == LandUse.xsr(p,m)

# 	# test house price function
# 	@test LandUse.q(0.0,p,m) >   LandUse.q(1.0,p,m)
# 	@test LandUse.q(0.0,p,m) >   m.qr
# 	@test LandUse.q(m.ϕ,p,m) ==  m.qr
# 	@test LandUse.q(1.0,p,m) ==  m.qr

# 	# test housing consumption
# 	@test LandUse.h(0.0,p,m) <   LandUse.h(m.ϕ,p,m)
# 	@test LandUse.h(1.0,p,m) ==  LandUse.h(m.ϕ,p,m)
# 	@test LandUse.h(1.0,p,m) ≈  p.γ * (m.wr + m.r + p.sbar - m.pr * p.cbar) / m.qr
# 	@test LandUse.h(0.0,p,m) ≈  p.γ * (m.wu0 + m.r + p.sbar - m.pr * p.cbar) / LandUse.q(0.0,p,m)
# 	@test LandUse.Srh(p,m) ≈ m.Lr * LandUse.γ(m.ϕ,m.ϕ,p) * (m.wr + m.r + p.sbar - m.pr * p.cbar) / m.ρr

# 	# equation (17)
# 	@test LandUse.D(l,p,m) ≈ (LandUse.χ(l,m.ϕ,p) * LandUse.q(l,p,m)^(1+LandUse.ϵ(l,m.ϕ,p))) / (p.γ * (LandUse.w(m.Lu,l,m.ϕ,p) + m.r + p.sbar - m.pr * p.cbar))

# 	@test LandUse.Yu(m,p) == p.θu * m.Lu
# 	# @test LandUse.D(l,p,m) ≈ LandUse.D2(l,p,m)
# end


# function t3()
# 	p = LandUse.Param(par = Dict( :ηm => 1.1 ))
	
# 	x,M,p = LandUse.run(p)

# 	tol =  p.S == 1 ? 0.02 : 0.02
# 	# rutol = 1.0e-2
# 	# @warn("utility at precision $tol, rural cons $rutol")
# 	# for it in 1:length(p.T)
# 	for it in 1:length(p.T)
# 		if it > 8
# 			tol = 1e-1
# 			@warn "period $it on requires low tol=$tol" maxlog=1
# 		end
# 		@testset "period $it" begin
# 			LandUse.setperiod!(p,it)

# 			@test M[it].Sr ≈ p.S - M[it].ϕ - M[it].Srh
# 			@test M[it].Lu ≈ M[it].iDensity
# 			@test p.ϵr == LandUse.ϵ(1.0,M[it].ϕ,p)
# 			@test M[it].ρr ≈ (LandUse.χ(1.0,M[it].ϕ,p)/(1+LandUse.ϵ(1.0,M[it].ϕ,p))) * M[it].qr^(1+LandUse.ϵ(1.0,M[it].ϕ,p))
# 			# @test M[it].ρr ≈ (1-p.α) * M[it].pr * p.θr * (p.α * (M[it].Lr / M[it].Sr)^((p.σ-1)/p.σ) + (1-p.α))^(1/(p.σ-1))
# 			@test p.L  ≈ M[it].Lu + M[it].Lr

# 			# equation (23)
# 			@test M[it].r * p.L ≈ M[it].iq + M[it].ρr * (p.S - M[it].ϕ)

# 			# test equation (24)
# 			# no longer true with nonlinear commcost
# 			# @test isapprox(M[it].r * p.L , M[it].ρr + M[it].wu0 * M[it].iτ, atol = tol)

# 			# test per capita income formulation - equivalent in both formulations
# 			# pcy(m,p) = M[it].r + wr(M[it].Lu,M[it].ϕ,p) * M[it].Lr / p.L + M[it].iy / p.L
# 			@testset "per capita aggregate income" begin
# 				@test LandUse.τ(M[it].ϕ,p) > 0.0
# 				# @test isapprox( LandUse.pcy(M[it],p) , M[it].ρr / p.L + M[it].wu0 * (1.0 - LandUse.τ(M[it].ϕ,M[it].ϕ,p) * M[it].Lr / p.L), atol = tol)
# 			end

# 			# test walras' law: also the rural goods market must clear at equilibrium:
# 			# equation (25)
# 			@testset "Walras' Law" begin
# 				@test p.ν * (1 - p.γ) * (LandUse.pcy(M[it],p) + p.sbar) + M[it].pr * p.cbar * (1.0 - p.ν * (1 - p.γ)) ≈  M[it].pr * LandUse.Yr(M[it],p) / p.L atol=tol
# 				# @test isapprox( p.ν * (1 - p.γ) * (LandUse.pcy(M[it],p) - M[it].pr * p.cbar) + M[it].pr * p.cbar , M[it].pr * LandUse.Yr(M[it],p) /	 p.L, atol = rutol)
# 				@test isapprox(p.ν * (1 - p.γ) * (LandUse.pcy(M[it],p) + p.sbar) + M[it].pr * p.cbar * (1.0 - p.ν * (1 - p.γ)) - M[it].pr * LandUse.Yr(M[it],p) / p.L  , 0.0, atol = tol)
# 				@test isapprox(LandUse.Rmk(M[it],p), 0.0,  atol = tol)

# 				# try with different formulation of total rural demand:
# 				# rural cons in city + rural cons for each rural worker == total rural good production
# 				@test isapprox(M[it].icr + M[it].Lr * LandUse.cr(M[it].ϕ,p,M[it]), LandUse.Yr(M[it],p), atol = tol)

# 			end

# 			# tight tol on utility
# 			l = rand()
# 			utol = 1e-9
# 			@test LandUse.D(l,p,M[it]) ≈ (LandUse.χ(l,M[it].ϕ,p) * LandUse.q(l,p,M[it])^(1+LandUse.ϵ(l,M[it].ϕ,p))) / (p.γ * (LandUse.w(M[it].Lu,l,M[it].ϕ,p) + M[it].r + p.sbar - M[it].pr * p.cbar))
# 			@test isapprox(LandUse.utility(1.0,p,M[it]), M[it].U, atol = utol)
# 			@test isapprox(LandUse.utility(0.1,p,M[it]), M[it].U, atol = utol)
# 			@test isapprox(LandUse.utility(M[it].ϕ,p,M[it]), M[it].U, atol = utol)
# 		end
# 	end
# end

# function tall()
# 	t1()
# 	t2()
# 	t3()
# end

@testset "LandUse.jl" begin
	include("model_test.jl")
	include("country_test.jl")
end
