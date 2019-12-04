using LandUse
using Test

@testset "LandUse.jl" begin
    @testset "struct change model" begin
		p = Param()
		m0 = CD0Model(p)
		StructChange_closure(F,x) = StructChange!(F,x,p,m0)
		r0 = LandUse.nlsolve(StructChange_closure, [1; 0.5])
    end
    @testset "model constructor updating" begin
		p = Param()
		m = Model(p)
		x = rand(6)

		@test isnan(m.Lu )
		@test isnan(m.wu0)
		@test isnan(m.wr )
		@test isnan(m.Srh)
		@test isnan(m.ur )
		
		update!(m,p,x)

		@test m.qr == x[1]   # land price in rural sector
		@test m.ϕ  == x[2]   # city size
		@test m.r  == x[3]   # land rent
		@test m.Lr == x[4]   # employment in rural sector
		@test m.pr == x[5]   # relative price rural good
		@test m.Sr == x[6]   # amount of land used in rural production

		@test m.Lu   == p.L - m.Lr   # employment in urban sector
		@test m.wu0  == LandUse.wu0(m.Lu,p)   # wage rate urban sector at city center (distance = 0)
		@test m.wr   == LandUse.wr(m.Lu,m.ϕ,p) # wage rate rural sector
		@test m.Srh  == LandUse.Srh(p,m)
		@test m.ur   == LandUse.ur(p,m)
    end
end
