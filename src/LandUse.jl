module LandUse
	using JSON
	using NLsolve
	using FastGaussQuadrature
	using Plots
	using LaTeXStrings

	const PEN = 10000.0  # penalty for nl solver

	import Base.show

	include("param.jl")
	include("structchange.jl")
	include("model.jl")
	include("plotter.jl")

	export Param, Model, CD0Model, StructChange!, solve!, update!

	function CD()
		p = Param()
		m = CD0Model(p)
		StructChange_closure(F,x) = StructChange!(F,x,p,m)
		r = nlsolve(StructChange_closure, [1; 0.5])
	end

	function model()
		p = LandUse.Param()
		m0 = LandUse.CD0Model(p)
		StructChange_closure(F,x) = LandUse.StructChange!(F,x,p,m0)
		r0 = LandUse.nlsolve(StructChange_closure, [1; 0.5])
		LandUse.update!(m0,p,r0.zero...)

		# use m0 values as starting values
		m = LandUse.Model(p)
		LandUse.update!(m,m0,p)
		# plot_static(m,p)
		(p,m)
	end


	function main(;pars=Dict())
		# generate starting values
		# p = LandUse.Param()
		p = LandUse.Param(par=pars)
		m0 = LandUse.CD0Model(p)
		StructChange_closure(F,x) = LandUse.StructChange!(F,x,p,m0)
		r0 = LandUse.nlsolve(StructChange_closure, [1; 0.5])
		LandUse.update!(m0,p,r0.zero...)

		# use m0 values as starting values
		m = LandUse.Model(p)
		LandUse.update!(m,m0,p)
		F_closure(F,x) = LandUse.solve!(F,x,p,m)
		r = LandUse.nlsolve(F_closure,[m.qr;m.ϕ;m.r;m.Lr;m.pr;m.Sr],iterations = 1000)
		println(r)
		res = r.zero
		println("qr = $(res[1])")
		println("ϕ = $(res[2])")
		println("r = $(res[3])")
		println("Lr = $(res[4])")
		println("pr = $(res[5])")
		println("Sr = $(res[6])")

		LandUse.update!(m,p,r.zero)
		return m

		# 
	end

	







end # module
