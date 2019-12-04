module LandUse
	using JSON
	using NLsolve
	using FastGaussQuadrature

	const PEN = 10000.0  # penalty for nl solver

	include("param.jl")
	include("structchange.jl")
	include("model.jl")

	export Param, Model, CD0Model, StructChange!, solve!, update!

	function CD()
		p = Param()
		m = CD0Model(p)
		StructChange_closure(F,x) = StructChange!(F,x,p,m)
		r = nlsolve(StructChange_closure, [1; 0.5])
	end

	function main()
		# generate starting values
		p = LandUse.Param()
		m0 = LandUse.CD0Model(p)
		StructChange_closure(F,x) = LandUse.StructChange!(F,x,p,m0)
		r0 = LandUse.nlsolve(StructChange_closure, [1; 0.5])
		LandUse.update!(m0,p,r0.zero...)

		# use m0 values as starting values
		m = Model(p)
		update!(m,m0,p)
		F_closure(F,x) = solve!(F,x,p,m)
		r = nlsolve(F_closure,[m.qr;m.ϕ;m.r;m.Lr;m.pr;m.Sr])

		# 
	end

	







end # module
