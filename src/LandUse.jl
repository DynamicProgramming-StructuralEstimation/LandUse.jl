module LandUse

	const PEN = 10000.0

greet() = print("Hello World!")
	using NLsolve

	function param()
		p = (
		gama    = 1/4,   # housing weight
		epsilon = 2,     # housing supply elasticity
		eta     = 0.0,   # agglomeration forces
		nu      = 0.015, # weight of rural good consumption on consumption composite
		cbar    = 0.2,   # agr cons subsistence level (0.14 for thailand)
		sbar    = 0.0,   # neg urba subsistence level
		theta2  = 1,     # urban sector TFP
		theta1  = 1,     # farm sector TFP
		alpha   = 0.75,  # labor weight on farm sector production function
		lambda  = 0,     # useless land: non-farm, non-urban land (forests, national parks...)
		tau     = 8,     # commuting costs
		chi1    = 1,     # "easiness" to convert land into housing in rural sector
		chi2    = 1,     # "easiness" to convert land into housing in urban sector
		L       = 1,     # total population (note: total land normalized to 1)
		psi     = 1,     # urban ammenities relative to rural ammenities (note: creates wedge in urban-rural wage)
		times   = 1870:15:2035,
		sigma   = 0.99)  # land-labor elasticity of substitution in farm production function
		p
	end

	function Eqsys!(F,x,p)

		q1   = x[1]   # land price in rural sector
		L1   = x[2]   # employment in rural sector

		if (q1 < 0) || (L1 < 0)
			F[1] = PEN
			F[2] = PEN
		else
			L2   = p.L-L1   # employment in urban sector
			w2   = p.theta2 # wage rate urban sector with no commuting costs
			w1   = w2     # wage rate rural sector
			S1   = (1-p.alpha)/p.alpha*w1/q1*L1
			r    = 1/p.L*q1*(1-p.lambda)
			p1   = w1/(p.alpha*p.theta1)*(L1/S1)^(1-p.alpha)
			phi  = (L2/p.chi2)*(p.gama*(w2+r-p1*p.cbar)/q1)
			S1h  = (L1/p.chi1)*(p.gama*(w1+r-p1*p.cbar)/q1)
			
			F[1] = (1-p.gama)*(1-p.nu)*(w1+r-p1*p.cbar)-p.theta2*L2
			F[2] = S1+S1h+phi-(1-p.lambda)
		end
		
	end


	function main()
		p = param()
		Eqsys_closure(F,x) = Eqsys!(F,x,p)
		r = nlsolve(Eqsys_closure, [1; 0.5])
	end








end # module
