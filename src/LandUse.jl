module LandUse
	using JSON
	using NLsolve
	using FastGaussQuadrature
	using Plots
	using LaTeXStrings
	using Interact

	const PEN = 100.0  # penalty for nl solver

	import Base.show

	include("param.jl")
	include("structchange.jl")
	include("model.jl")
	include("plotter.jl")
	include("interact.jl")

	export Param, Model, CD0Model, StructChange!, solve!, update!

	function CD()
		p = Param()
		m = CD0Model(p)
		StructChange_closure(F,x) = StructChange!(F,x,p,m)
		r = nlsolve(StructChange_closure, [1; 0.5])
	end




	function model(;pars=Dict())
		p = LandUse.Param(par=pars)
		m0 = LandUse.CD0Model(p)
		CD0_closure(F,x) = LandUse.solve!(F,x,p,m0)
		r0 = LandUse.nlsolve(CD0_closure, [1; 0.5])
		LandUse.update!(m0,p,r0.zero...)

		# first set of starting values
		x0 = [m0.ρr;0.00005;m0.r;m0.Lr;m0.pr;m0.Sr]  # manually overriding to a smaller city

		# use StructChange values as starting values for a Cobb-Douglas pref model
		m1 = LandUse.FModel(p)  # create a fixed elasticity model
		F_closure(F,x) = LandUse.solve!(F,x,p,m1)
		r = LandUse.nlsolve(F_closure,x0,iterations = 1000)

		LandUse.update!(m,m0,p)

		# r = LandUse.nlsolve(F_closure,[m.ρr;m.ϕ;m.r;m.Lr;m.pr;m.Sr],iterations = 2)
		r = LandUse.mcpsolve(F_closure,[0.01,0.01,0.01,0.01,0.01,0.01], [1.0,1.0,1.0,1.0,100.0,1.0],[m.ρr;m.ϕ;m.r;m.Lr;m.pr;m.Sr],iterations = 1000)
		println(r)
		if !converged(r)
			@warn("model has not converged!")
		end
		# print(r)
		plot_static(m,p)
		# (p,m)
	end

	"function to generate starting values for all years"
	function get_starts()
		# 1. initialize parameter
		p = Param()
		startvals = Vector{Float64}[]  # an empty array of vectors

		# 2. For each time period
		for it in 1:length(p.T)
			setperiod!(p, it)   # set period on param to it

			# if first year
			if it == 1
				# 3. solve structural change model without space
				# --> use solution to that system, but replace ϕ with something very small.
				# --> produces starting value x00
				m0 = LandUse.CD0Model(p)
				r0 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,m0), [1; 0.5])
				LandUse.update!(m0,p,r0.zero...)
				x00 = [m.ρr; 0.00005; m.r; m.Lr; m.pr; m.Sr]   # set very small city!

				# 4. solve general model with fixed elasticity starting from x00
				# --> closed form solutions for integrals
				# --> produces starting value x0
				m1 = LandUse.FModel(p)  # create a fited elasticity model
				r1 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,m1),x0,iterations = 1000)
				if converged(r1)
					push!(startvals, r1.zero)
				else
					error("CD0Model not converged")
				end

				LandUse.update!(m1,m0,p)
			else  # in other years just start at previous solution
				r1 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,m1),startvals[it-1],iterations = 1000)
				if converged(r1)
					push!(startvals, r1.zero)
				else
					error("FModel not converged")
				end
			end
		end
		return startvals
	end

	



			# 5. solve general model with fixed elasticity from x0
			# overwrite x0 with results
			# store this x0 somewhere indexed by the current index for theta, i

		# end i


	end

	function main()

		# 1. initialize parameter

		# 2. For each value i of θ

			# if first value of θ
				# 3. solve structural change model without space
				# --> use solution to that system, but replace ϕ with something very small.
				# --> starting value x00

				# 4. solve general model with fixed elasticity from x00
				# --> closed form solutions for integrals
				# --> starting value x0
			# end if

			# 5. solve general model with fixed elasticity from x0
			# overwrite x0 with results
			# store this x0 somewhere indexed by the current index for theta, i

		# end i

		# 6. For each value i of θ

			# if i == 1

				# adaptive search for increasing slope in epsilon function
				# desired final slope is ϵtarget
				# this finds a feasible solution at 1860 values

				# for each ϵ in range(0,stop = ϵtarget, length = 10)
					# if first value
						# resolve general model with zero slope in epsilon. starting value x0[1] from above
						# i.e. you take the values for 1860
						# store as x 
					# else (remaining values)
						# solve general model with current epsilon slope but starting from previous solution in x 
						# then store solution in x 
					# end if 
				# end adaptive loop for ϵ

			# else (later periods)

				# solve general model with final desired slope on epsilon function
				# each time taking the solution of the previous iteration, x, as the starting value.
				# save solution in a suitable spot.

		# end i

	end


	function main2(;pars=Dict())
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
