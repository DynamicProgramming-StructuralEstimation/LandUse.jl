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
	include("model.jl")
	include("structchange.jl")
	include("plotter.jl")
	include("interact.jl")

	export Param, Model, CD0Model, StructChange!, solve!, update!





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

	"""
		get_starts()

	Generate starting values for all years.
	"""
	function get_starts()
		# 1. initialize parameter
		# p = Param(par=Dict(:ϵr => 0.0, :ϵs => 0.0))
		p = Param()
		fm = LandUse.FModel(p)  # create a fixed elasticity model
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
				r0 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,m0), [1; 0.5; 0.0])
				LandUse.update!(m0,p,r0.zero[1],r0.zero[2])
				x00 = [m0.ρr; 0.00005; m0.r; m0.Lr; m0.pr; m0.Sr; m0.U]   # set very small city!

				# 4. solve general model with fixed elasticity starting from x00
				# --> closed form solutions for integrals
				# --> produces starting value x0
				# r1 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,fm),x00,iterations = 1000, autodiff = :forward)
				# r1 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,fm),x00,iterations = 100, method = :trust_region,show_trace = true, extended_trace = true)
				r1 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,fm),x00,iterations = 100)
				if converged(r1)
					push!(startvals, r1.zero)
				else
					error("first FModel not converged")
				end

			else  # in other years just start at previous solution
				lb = zeros(7)
				ub = fill(Inf,7)
				r1 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,fm),startvals[it-1],iterations = 100)
				# r1 = LandUse.mcpsolve((F,x) -> LandUse.solve!(F,x,p,fm),lb,ub, startvals[it-1],iterations = 10000, factor = 0.1)
				if converged(r1)
					push!(startvals, r1.zero)
				else
					print(r1)
					error("FModel not converged")
				end
			end
		end
		return startvals

	end


	"""
		adapt_ϵ(x0::Vector{Float64})

	Adaptively increase slope coefficient ``s`` in elasticity of housing supply function [`ϵ`](@ref).
	Starts from the first period solution of [`FModel`](@ref).
	"""
	function adapt_ϵ(x0::Vector{Float64})

		p = Param()
		m = GModel(p)

		startvals = Vector{Float64}[]  # an empty array of vectors
		push!(startvals, x0)  # put 1860 solution for flat epsilon function

		# range of elasticity slope values
		ϵs = range(0,stop = p.ϵsmax, length = p.ϵnsteps)

		for (i,ϵ) in enumerate(ϵs)
			setfield!(p, :ϵs, ϵ)  # set current value for elaticity function slope

			r1 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,m),startvals[i],iterations = 1000)
			if converged(r1)
				push!(startvals, r1.zero)
			else
				error("adaptive search not converged for ϵ = $ϵ")
			end
		end
		return (startvals,p)
	end

	function run()

		x0 = get_starts()   # a T-array of starting vectors

		(x1,p) = adapt_ϵ(x0[1])  # adaptive search for higher epsilon in first period only

		x,M = get_solutions(x1[end],p)  # get general model solutions

		(x,M,p) # solutions, models, and parameter

	end

	function matlab_bm()
		println("starting values at same parameter vector:")
		display(vcat(LandUse.get_starts()'...))

		println("timing full model run")
		@time (x,M,p) = run();

	end


	"""
		get_solutions()

	Compute general model solutions for all years. Starts from solution
	obtained for `t=1` and desired slope on elasticity function via [`adapt_ϵ`](@ref)
	"""
	function get_solutions(x0::Vector{Float64},p::Param)
		m = [GModel(p) for it in 1:length(p.T)] # create a general elasticity model for each period

		sols = Vector{Float64}[]  # an empty array of vectors
		push!(sols, x0)  # first solution is obtained via `adapt_ϵ`

		# update t=1 model
		setperiod!(p, 1)   # set period on param to it=1
		update!(m[1],p,x0)

		# 2. For all periods starting at 2
		for it in 2:length(p.T)
			setperiod!(p, it)   # set period on param to it

			r1 = nlsolve((F,x) -> solve!(F,x,p,m[it]),
				                     sols[it-1],iterations = 1000)
			# r1 = LandUse.mcpsolve((F,x) -> LandUse.solve!(F,x,p,m),
			# 	                                        [0.01,0.01,0.01,0.01,0.01,0.01],
			# 	                                        [Inf,1.0,Inf,1.0,Inf,1.0],
			# 	                                        sols[it-1],iterations = 1000)
			if converged(r1)
				push!(sols, r1.zero)
				update!(m[it],p,r1.zero)
			else
				error("General Model not converged in period $it")
			end
		end
		# final update

		return (sols, m)
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
