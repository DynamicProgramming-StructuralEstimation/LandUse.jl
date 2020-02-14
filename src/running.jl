

"""
	get_solutions()

Compute general model solutions for all years. Starts from solution
obtained for `t=1` and desired slope on elasticity function via [`adapt_ϵ`](@ref)
"""
function get_solutions2(x0::Vector{Float64},p::Param)
	# m = [Region(p) for it in 1:length(p.T)] # create a general elasticity model for each period
	m = Region[]
	sols = Vector{Float64}[]  # an empty array of vectors
	push!(sols, x0)  # first solution is obtained via `adapt_ϵ`

	# update t=1 model
	# setperiod!(p, 1)   # set period on param to it=1
	# update!(m[1],p,x0)

	# 2. For all periods
	for it in 1:length(p.T)
		# println("period $it")
		setperiod!(p, it)   # set period on param to it
		m0 = Region(p)

		r1 = LandUse.nlsolve((F,x) -> LandUse.solve2!(F,x,p,m0),
			                     sols[it],iterations = 1000)
		# r1 = LandUse.mcpsolve((F,x) -> LandUse.solve!(F,x,p,m),
		# 	                                        [0.01,0.01,0.01,0.01,0.01,0.01],
		# 	                                        [Inf,1.0,Inf,1.0,Inf,1.0],
		# 	                                        sols[it-1],iterations = 1000)
		if converged(r1)
			push!(sols, r1.zero)
			update2!(m0,p,r1.zero)
			@assert abs(Rmk(m0,p)) < 1e-7   # walras' law
			push!(m,m0)
		else
			error("General Model not converged in period $it")
		end
	end
	# final update

	return (sols, m)
end

function get_solutions(x0::Vector{Float64},p::Param)
	# m = [Region(p) for it in 1:length(p.T)] # create a general elasticity model for each period
	m = Region[]
	sols = Vector{Float64}[]  # an empty array of vectors
	push!(sols, x0)  # first solution is obtained via `adapt_ϵ`

	# update t=1 model
	# setperiod!(p, 1)   # set period on param to it=1
	# update!(m[1],p,x0)

	# 2. For all periods
	for it in 1:length(p.T)
		# println("period $it")
		setperiod!(p, it)   # set period on param to it
		m0 = Region(p)

		r1 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,m0),
			                     sols[it],iterations = 1000)
		# r1 = LandUse.mcpsolve((F,x) -> LandUse.solve!(F,x,p,m),
		# 	                                        [0.01,0.01,0.01,0.01,0.01,0.01],
		# 	                                        [Inf,1.0,Inf,1.0,Inf,1.0],
		# 	                                        sols[it-1],iterations = 1000)
		if converged(r1)
			push!(sols, r1.zero)
			update!(m0,p,r1.zero)
			@assert abs(Rmk(m0,p)) < 1e-7   # walras' law
			push!(m,m0)
		else
			error("General Model not converged in period $it")
		end
	end
	# final update

	return (sols, m)
end

function run(;par = Dict())

	x0 = get_starts(par=par)   # a T-array of starting vectors

	(x1,p) = adapt_ϵ(x0[1],par=par)  # adaptive search for higher epsilon in first period only

	x,M = get_solutions(x1[end],p)  # get general model solutions

	(x,M,p) # solutions, models, and parameter

end

function run2(;par = Dict())

	x0 = get_starts2(par=par)   # a T-array of starting vectors

	(x1,p) = adapt_ϵ2(x0[1],par=par)  # adaptive search for higher epsilon in first period only

	x,M = get_solutions2(x1[end],p)  # get general model solutions

	(x,M,p) # solutions, models, and parameter

end

function matlab_bm()
	println("starting values at same parameter vector:")
	display(vcat(LandUse.get_starts()'...))

	println("timing full model run")
	@time (x,M,p) = run();

end
