

xtrace(x::NLsolve.SolverResults) = vcat([(x.trace[i].metadata["x"])' for i in 1:x.iterations]...)
ftrace(x::NLsolve.SolverResults) = vcat([(x.trace[i].metadata["f(x)"])' for i in 1:x.iterations]...)


function solve_once(p::Param,m0::Model,x0::Vector{Float64})
	nlsolve((F,x) -> solve!(F,x,p,m0),x0,iterations = p.iters,store_trace = p.trace, extended_trace = p.trace)
end

function reduce_θ_step!(p::Param)
	# reduce step size
	p.θstep -= 0.05
	if p.θstep > 0 || error("reached neg step") end
end

function x0grid(m,p,x0,it)
	# make a grid over potential starting vectors: +- 5 % of current value
	# m.r    = x[1]   # land rent
	# m.Lr   = x[2]   # employment in rural sector
	# m.pr   = x[3]   # relative price rural good
	# m.Sr   = x[4]   # amount of land used in rural production

	r1 = solve_once(p,m,x0)
	if converged(r1)
		return r1.zero
	else
		println(it)
		n = 5
		ra = range(0.96,1.04,length = n)
		rp = range(1.00,1.05,length = n)
		rn = range(0.95,1.0,length = n)
		grid = vec(collect(Base.product(x0[1] .* rp, x0[2] .* rp,x0[3] .* rn,x0[4] .* rn)))

		for i in grid
			println([i...])
			r = solve_once(p,m,[i...])
			if converged(r)
				return r.zero
			end
			if p.trace
				traceplot(r,it)
			end
		end
		error("no solution found in entire grid")
	end
end

function trysolve(m,p,x0,it)

	# solution vector
	solutions = Vector[]
	push!(solutions,x0)
	ps   = Tuple{Float64,Float64}[]


	# targets
	θr1 = p.θr
	θu1 = p.θu
	r1 = solve_once(p,m,solutions[end])

	if converged(r1)
		return r1.zero
	else
		# get previous values
		setperiod!(p,it-1)
		θr0 = p.θr
		θu0 = p.θu
		println("enter while in period $it")

		while θu0!=θu1 || θr0!=θr1
			for current in (:urban, :rural)
				println(current)
				# current tuple (thetau, thetar)
				done = false
				step = 0.1
				while (!done)
					if current == :urban

						# p.θu = θu0 + step*(θu1-θu0)
						p.θu = θu0 + step
						# println("thetau0 = $θu0")
						# println("thetar0 = $θr0")
						println("thetau = $(p.θu)")
						r = solve_once(p,m,solutions[end])
						if converged(r)
							push!(solutions,r.zero)
							push!(ps,(p.θu,p.θr))
							done = true
							println(r)
							θu0 = p.θu
							println("done!")
							if p.trace
								traceplot(r,it)
							end
						else
							step = step * 0.9
						end


						# println("urban step = $step")
					else
						# p.θr = θr0 + step*(θr1-θr0)
						p.θr = θr0 + step

						# println("thetau0 = $θu0")
						# println("thetar0 = $θr0")
						println("thetar = $(p.θr)")
						r = solve_once(p,m,solutions[end])
						if converged(r)
							push!(solutions,r.zero)
							push!(ps,(p.θu,p.θr))
							done = true
							θr0 = p.θr
							println(r)


							println("done!")
							if p.trace
								traceplot(r,it)
							end
						else
							step = step * 0.9
						end
						# println("rural step = $step")
					end
				end
			end
		end
		return solutions[end]
	end

end

"""
in a given period and given a target param, this step-wise increases thetau
"""
function θstepper(p::Param,it::Int,m::Model,x0::Vector)

	sols = Vector[]

	push!(sols,x0)

	# target thetas
	if abs(p.θu - p.θr) > 0.78
		println("period $it")
		println("diff=$(p.θu - p.θr)")
		# go back to previous period, take p.θu, and set p.θu = p.θu_old + half the distance to the new obtained

		θus = range(p.θut[it-1],stop = p.θu, length = p.nsteps)
		θrs = range(p.θrt[it-1],stop = p.θr, length = p.nsteps)

		for i in 1:length(θus)
			p.θu = θus[i]
			p.θr = θrs[i]
			println("thetau = $(p.θu),thetar = $(p.θr)")
			r1 = solve_once(p,m,sols[i])
			if p.trace
				traceplot(r1,it)
			end
			if converged(r1)
				push!(sols, r1.zero)
			else
				println("last solution: $(sols[end])")

				error("adaptive searchrrrrrr not converged for thetau = $(p.θu),thetar = $(p.θr)")
			end
		end

		# if that converges, go to final values

		# else make step size smaller
	else
		r1 = solve_once(p,m,sols[1])
		if p.trace
			traceplot(r1,it)
		end
		if converged(r1)
			push!(sols, r1.zero)
		else
			error("not converged")
		end
	end


	return sols[end]
end

"""
	get_solutions()

Compute general model solutions for all years. Starts from solution
obtained for `t=1` and desired slope on elasticity function via [`adapt_ϵ`](@ref)
"""
function get_solutions(T::Type,x0::Vector{Float64},p::Param)
	# m = [Region(p) for it in 1:length(p.T)] # create a general elasticity model for each period
	m = T[]
	sols = Vector{Float64}[]  # an empty array of vectors
	push!(sols, x0)  # first solution is obtained via `adapt_ϵ`

	# update t=1 model
	# setperiod!(p, 1)   # set period on param to it=1
	# update!(m[1],p,x0)

	# 2. For all periods
	for it in 1:length(p.T)
	# for it in 1:20
		println("period $it")
		setperiod!(p, it)   # set period on param to it
		m0 = T(p)
		# println("tau = $(p.τ)")

		# x0 = nlopt_solve(p=p,x0 = sols[it])
		# if (x0[3] == :ROUNDOFF_LIMITED) | (x0[3] == :SUCCESS)
		# 	push!(sols, x0[2])
		# 	update!(m0,p,x0[2])
		# 	push!(m,m0)
		# else
		# 	error("type $T Model not converged in period $it")
		# end
		# m.r    = x[1]   # land rent
		# m.Lr   = x[2]   # employment in rural sector
		# m.pr   = x[3]   # relative price rural good
		# m.Sr   = x[4]   # amount of land used in rural production

		# nlopt solution
		# x0 = nlopt_solve(m0,p,sols[it])
		# if (x0[3] == :ROUNDOFF_LIMITED) | (x0[3] == :SUCCESS)
		# 	push!(sols, x0[2])
		# 	update!(m0,p,x0[2])
		# 	push!(m,m0)
		# else
		# 	error("type $T Model not converged in period $it")
		# end

		# x = θstepper(p,it,m0,sols[it])
		# push!(sols, x)
		# update!(m0,p,x)
		# push!(m,m0)

		# x = trysolve(m0,p,sols[it],it)
		# x = x0grid(m0,p,sols[it],it)


		# nlsolve solution
		r1 = solve_once(p,m0,sols[1])
		if p.trace
			traceplot(r1,it)
		end
		if converged(r1)
			push!(sols, r1.zero)
			update!(m0,p,r1.zero)
			# println("rmk = $(abs(Rmk(m0,p)))")
			# @assert abs(Rmk(m0,p)) < 1e-7   # walras' law
			push!(m,m0)
		else
			println(r1)
			error("not converged")
		end
	end
	return (sols, m)
end


# function run(;par = Dict())
function run(T::Type,p::Param)
	# x0 = get_starts(par=par)
	# x0 = get_starts(p)   # a T-array of starting vectors
	# r0 = nlsolve_starts(startval(),p=p)   # a nlsolve result object
	r0 = startval(p)
	x0 = [r0...]

	# if T == Urban
	# 	x0[1] = x0[1][1:3]
	# end

	# (x1,p) = adapt_ϵ(x0[1],par=par)
	# println("x0 = $(x0[1])")

	setperiod!(p,1)  # go back to period 1
	# if p.ϵsmax == 0.0
	# 	# p.ϵs = 0.0
	# 	x1 = x0[1]
	# 	@assert p.ϵs == 0.0
	# else
	# 	(xt,p) = adapt_ϵ(T(p),p,x0[1])  # adaptive search for higher epsilon in center first period only
	# 	x1 = xt[1]
	# end

	# if T == Urban
	# 	x1[end] = x1[end][1:3]
	# end
	# println("x1 = $(x1[end])")

	x,M = get_solutions(T,x0,p)  # get general model solutions

	(x,M,p) # solutions, models, and parameter

end

function get_urban_sols(x0::Vector{Float64},p::Param)
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
			# println("rmk = $(abs(Rmk(m0,p)))")
			@assert abs(Rmk(m0,p)) < 1e-7   # walras' law
			push!(m,m0)
		else
			error("General Model not converged in period $it")
		end
	end
	# final update

	return (sols, m)
end

function runm()
	run(Region,Param())
end


function matlab_bm()
	println("starting values at same parameter vector:")
	display(vcat(LandUse.get_starts()'...))

	println("timing full model run")
	@time (x,M,p) = run();

end

function make_space_gif()
	x,C,cpar,par = LandUse.runk()
	anim_space(C,par)
end

function make_ts_space_gif()
	x,C,cpar,par = LandUse.runk()
	plot_ts_xsect(C,par,1)
end

function plotsingle()
	p1  = Param() # baseline param: high cbar and low sbar
	x,M,p0  = run(Region,p1)
	LandUse.ts_plots(M,p1)
end

function run1()
	p = Param()
	# x0 = get_starts(par=par)  # nlopt
	r0 = nlsolve_starts(p=p)    # nlsolve
	x0 = [r0.zero]

	# if T == Urban
	# 	x0[1] = x0[1][1:3]
	# end

	# (x1,p) = adapt_ϵ(x0[1],par=par)
	# println("x0 = $(x0[1])")

	setperiod!(p,1)  # go back to period 1
	if p.ϵsmax == 0.0
		# p.ϵs = 0.0
		x1 = x0[1]
		@assert p.ϵs == 0.0
	else
		(xt,p) = adapt_ϵ(T(p),p,x0[1])  # adaptive search for higher epsilon in center first period only
		x1 = xt[1]
	end

	# if T == Urban
	# 	x1[end] = x1[end][1:3]
	# end
	# println("x1 = $(x1[end])")

	p.T = p.T[1:2]
	x,M = get_solutions(Region,x1,p)  # get general model solutions

	(x,M,p) # solutions, models, and parameter

end
