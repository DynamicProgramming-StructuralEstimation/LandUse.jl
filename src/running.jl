

xtrace(x::NLsolve.SolverResults) = vcat([(x.trace[i].metadata["x"])' for i in 1:x.iterations]...)
ftrace(x::NLsolve.SolverResults) = vcat([(x.trace[i].metadata["f(x)"])' for i in 1:x.iterations]...)


function solve_once(p::Param,m0::Model,x0::Vector{Float64})
	nlsolve((F,x) -> solve!(F,x,p,m0),x0,iterations = p.iters,store_trace = p.trace, extended_trace = p.trace)
end


"""
run model for all time periods
"""
function run(p::Param; jump = true)

	setperiod!(p,1)
	x0 = startval(p)
	sols = NamedTuple[]
	push!(sols,x0)
	M = Region[]

	for it in 1:length(p.T)
		# println(it)
		setperiod!(p,it)
		m = Region(p)
		if jump
			x = jm(p,m,sols[it])
			push!(sols,x)
			update!(m,p,[x...])
		else
			x = solve_once(p,m,[sols[it]...])
			if converged(x)
				push!(sols,(; zip(keys(x0), x.zero)...))
				update!(m,p,x.zero)
			else
				error("nlsolve not converged in period $it")
			end
		end
		push!(M,m)
	end

	(sols[2:end],M,p) # solutions, models, and parameter

end

function runm(; jump = true)
	run(Param(), jump = jump)
end
function plot1()
	x,M,p = run(Param())
	ts_plots(M,p)
end




# archive


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


function solve1(p::Param)
	x = stmodel(p)
	# p.taum = 1
	# p.taul = 1
	m = Region(p)
	solve_once(p,m,[x...])
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
