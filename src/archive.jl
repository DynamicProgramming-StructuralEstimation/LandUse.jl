
# archive



# new strategy to find starting values:
# use a constrained optimizer to avoid x < 0
# this works for the fixed model so far.
# julia> x = LandUse.nlopt_solve(par = Dict(:L => 100.0, :S => 20.0))
# (1.0, [0.04796005781257872, 0.4480817574692655, 0.3219286166291397, 14.974867748022646], :ROUNDOFF_LIMITED)


function NLopt_wrap(result::Vector, x::Vector, grad::Matrix,m::Model,p::Param)
	if length(grad) > 0
		# not implemented
	end
	solve!(result,x,p,m)
end


function nlopt_solve(m::Model,p::Param,x0::Vector{Float64})
	opt = Opt(:LN_COBYLA,length(x0))
	opt.lower_bounds = fill(0.001,length(x0))
	f0(x::Vector, grad::Vector) = 1.0
	opt.min_objective = f0  # fix at a constant function
	equality_constraint!(opt,(r,x,g) -> NLopt_wrap(r,x,g,m,p), fill(1e-9,length(x0)))
	# m.r    = x[1]   # land rent
	# m.Lr   = x[2]   # employment in rural sector
	# m.pr   = x[3]   # relative price rural good
	# m.Sr   = x[4]   # amount of land used in rural production
	# (optf,optx,ret) = optimize(opt, [0.1 * p.S, 0.5 * p.L, 0.3775, 0.545 * p.S])
	opt.upper_bounds = [Inf,1.0,x0[3]*1.5,1.0]

	(optf,optx,ret) = optimize(opt, x0)
end
#
# function nlopt_solve(;p = Param(),x0=nothing)
# 	fm = Region(p)
# 	opt = Opt(:LN_COBYLA,4)
# 	# opt = Opt(:LN_NELDERMEAD,4)
# 	opt.lower_bounds = fill(0.001,4)
# 	# opt.upper_bounds = [1.0,1.0,1.0,1.0]
# 	f0(x::Vector, grad::Vector) = 1.0
# 	opt.min_objective = f0  # fix at a constant function
# 	equality_constraint!(opt,(r,x,g) -> NLopt_wrap(r,x,g,fm,p), fill(1e-9,4))
# 	# m.r    = x[1]   # land rent
# 	# m.Lr   = x[2]   # employment in rural sector
# 	# m.pr   = x[3]   # relative price rural good
# 	# m.Sr   = x[4]   # amount of land used in rural production
# 	# (optf,optx,ret) = optimize(opt, [0.1 * p.S, 0.5 * p.L, 0.3775, 0.545 * p.S])
# 	if isnothing(x0)
# 		# x0 = [0.1 * p.S, 0.5 * p.L, 0.3775, 0.545 * p.S]   # an almost arbitrary point in the interior domain.
# 		x0 = [0.24, 0.68, 1.1, 0.86]   # an almost arbitrary point in the interior domain.
# 	else
# 		@assert length(x0) == ndims(opt)
# 		opt.upper_bounds = [Inf,1.0,x0[3]*1.5,1.0]
#
# 	end
# 	(optf,optx,ret) = optimize(opt, x0)
# end

function nlsolve_starts(x0;p = Param())
	m = Region(p)
	r1 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,m),
							 x0,iterations = 10000,store_trace = p.trace, extended_trace = p.trace)
	if !converged(r1)
		error("no starting values found")
	end
	r1
end








"""
	get_starts(p::Param)

Generate starting values for all years from the fixed ϵ model `FModel`.

In year 1:

1. Construct a `FModel` from a starting value `[0.1 * p.S, 0.5 * p.L, 0.3775, 0.545 * p.S]`
2. In subsequent years, take solution from 1. of previous year as starting value

"""
# function get_starts(;par = Dict())
function get_starts(p::Param)
	# fm = LandUse.FModel(p)  # create a fixed elasticity model
	# p = Param(par=par)
	starts = Vector{Float64}[]  # an empty array of vectors

	# 2. For each time period
	for it in 1:length(p.T)
		setperiod!(p, it)   # set period on param to it

		# if first year
		if it == 1
			x0 = nlopt_solve(startval(),p=p)
			if (x0[3] == :ROUNDOFF_LIMITED) | (x0[3] == :SUCCESS)
				# update2!(fm,p,x0[2])
				push!(starts, x0[2])
				# println("rural market clears with $(Rmk(fm,p))")
			else
				# p1 = plot(fm.Ftrace',ylim = (-5,5),title = "Ftrace")
				# p2 = plot(fm.xtrace',title = "xtrace",label = ["rho" "phi" "r" "Lr" "pr" "Sr"])
				# pl = plot(p1,p2,layout = (1,2))
				# savefig(pl,joinpath(@__DIR__,"..","images","Fmodel_trace.png"))
				error("first FModel not converged")
			end

		else  # in other years just start at previous solution
			x0 = nlopt_solve(p=p,x0 = startvals[it-1])
			if (x0[3] == :ROUNDOFF_LIMITED) | (x0[3] == :SUCCESS)
				# update2!(fm,p,x0[2])
				push!(starts, x0[2])
				# println("rural market clears with $(Rmk(fm,p))")
			else
				print(x0)
				error("FModel not converged in period $it")
			end
		end
	end
	return starts
end



"""
	adapt_ϵ(p::Param,x0::Vector{Float64})

Adaptively increase slope coefficient ``s`` in elasticity of housing supply function [`ϵ`](@ref).
Starts from the first period solution of [`FModel`](@ref).
"""
function adapt_ϵ(m::Model,p::Param,x0::Vector{Float64})
	# setperiod!(p,1)  # start in year one again
	# m = Region(p)

	startvals = Vector{Float64}[]  # an empty array of vectors
	push!(startvals, x0)  # put 1860 solution for flat epsilon function

	# range of elasticity slope values
	ϵs = range(0,stop = p.ϵsmax, length = p.ϵnsteps)[2:end]

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
