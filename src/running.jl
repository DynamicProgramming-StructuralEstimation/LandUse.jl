

xtrace(x::NLsolve.SolverResults) = vcat([(x.trace[i].metadata["x"])' for i in 1:x.iterations]...)
ftrace(x::NLsolve.SolverResults) = vcat([(x.trace[i].metadata["f(x)"])' for i in 1:x.iterations]...)


function solve_once(p::Param,m0::Model,x0::Vector{Float64})
	nlsolve((F,x) -> solve!(F,x,p,m0),x0,iterations = p.iters,store_trace = p.trace, extended_trace = p.trace)
end


"""
run Single region model for all time periods
"""
function run(p::Param; jump = true, estimateθ = true)

	setperiod!(p,1)
	x0 = startval(p)

	# m = Region(p)
	# x = jm(p,m,x0)
	# px = x0.pr / p.moments[1,:P_rural]
	# transform!(p.moments,:P_rural => (x -> x .* px) => :P_rural)

	sols = NamedTuple[]
	push!(sols,x0)
	M = Region[]

	for it in 1:length(p.T)
		# println(it)
		setperiod!(p,it)
		m = Region(p)
		if jump
			x = jm(p,m,sols[it], estimateθ = estimateθ)
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
		if it == 1
			# adjust relative price in data to first period solution
			px = x.pr / p.moments[1,:P_rural]
			transform!(p.moments,:P_rural => (x -> x .* px) => :P_rural)
		end

	end

	(sols[2:end],M,p) # solutions, models, and parameter

end

"""
run Multi-region model for all time periods
"""
function runk(;par = Dict(:K => 2, :kshare => [0.5,0.5], :factors => [1.0,1.0]))

	# get single city solution in first period
	p = LandUse.Param(par = par, use_estimatedθ = true)
	@assert p.K > 1

	setperiod!(p,1)
	x0 = startval(p)
	m = Region(p)
	x0 = jm(p,m,x0)
	update!(m,p,[x0...])


	# starting value for country solution
	x = Float64[]
	push!(x, m.Lr / m.Sr)
	push!(x, m.r)
	push!(x, m.pr)
	for ik in 1:p.K
		push!(x,m.Sr)
	end
	for ik in 1:p.K
		push!(x,m.Lu)
	end

	sols = Vector{Float64}[]  # an empty array of solutions
	push!(sols,x)

	C = Country[]  # an emtpy array of countries

	for it in 1:length(p.T)
		# println(it)
		setperiod!(p,it)
		c = Country(p)
		x = jc(c,sols[it])
		push!(sols,x)
		update!(c,x)

		push!(C,c)
	end

	(sols,C,p) # solutions, models, and parameter

end

function runm(; jump = true, estimateθ = true)
	run(Param(), jump = jump,estimateθ = estimateθ)
end
function plot1()
	x,M,p = run(Param())
	ts_plots(M,p)
end
function export_thetas()
	x,M,p = runm()
	latex_param()
	d = DataFrame(year = collect(p.T), thetau = p.θut, thetar = p.θrt, pr = [M[it].pr for it in 1:length(M)])
	CSV.write(joinpath(dbtables,"export_theta_pr.csv"),d)

	x0 = LandUse.startval(Param())
	CSV.write(joinpath(dbtables,"export_x0.csv"),DataFrame([x0]))
end
