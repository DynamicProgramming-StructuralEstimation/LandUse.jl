

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
