

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

function nlopt_solveFM(;p = Param(),x0=nothing)
	fm = FModel(p)
	opt = Opt(:LN_COBYLA,4)
	# opt = Opt(:LN_NELDERMEAD,4)
	opt.lower_bounds = fill(0.001,4)
	# opt.upper_bounds = [1.0,1.0,1.0,1.0]
	f0(x::Vector, grad::Vector) = 1.0
	opt.min_objective = f0  # fix at a constant function
	equality_constraint!(opt,(r,x,g) -> NLopt_wrap(r,x,g,fm,p), fill(1e-9,4))
	# m.r    = x[1]   # land rent
	# m.Lr   = x[2]   # employment in rural sector
	# m.pr   = x[3]   # relative price rural good
	# m.Sr   = x[4]   # amount of land used in rural production
	# (optf,optx,ret) = optimize(opt, [0.1 * p.S, 0.5 * p.L, 0.3775, 0.545 * p.S])
	if isnothing(x0)
		x0 = [0.1 * p.S, 0.5 * p.L, 0.3775, 0.545 * p.S]   # an almost arbitrary point in the interior domain.
	else
		@assert length(x0) == ndims(opt)
	end
	(optf,optx,ret) = optimize(opt, x0)
end








"""
	get_starts(p::Param)

Generate starting values for all years from the fixed ϵ model `FModel`.

In year 1:

1. Construct a `FModel` from a starting value `[0.1 * p.S, 0.5 * p.L, 0.3775, 0.545 * p.S]`
2. In subsequent years, take solution from 1. of previous year as starting value

"""
function get_starts(p::Param)
	# fm = LandUse.FModel(p)  # create a fixed elasticity model
	startvals = Vector{Float64}[]  # an empty array of vectors

	# 2. For each time period
	for it in 1:length(p.T)
		setperiod!(p, it)   # set period on param to it

		# if first year
		if it == 1
			x0 = nlopt_solveFM(p=p)
			if (x0[3] == :ROUNDOFF_LIMITED) | (x0[3] == :LIMITED)
				# update2!(fm,p,x0[2])
				push!(startvals, x0[2])
				# println("rural market clears with $(Rmk(fm,p))")
			else
				# p1 = plot(fm.Ftrace',ylim = (-5,5),title = "Ftrace")
				# p2 = plot(fm.xtrace',title = "xtrace",label = ["rho" "phi" "r" "Lr" "pr" "Sr"])
				# pl = plot(p1,p2,layout = (1,2))
				# savefig(pl,joinpath(@__DIR__,"..","images","Fmodel_trace.png"))
				error("first FModel not converged")
			end

		else  # in other years just start at previous solution
			x0 = nlopt_solveFM(p=p,x0 = startvals[it-1])
			if (x0[3] == :ROUNDOFF_LIMITED) | (x0[3] == :LIMITED)
				# update2!(fm,p,x0[2])
				push!(startvals, x0[2])
				# println("rural market clears with $(Rmk(fm,p))")
			else
				print(x0)
				error("FModel not converged in period $it")
			end
		end
	end
	return startvals
end



"""
	adapt_ϵ(p::Param,x0::Vector{Float64})

Adaptively increase slope coefficient ``s`` in elasticity of housing supply function [`ϵ`](@ref).
Starts from the first period solution of [`FModel`](@ref).
"""
function adapt_ϵ(p::Param,x0::Vector{Float64})

	m = Region(p)

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
