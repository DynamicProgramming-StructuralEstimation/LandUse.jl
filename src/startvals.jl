
function nocommute!(F::Vector,x::Vector,p::Param)
	gamma2 = p.γ / (1 + p.ϵr)

	ρ = x[1]
	Lr = x[2]  

	Lu = p.L - Lr
	wu = p.θu*Lu^p.η
	wr = wu
	Sr = (((1 - p.α)/ p.α) * wr / ρ)^p.σ * Lr # farm land input
	r  = ρ * (p.S - p.λ) / p.L
	pr = wr / (p.α * p.θr) * (p.α + (1-p.α) * (Lr / Sr)^((1-p.σ)/p.σ))^(1/(1-p.σ))
	ϕ  = gamma2 * (wu + r - pr * p.cbar + p.sbar ) * Lu / (ρ) 
	Srh= gamma2 * (wr + r - pr * p.cbar + p.sbar ) * Lr / (ρ) 

	F[1] = (1 - p.γ) * (1 - p.ν) * (wr + r - pr * p.cbar + p.sbar ) + p.ϵr * ρ * (Srh + ϕ) - p.sbar * p.L - Lu * p.θu
	F[2] = Sr + Srh + ϕ - (p.S - p.λ)

end

function stmodel(p::Param)
	x0 = [1.0,0.5]
	x00 = nlsolve((F,x) -> nocommute!(F,x,p),x0)
	gamma2 = p.γ / (1 + p.ϵr)

	ρ  = x00.zero[1]
	Lr = x00.zero[2]  
	Lu = p.L - Lr
	wu = p.θu*Lu^p.η
	wr = wu
	Sr = (((1 - p.α)/ p.α) * wr / ρ)^p.σ * Lr # farm land input
	r  = ρ * (p.S - p.λ) / p.L
	pr = wr / (p.α * p.θr) * (p.α + (1-p.α) * (Lr / Sr)^((1-p.σ)/p.σ))^(1/(1-p.σ))
	ϕ  = gamma2 * (wu + r - pr * p.cbar + p.sbar ) * Lu / (ρ) 
	Srh= gamma2 * (wr + r - pr * p.cbar + p.sbar ) * Lr / (ρ) 

	(r = r, Lr = Lr, pr = pr, Sr = Sr)

end

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

startval(p::Param) = stmodel(p)

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
