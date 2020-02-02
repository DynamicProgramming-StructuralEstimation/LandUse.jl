module LandUse

	# dependendencies
	using JSON
	using NLsolve
	using FastGaussQuadrature
	using Plots
	using LaTeXStrings
	using Interact
	using DataFrames
	using StatsPlots
	using DataFramesMeta

	# constants
	const PEN = 100.0  # penalty for nl solver
	const LU_CONST = 1.0

	# imports
	import Base.show

	# our code
	include("param.jl")
	include("model.jl")
	include("country.jl")
	include("startvals.jl")
	include("running.jl")
	include("plotter.jl")
	include("interact.jl")

end # module
