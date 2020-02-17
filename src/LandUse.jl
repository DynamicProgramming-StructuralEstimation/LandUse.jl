module LandUse

	# dependendencies
	using JSON
	using NLsolve
	using FastGaussQuadrature
	using Plots
	using LaTeXStrings
	# using Interact
	using DataFrames
	using StatsPlots
	using DataFramesMeta
	using Printf
	using NLopt

	# constants
	const PEN = 100.0  # penalty for nl solver
	const dbpath = joinpath(ENV["HOME"],"Dropbox","research","LandUse","output","model")

	# imports
	import Base.show, Base.convert

	# our code
	include("param.jl")
	include("model.jl")
	include("country.jl")
	include("startvals.jl")
	include("running.jl")
	include("plotter.jl")
	# include("interact.jl")

	global Xtrace = Vector{Float64}[]
	global Ftrace = Vector{Float64}[]



end # module
