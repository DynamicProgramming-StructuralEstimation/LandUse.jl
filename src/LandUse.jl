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
	using LineSearches

	# constants
	const PEN = 100.0  # penalty for nl solver
	const dbpath = joinpath(ENV["HOME"],"Dropbox","research","LandUse","output","model")

	# imports
	import Base.show, Base.convert

	# globals
	Xtrace = Vector{Float64}[]
	Ftrace = Vector{Float64}[]

	# our code
	include("param.jl")
	include("model.jl")
	include("country.jl")
	include("startvals.jl")
	include("running.jl")
	include("plotter.jl")
	# include("interact.jl")





end # module
