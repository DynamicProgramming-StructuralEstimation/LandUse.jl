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
	using Formatting

	# constants
	const PEN = 100.0  # penalty for nl solver
	const dbpath = joinpath(ENV["HOME"],"Dropbox","research","LandUse","output","model")
	const dbplots = joinpath(dbpath,"plots")
	const dbtables = joinpath(dbpath,"tables")
	const originalθ = [0.32, 0.33, 0.34, 0.36, 0.38, 0.41, 0.48, 0.7, 1.35, 2.3, 3, 4.5, 5, 5.5]  # those number from initial matlab code

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
	include("experiments.jl")
	# include("interact.jl")





end # module
