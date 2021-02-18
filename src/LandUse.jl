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
	using Printf
	using NLopt
	using LineSearches
	using Formatting
	using Roots
	using ColorSchemes
	using CSV
	using SmoothingSplines
	using LaTeXTabulars
	using QuantEcon: smooth
	using JuMP
	using Ipopt
	using Interpolations: interpolate, Gridded, Linear
	using Impute
	using Blink
	using BlackBoxOptim
	using Serialization
	using Dates
	using LsqFit

	# constants
	const PEN = 100.0  # penalty for nl solver
	user = splitdir(homedir())[end]
	isflo = (user == "florian.oswald") || (user == "74097")
	const dbpath = isflo ? joinpath(ENV["HOME"],"Dropbox","research","LandUse") : "/home/oswald/LandUseDropbox"
	const dbplots = joinpath(dbpath,"output","model","plots")
	const dbdataplots = joinpath(dbpath,"output","data","plots")
	const dbtables = joinpath(dbpath,"output","model","tables")
	const CTRY_MAXTRY = 100

	# imports
	import Base.show, Base.convert

	# globals
	C_TRIED = [0]

	# our code
	include("param.jl")
	include("model.jl")
	include("urban.jl")
	include("country.jl")
	include("startvals.jl")
	include("running.jl")
	include("plotter.jl")
	include("experiments.jl")
	include("interact.jl")
	include("jump.jl")
	include("estimation.jl")

	cpslides(name) = cp(joinpath(@__DIR__,"..","tex","slides","COT_slides.pdf"),
	                   joinpath(dbpath,"slides","flo-slides","COT_slides-$name.pdf"),
					   force = true)


end # module
