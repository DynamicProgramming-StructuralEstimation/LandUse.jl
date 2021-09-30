module LandUse

	# dependendencies
	using Base: root_module
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
	using LineSearches
	using Formatting
	using Roots
	using ColorSchemes
	using CSV
	using SmoothingSplines
	using LaTeXTabulars
	using QuantEcon: smooth, gridmake
	using JuMP
	using Ipopt
	using Interpolations: interpolate, Gridded, Linear
	using Impute
	using Blink
	using BlackBoxOptim
	using Serialization
	using Dates
	using LsqFit
	using SharedArrays
	using ProgressMeter
	using Distributed
	using DelimitedFiles
	using BSON
	using Flux
	using GLM
	using FileIO
	using StatsBase
	using Statistics

	# constants
	const PEN = 100.0  # penalty for nl solver

	# set up file paths
	isCI = if haskey(ENV,"CI")
			true 
		else
			false
		end

	# root dir
	root = joinpath(@__DIR__,"..")
	indir = joinpath(root,"in")
	intables = joinpath(indir,"tables")

	user = splitdir(homedir())[end]
	host = gethostname()
	isflo = (user == "florian.oswald") || (user == "74097")
	dbpath = if isflo
					joinpath(ENV["HOME"],"Dropbox","research","LandUse")
				elseif (contains(host,"cnode") || contains(host,"malbec"))
					"/home/oswald/LandUseDropbox"
				elseif host == "scpo-floswald"
					"/home/floswald/LandUseDropbox"	
				elseif isCI
					joinpath(@__DIR__,"..","in")
				end				
	dbplots = joinpath(dbpath,"output","model","plots")
	dboutdata = joinpath(dbpath,"output","data")
	dbdataplots = joinpath(dbpath,"output","data","plots")
	dbtables = joinpath(dbpath,"output","model","tables")
	CTRY_MAXTRY = 100

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
	include("learning.jl")
	include("data.jl")

	cpslides(name) = cp(joinpath(root,"tex","slides","COT_slides.pdf"),
	                   joinpath(dbpath,"slides","flo-slides","COT_slides-$name.pdf"),
					   force = true)


end # module
