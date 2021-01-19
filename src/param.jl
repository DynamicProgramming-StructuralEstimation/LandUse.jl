mutable struct Param
	γ     :: Float64 # housing weight
	ϵr     :: Float64 # housing supply elasticity in rural sector
	ϵs    :: Float64  # current slope of elasticity function
	ϵsmax :: Float64  # max slope of elasticity function
	nsteps :: Int  # num steps in thetau search
	η     :: Float64 # agglomeration forces
	ν     :: Float64 # weight of rural good consumption on consumption composite
	cbar  :: Float64 # agr cons subsistence level (0.14 for thailand)
	sbar  :: Float64 # neg urba subsistence level

	# productivity setup
	θu_g  :: Float64 # constant growth factors by sector
	θr_g  :: Float64 # constant growth factors by sector
	θut   :: Vector{Float64}  # idiosyncratic level of urban producitivity in each period
	θrt   :: Vector{Float64}  # idiosyncratic level of rural producitivity in each period
	θr    :: Float64 # current period rural sector TFP = θut[t]
	θu    :: Float64 # current period urban sector TFP

	# commuting cost setup
	ηm     :: Float64   # speed elasticity of fixed cost
	ηl     :: Float64   # location elasticity of fixed cost
	ηw     :: Float64   # wage elasticity of fixed cost
	cτ     :: Float64   # efficiency of transport technology
	ζ      :: Float64   # valuation of commuting time in terms of forgone wages
	a      :: Float64   # implied combination of above parameters
	tauw     :: Float64   # exponent on wu (equation 5 in commutingtech)
	taul     :: Float64   # exponent on l (equation 5 in commutingtech)

	α     :: Float64 # labor weight on farm sector production function
	λ     :: Float64 # useless land: non-farm, non-urban land (forests, national parks...)
	L     :: Float64 # total population
	Lt     :: Vector{Float64} # total population by year
	T     :: StepRange{Int64,Int64}
	t     :: Int64
	it     :: Int64
	σ     :: Float64 # land-labor elasticity of substitution in farm production function
	Ψ     :: Float64  # urban ammenities rel to rural
	int_nodes :: Int  # number of integration nodes
	iweights :: Vector{Float64}  # int weights
	inodes    :: Vector{Float64}  # points where to evaluate integrand (inodes scaled into [0,ϕ])
	S     :: Float64  # area of region
	ϕ1  :: Float64   # first period fringe
	ϕ1x  :: Float64   # fraction of first fringe that defines "center"

	# Country setup
	K :: Int # number of regions
	kshare :: Vector{Float64} # share of space of each region
	factors :: Vector{Float64} # growth factor offsets
	# kθu :: Dict  # collection of θu's for each region for each period
	# kθr :: Dict

	trace :: Bool  # whether to trace solver
	iters :: Int  # max iterations
	ma    :: Int64  # moving average window size
	mag    :: Float64  # assumed growth for extraploating producivities.

	moments :: DataFrame
	thetas :: DataFrame

	function Param(;par=Dict(),use_estimatedθ = false)
        f = open(joinpath(dirname(@__FILE__),"params.json"))
        j = JSON.parse(f)
        close(f)
        this = new()

		# read data from json file
        for (k,v) in j
			# if v["type"] == "region"
	            if v["value"] isa Vector{Any}
	                if k == "T"
	                	vv = v["value"]
		                setfield!(this,Symbol(k),vv[1]:vv[2]:vv[3])
		            else
		                setfield!(this,Symbol(k),convert(Vector{Float64},v["value"]))
		            end
	            else
	                setfield!(this,Symbol(k),v["value"])
	            end
			# elseif v["type"] == "numerical"
			# 	setfield!(this,Symbol(k),v["value"])
			# end
		end

		# set some defaults
		T = length(this.T)
		this.t = 1815
		this.it = 1
		this.S = 1.0  # set default values for space and population
		this.ϕ1 = NaN
		this.ϕ1x = 0.5


		# read data from disk
		# this.thetas = select(CSV.read(joinpath(LandUse.dbtables,"thetas_data.csv"), DataFrame), :year , :stheta_rural => :thetar, :stheta_urban => :thetau)
		moments = CSV.read(joinpath(LandUse.dbtables,"data-moments.csv"), DataFrame)
		this.thetas = select(moments, :year , :stheta_rural => :thetar, :stheta_urban => :thetau)

		if use_estimatedθ
			this.thetas = CSV.read(joinpath(dbtables,"export_theta_pr.csv"), DataFrame)
		end

		this.moments = moments[ moments.year .∈ Ref(this.T), : ]

		Lt = CSV.read(joinpath(LandUse.dbtables,"population.csv"), DataFrame)
		Lt = Lt[ Lt.year .∈ Ref(this.T), : ]
		this.Lt = Lt.population ./ Lt.population[1]
		# this.Lt = exp.(collect(range(log(1.0),log(2.42),length = T)))
		this.L = this.Lt[1]

		this.moments = innerjoin(this.moments, Lt, on = :year)

		# bring on scale we know that works
		# thetas.thetar .= thetas.thetar .* 0.32
		# thetas.thetau .= thetas.thetau .* 0.32
		# # pick out correct years
		tt = this.thetas[ this.thetas.year .∈ Ref(this.T),  : ]
		this.θut = tt.thetau
		this.θrt = tt.thetar

		# override parameters from dict par, if any
        if length(par) > 0
            for (k,v) in par
				if hasfield(Param,k)
					if (k == :θut)
						if length(v) > 1
							setfield!(this,k,v)
						else
							this.θut = [v ; Float64[growθ(v, [this.θu_g for i in 2:it]) for it in 2:T]]
						end
					elseif (k == :θrt)
						if length(v) > 1
							setfield!(this,k,v)
						else
							this.θrt = [v ; Float64[growθ(v, [this.θr_g for i in 2:it]) for it in 2:T]]
						end
					# elseif (k == :τ0t)
					# 	if length(v) > 1
					# 		setfield!(this,k,v)
					# 	else
					# 		this.τ0t = [v ; Float64[growθ(v, [this.τ0_g for i in 2:it]) for it in 2:T]]
					# 	end
					else
                		setfield!(this,k,v)
					end
				end
            end
			# if (:τ0t ∉ collect(keys(par)))
			# 	this.τ0t = [this.τ0t[1] ; Float64[growθ(this.τ0t[1], [this.τ0_g for i in 2:it]) for it in 2:T]]
			# end
        else
			# this.τ0t = [this.τ0t[1] ; Float64[growθ(this.τ0t[1], [this.τ0_g for i in 2:it]) for it in 2:T]]
		end

		this.θu = this.θut[1]
		this.θr = this.θrt[1]
		# this.τ  = this.τ0t[1]

		# set epsilon slope
		if this.ϵs != this.ϵsmax
			@warn "you set ϵs not equal to ϵsmax. I take ϵsmax => ϵs."
			this.ϵs = this.ϵsmax
		end

        if this.η != 0
        	@warn "current wage function hard coded \n to LU_CONST=$LU_CONST. Need to change for agglo effects!"
        end

		if this.ηm < 0 error("ηm < 0 violated") end
		if this.ζ > 1.0 || this.ζ < 0.0 error("ζ ∈ (0,1) violated") end
		if this.ηl > 1.0 || this.ηl < 0.0 error("ηl ∈ (0,1) violated") end
		# derived parameters
		this.a = this.cτ
		# this.a = ((1 + this.ηm) / this.ηm) * this.cτ^(1 / (1 + this.ηm)) * (2 * this.ζ)^((this.ηm + this.ηl) / (1 + this.ηm))
		# this.taum = 0.571
		this.tauw = (this.ηm + this.ηw) / (1+this.ηm)
		this.taul = (this.ηm + this.ηl) / (1+this.ηm)
		this.inodes, this.iweights = gausslegendre(this.int_nodes)

		if (isnan(this.ϕ1x) || (this.ϕ1x <= 0)) error("invalid value for ϕ1x: $(this.ϕ1x)") end

		# this.taul = 0.571
		# println("taum = $(this.taum)")
		# println("taul = $(this.taul)")
		# this.ew = (-1)/(1+this.ηm)
		# this.el = (this.ηm + this.ηl)/(1+this.ηm)
		# as ηm goes to infinity the transport cost goes to 2 ζ w l


    	return this
	end
end
function getgrowth(p::Param,s::Symbol,g::Float64)
	x = getfield(p,s)
	[x[1] ; Float64[growθ(x[1], [g for i in 2:it]) for it in 2:length(p.T)]]
end

function backoutη(p::Param)
	function eq!(F,x,p::Param)
		#x = [ηm, ηl]
		F[1] = p.taul - ((x[1] + x[2])/(1 + x[1]))
		F[2] = p.taum - ((x[1])/(1 + x[1]))
	end
	r = nlsolve((F,x) -> eq!(F,x,p),ones(2))
	if converged(r)
		return(Dict(:ηm => r.zero[1], :ηl => r.zero[2]))
	else
		error("not converged")
	end
end


"print default param to latex table"
function latex_param()
	f = open(joinpath(dirname(@__FILE__),"params.json"))
	j = JSON.parse(f)
	close(f)

	getline(x) = [latexstring(x["symbol"]), x["description"], x["value"]]

	latex_tabular(joinpath(dbtables,"params.tex"), Tabular("l l D{.}{.}{5.5}@{}"), [
	   Rule(:top),
       ["Parameter", "Description", MultiColumn(1,:c,"value")],
       Rule(:mid),
	   getline(j["S"]),
	   getline(j["L"]),
	   getline(j["γ"]),
	   getline(j["σ"]),
	   getline(j["α"]),
	   getline(j["ν"]),
	   getline(j["cbar"]),
	   getline(j["sbar"]),
	   getline(j["ηw"]),
	   getline(j["ηl"]),
	   getline(j["ηm"]),
	   getline(j["cτ"]),
	   getline(j["ζ"]),
	   getline(j["ϵr"]),
       Rule(:bottom)
	   ]
	)

end


function show(io::IO, ::MIME"text/plain", p::Param)
    print(io,"LandUse Param:\n")
	print(io,"      γ       : $(p.γ   )\n")
	print(io,"      ϵr      : $(p.ϵr  )\n")
	print(io,"      ϵs      : $(p.ϵs  )\n")
	print(io,"      ϵsmax   : $(p.ϵsmax  )\n")
	print(io,"      nsteps : $(p.nsteps  )\n")
	print(io,"      η       : $(p.η   )\n")
	print(io,"      ν       : $(p.ν   )\n")
	print(io,"      cbar    : $(p.cbar)\n")
	print(io,"      sbar    : $(p.sbar)\n")
	print(io,"      θr      : $(p.θr  )\n")
	print(io,"      θu      : $(p.θu  )\n")
	print(io,"      α       : $(p.α   )\n")
	print(io,"      λ       : $(p.λ   )\n")
	print(io,"      cτ       : $(p.cτ   )\n")
	print(io,"      L       : $(p.L   )\n")
	print(io,"      S       : $(p.S   )\n")
	print(io,"      T       : $(p.T   )\n")
	print(io,"      t       : $(p.t   )\n")
	print(io,"      σ       : $(p.σ   )\n")
	print(io,"      K       : $(p.K  )\n")
end


function smooth_p_check(periods)
	d = DataFrame(CSV.File(joinpath(LandUse.dbpath,"data","nico-output","FRA_model.csv")))
	d = d[.!ismissing.(d.P_rural),:]
	return d
	z = transform(d[!,[:year,:P_rural]], :P_rural => (y -> LandUse.smooth(collect(skipmissing(y)),periods)))
	@df z plot(:year, [:P_rural , :P_rural_function], seriestype = [:path],lw = 2,
			   title = "smoothing over $periods periods",
			   labels = ["" "Smoothed"])
end

"""
Takes raw csv data and prepares to be used in model.
writes to csv input files.
"""
function prepare_data(p::Param; digits = 9)

	# population data for each year (needs interpolation)
	pop = DataFrame(CSV.File(joinpath(LandUse.dbpath,"data","France-population.csv")))
	append!(pop, DataFrame(year = 2050, population = 74.0))  # https://www.insee.fr/fr/statistiques/2859843
	itp = interpolate((pop.year,), pop.population, Gridded(Linear()))
	popd = DataFrame(year = 1815:2050, population = itp(1815:2050))
	CSV.write(joinpath(dbtables,"population.csv"),popd)


	# productivity and employment moments

	dt = p.T
	ma = p.ma
	growth = p.mag

	d = DataFrame(CSV.File(joinpath(LandUse.dbpath,"data","nico-output","FRA_model.csv")))
	# @df d plot(:year, :P_rural)

	# plot(d.year, d.theta_rural)
	# plot(rand(10))
	x0 = @linq d |>
		where((:year .<= 2015) .& (:year .>= dt.start))   # rural data stops in 2015



	x = select(x0, :year, :theta_rural, :theta_urban,:Employment_rural, r"SpendingShare", :P_rural)


	# normalize year 1 to 1.0
	# x = transform(x, :theta_rural => (x -> x ./ x[1]) => :theta_rural, :theta_urban => (x -> x ./ x[1]) => :theta_urban)

	#impute
	Impute.interp!(x)

	# moving average smoother
	transform!(x, :theta_rural => (y -> smooth(collect(skipmissing(y)),ma)) => :stheta_rural,
	              :theta_urban => (y -> smooth(collect(skipmissing(y)),ma)) => :stheta_urban,
	              :P_rural => (y -> LandUse.smooth(collect(skipmissing(y)),61)) => :P_rural)

	if maximum(x.year) < dt.stop
		# append the future to end of data
		x1 = vcat([DataFrame(Dict(zip(names(x), [missing for i in 1:length(names(x))]))) for i in 1:length(2016:dt.stop)]...)
		x1.year = 2016:dt.stop
		allowmissing!(x)
		append!(x,x1)


		# interpolate forward
		Impute.locf!(x)
	end

	# from year 2000 onwards, replace rural with `growth` percent growth
	x[:, :stheta_rural] .= ifelse.(x.year .> 2000,vcat(zeros(sum(x.year .<= 2000)),[x[x.year .== 2000,:stheta_rural] * growth^i for i in 1:sum(x.year .> 2000)]...),x.stheta_rural)
    x[:, :stheta_urban] .= ifelse.(x.year .> 2000,vcat(zeros(sum(x.year .<= 2000)),[x[x.year .== 2000,:stheta_urban] * growth^i for i in 1:sum(x.year .> 2000)]...),
   							   x.stheta_urban)

    # round to 5 digits
	x[:, [:stheta_rural, :stheta_urban]] .= mapcols(col -> round.(collect(skipmissing(col)),digits = digits),x[:, [:stheta_rural, :stheta_urban]])

	x = transform(x, :theta_rural => (x -> x ./ x[1]) => :theta_rural, :theta_urban => (x -> x ./ x[1]) => :theta_urban)
	x = transform(x, :stheta_rural => (x -> x ./ x[1]) => :stheta_rural, :stheta_urban => (x -> x ./ x[1]) => :stheta_urban)


    p1 = @df x plot(:year, [:theta_rural :stheta_rural],leg=:left, lw = 2)
	savefig(p1, joinpath(dbplots,"smooth-theta-rural.pdf"))
	p2 = @df x plot(:year, [:theta_urban :stheta_urban],leg=:left, lw = 2)
	savefig(p2, joinpath(dbplots,"smooth-theta-urban.pdf"))

	p3 = @df x plot(:year, [:stheta_rural :stheta_urban],leg=:left, lw = 2)
	savefig(p3, joinpath(dbplots,"smooth-thetas-model.pdf"))

	p4 = @df x scatter(:year, :theta_rural, title = "Rural", label = "data",leg=:left)
	plot!(p4, x.year, x.stheta_rural, lw = 2, label = "model")

	p5 = @df x scatter(:year, :theta_urban, title = "Urban", label = "data",leg=:left)
	plot!(p5, x.year, x.stheta_urban, lw = 2, label = "model")

	p6 = plot(p4,p5,layout = (2,1),size = (800,500))
	savefig(p6, joinpath(dbplots,"smooth-thetas-data-model.pdf"))

	CSV.write(joinpath(dbtables,"data-moments.csv"),x)
	#
	#
	# (ret, pl)
	(p1,p2,p3,p6)
end

function plot_shares()
	d = DataFrame(CSV.File(joinpath(LandUse.dbpath,"data","nico-output","FRA_model.csv")))
	x = @linq d |>
		where((:year .< 2018) .& (:year .> 1895))

	x = select!(x,:year,:SpendingShare_Rural => :Rural, :SpendingShare_Urban => :Urban, :SpendingShare_Housing => :Housing)
	x = stack(x, Not(:year))

	pl = @df x plot(:year, :value, group = :variable,
	           linewidth = 2,ylab = "Spending Share", leg = :bottomright)

   savefig(pl, joinpath(dbdataplots,"spending-shares.pdf"))



end

# set period specific values
function setperiod!(p::Param,i::Int)
	setfield!(p, :θr, p.θrt[i])   # this will be constant across region.
	setfield!(p, :θu, p.θut[i])   # in a country setting, we construct the growth rate differently for each region.
	setfield!(p, :L, p.Lt[i])

	# setfield!(p, :θr, i == 1 ? p.θr0 : growθ(p.θr0,p.θrg[1:(i-1)]))   # this will be constant across region.
	# setfield!(p, :θu, i == 1 ? p.θu0 : growθ(p.θu0,p.θug[1:(i-1)]))   # in a country setting, we construct the growth rate differently for each region.
	setfield!(p,  :t , p.T[i] )
	setfield!(p,  :it , i )
end

growθ(θ0,g::Vector{Float64}) = θ0 * prod(g)






function setperiod!(p::Vector{Param},i::Int)
	for ip in eachindex(p)
		setperiod!(p[ip],i)
	end
end


function setfields!(p::Vector{Param},name::Symbol,x)
	for ip in eachindex(p)
		setfield!(p[ip],name,x)
	end
end

# ϵfun(d,s,ϕ,ϵtarget) = ϵtarget * exp(-s * max(ϕ-d,0.0))
function ϵfun_tmp(d,s,ϕ,p::Param)
	setfield!(p, :ϵs, s)
	ϵ(d,ϕ,p)
end




function plot_ϵfun(ϵsmax;ϕ = 0.7,ϵtarget = 4)
	p = Param(par = Dict(:ϵsmax => ϵsmax))
	xr = range(0.0,stop = 1.0,length=200)
	sr = range(0.0,p.ϵsmax,length = 10)
	# y = [(d,s) -> ϵfun(d,ϕ,p) for d in xr, s in sr]
	y = [ϵfun_tmp(d,s,ϕ,p) for d in xr, s in sr]

	p = plot(xr,y, labels=reshape(["s=$i" for i in sr],1,length(sr)),title = L"\epsilon(\phi) \exp(-s (\phi - d))", xlabel = "distance to center", ylabel = "Elasticity")
	vline!(p,[ϕ], color = :red, label = "")
	savefig(p, joinpath(@__FILE__,"..","..","images","epsfun.pdf"))
	p
end



"""
A Country-wide Parameter struct
"""
mutable struct CParam
	L     :: Float64 # total population
	Lt     :: Vector{Float64} # total population
	S     :: Float64 # total space
	K     :: Int  # number of Regions
	kshare    :: Vector{Float64}  # region k's share of total space


	function CParam(;par=Dict())
        f = open(joinpath(dirname(@__FILE__),"params.json"))
        j = JSON.parse(f)
        close(f)
        this = new()

        for (k,v) in j
			if v["type"] == "country"
	            if v["value"] isa Vector{Any}
	                setfield!(this,Symbol(k),convert(Vector{Float64},v["value"]))
	            else
	                setfield!(this,Symbol(k),v["value"])
	            end
			end
        end

        if length(par) > 0
            # override parameters from dict par
            for (k,v) in par
				if hasfield(CParam, k)
                	setfield!(this,k,v)
				end
            end
        end
        if length(this.kshare) != this.K
        	throw(ArgumentError("your settings for number of regions are inconsistent"))
        end
        if (this.K > 1) & !(sum(this.kshare) ≈ 1.0)
        	throw(ArgumentError("Shares of regions space must sum to 1.0"))
        end
    	return this
	end
end

function show(io::IO, ::MIME"text/plain", p::CParam)
    print(io,"LandUse Country Param:\n")
	print(io,"      L       : $(p.L   )\n")
	print(io,"      S       : $(p.S   )\n")
	print(io,"      K       : $(p.K   )\n")
	print(io,"      kshare  : $(p.kshare   )\n")
end

"convert a country param to one param for each region"
function convert(c::CParam; par = Dict())
	if length(par) > 0
		p = [Param(par = par[ik]) for ik in 1:c.K]
	else
		p = [Param(par = par) for ik in 1:c.K]
	end
	# do adjustment of parameters one by one
	for ik in 1:c.K
		# p[ik].L = c.kshare[ik] * c.L  # by default allocate this share of people to k
		# p[ik].S = c.kshare[ik] * c.S  # true by definition
		# these entries are only used to compute a one-region country (i.e. first step to get starting values)
		# we allocated total space and population to each country.
		p[ik].L = c.L  # by default allocate this share of people to k
		p[ik].S =  c.S  # true by definition

		# transform aggregate productivity series into region equivalentes

	end
	return p
end
