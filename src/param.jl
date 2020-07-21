mutable struct Param
	γ     :: Float64 # housing weight
	ϵr     :: Float64 # housing supply elasticity in rural sector
	ϵs    :: Float64  # current slope of elasticity function
	ϵsmax :: Float64  # max slope of elasticity function
	ϵnsteps :: Int  # num steps in elasticity search
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
	cτ     :: Float64   # efficiency of transport technology
	ζ      :: Float64   # valuation of commuting time in terms of forgone wages
	τ      :: Float64   # implied combination of above parameters
	ew     :: Float64   # exponent on wu (equation 5 in commutingtech)
	el     :: Float64   # exponent on l (equation 5 in commutingtech)


	α     :: Float64 # labor weight on farm sector production function
	λ     :: Float64 # useless land: non-farm, non-urban land (forests, national parks...)
	χr    :: Float64 # "easiness" to convert land into housing in rural sector
	# χu    :: Float64 # "easiness" to convert land into housing in center of urban sector
	L     :: Float64 # total population
	T     :: StepRange{Int64,Int64}
	t     :: Int64
	σ     :: Float64 # land-labor elasticity of substitution in farm production function
	c0    :: Float64  # construction cost function intercept
	c1    :: Float64  # construction cost function gradient
	c2    :: Float64  # construction cost function quadratic term
	Ψ     :: Float64  # urban ammenities rel to rural
	int_nodes :: Int  # number of integration nodes
	S     :: Float64  # area of region
	ρrbar :: Float64  # fixed rural land value for urban model

	trace :: Bool  # whether to trace solver

	function Param(;par=Dict())
        f = open(joinpath(dirname(@__FILE__),"params.json"))
        j = JSON.parse(f)
        close(f)
        this = new()

		# read data from json file
        for (k,v) in j
			if v["type"] == "region"
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
			elseif v["type"] == "numerical"
				setfield!(this,Symbol(k),v["value"])
			end
        end

		# set some defaults
		T = length(this.T)
		this.t = 1
		this.S = 1.0  # set default values for space and population
		this.L = 1.0
		thetas = smooth_θ(this.T)
		this.θut = thetas[:θu]
		this.θrt = thetas[:θr]

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
		this.τ = ((1 + this.ηm) / this.ηm) * this.cτ^(1 / (1 + this.ηm)) * (2 * this.ζ)^((this.ηm + this.ηl) / (1 + this.ηm))
		this.ew = (-1)/(1+this.ηm)
		this.el = (this.ηm + this.ηl)/(1+this.ηm)


    	return this
	end
end
function getgrowth(p::Param,s::Symbol,g::Float64)
	x = getfield(p,s)
	[x[1] ; Float64[growθ(x[1], [g for i in 2:it]) for it in 2:length(p.T)]]
end


"print default param to latex table"
function latex_param()
	f = open(joinpath(dirname(@__FILE__),"params.json"))
	j = JSON.parse(f)
	close(f)

	getline(x) = [x["symbol"], x["symbol"], x["value"]]

	latex_tabular(joinpath(dbtables,"params.tex"), Tabular("l l D{.}{.}{5.5}@{}"), [
	   Rule(:top),
       ["name", "description", MultiColumn(1,:c,"value")],
       Rule(:mid),
	   getline(j["L"]),
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
	print(io,"      ϵnsteps : $(p.ϵnsteps  )\n")
	print(io,"      η       : $(p.η   )\n")
	print(io,"      ν       : $(p.ν   )\n")
	print(io,"      cbar    : $(p.cbar)\n")
	print(io,"      sbar    : $(p.sbar)\n")
	print(io,"      θr      : $(p.θr  )\n")
	print(io,"      θu      : $(p.θu  )\n")
	print(io,"      α       : $(p.α   )\n")
	print(io,"      λ       : $(p.λ   )\n")
	print(io,"      τ       : $(p.τ   )\n")
	print(io,"      χr      : $(p.χr  )\n")
	print(io,"      L       : $(p.L   )\n")
	print(io,"      S       : $(p.S   )\n")
	print(io,"      T       : $(p.T   )\n")
	print(io,"      t       : $(p.t   )\n")
	print(io,"      σ       : $(p.σ   )\n")
	print(io,"      c0      : $(p.c0  )\n")
	print(io,"      c1      : $(p.c1  )\n")
	print(io,"      c2      : $(p.c2  )\n")
	print(io,"      Ψ       : $(p.Ψ   )\n")
end

function smooth_θ(dt::StepRange)
	d = DataFrame(CSV.File(joinpath(LandUse.dbpath,"data","nico-output","FRA_model.csv")))
	x = @linq d |>
		where((:year .<= dt.stop) .& (:year .>= dt.start))

	# normalize to year 1
	r = select(x,:year,:theta_rural => (x -> x ./ x[1]) => :theta)
	u = select(x,:year,:theta_urban => (x -> x ./ x[1]) => :theta)

	r = r[completecases(r),:]
	disallowmissing!(r)
	# from year 2000 onwards, first do a moving average smoothing
	iy = findfirst(dt .>= 2000)
	ir = r.theta[iy:end]
	ma(v,n) = [mean(v[i:(i+n-1)]) for i in 1:(length(v)-(n-1))]
	r.theta[iy:(end-3)] = ma(ir,4)
	r.theta[(end-3):end] .= r.theta[(end-3)]
	sr = SmoothingSplines.fit(SmoothingSpline, convert(Array{Float64},r.year), r.theta, 250.0)

	u = u[completecases(u),:]
	disallowmissing!(u)
	su = SmoothingSplines.fit(SmoothingSpline, convert(Array{Float64},u.year), u.theta, 250.0)

	pu = SmoothingSplines.predict(su,convert(Array{Float64},dt))
	pr = SmoothingSplines.predict(sr,convert(Array{Float64},dt))

	# manually set pr to flat line after last data point
	pr[dt .>= 2015] .= pr[findfirst(dt .>= 2015)]



	ret = Dict(:θr => pr, :θu => pu )
	plu = scatter(u.year, u.theta ./ u.theta[1],title = "Urban Productivity",leg=false, ylabel = "$(dt.start) = 1")
	plot!(plu,dt,ret[:θu],m=(3,:auto,:red))
	savefig(plu, joinpath(dbplots,"smooth-thetau.pdf"))

	plr = scatter(r.year, r.theta ./ r.theta[1],title = "Rural Productivity",leg=false, ylabel = "$(dt.start) = 1")
	plot!(plr,dt,ret[:θr],m=(3,:auto,:red))
	savefig(plr, joinpath(dbplots,"smooth-thetar.pdf"))

	pls = plot(plu,plr,layout = (1,2),link = :y)
	savefig(pls, joinpath(dbplots,"smooth-thetas.pdf"))


	ret[:θu] = ret[:θu] ./ ret[:θu][1]
	ret[:θr] = ret[:θr] ./ ret[:θr][1]

	ret
end



# function theta_rec(x::Float64,cp::CParam,i::Int)
# 	if i==1
# 		cp.θg[i] * x
# 	else
# 		theta_rec(x,cp,i-1)
# 	end
# end

# set period specific values
function setperiod!(p::Param,i::Int)
	setfield!(p, :θr, p.θrt[i])   # this will be constant across region.
	setfield!(p, :θu, p.θut[i])   # in a country setting, we construct the growth rate differently for each region.
	# setfield!(p, :τ, p.τ0t[i])

	# setfield!(p, :θr, i == 1 ? p.θr0 : growθ(p.θr0,p.θrg[1:(i-1)]))   # this will be constant across region.
	# setfield!(p, :θu, i == 1 ? p.θu0 : growθ(p.θu0,p.θug[1:(i-1)]))   # in a country setting, we construct the growth rate differently for each region.
	setfield!(p,  :t , p.T[i] )
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
