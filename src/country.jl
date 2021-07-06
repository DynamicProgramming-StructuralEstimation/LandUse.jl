

"""
	Country

A `Country` is a collection of single `Regions`. There is free migration
across regions (same utilty everywhere). What characterises a Region
are different θs.

## shared characteristics

* total population
* rural wage
* rent and food price
"""
mutable struct Country
	R  :: Vector{Region}   # a set of regions
	pp  :: Vector{Param}   # a set of Params
	K  :: Int              # number of regions
	wr :: Float64          # a global rural wage
	pr :: Float64          # a global price of rural good
	ρr :: Float64          # a global land price
	r  :: Float64          # country-wide per capita rental income
	LS  :: Float64          # Lr/Sr for region 1
	Sk :: Vector{Float64}  # total space for each region
	L  :: Float64  # total population
	S  :: Float64 # total space
	T  :: StepRange{Int,Int}

	"""
		Country constructor

	want to give single p and get back
	differnet theta series for each country
	"""
	function Country(p::Param;istest = false)
		this = new()
		@assert p.K == length(p.factors)
		@assert p.K == length(p.gs)
		@assert p.K == length(p.kshare)
		this.K = p.K

		# create K copies of parameter
		this.pp = Param[deepcopy(p) for _ in 1:p.K]

		# modify θus for each
		@warn "modifying θu data in periods 2 and 3 to be == 1.0" maxlog=1

		for ik in 1:p.K
			if !istest 
				p.θu = max(p.θu, 1.0) 
			end
			this.pp[ik].θu =  p.θu * p.factors[ik] * exp(p.gs[ik] * (p.it-1))
		end

		this.R = [Region(this.pp[ik]) for ik in 1:p.K]


		this.pr = NaN
		this.ρr = NaN
		this.r  = NaN
		this.LS  = NaN
		this.L  = p.Lt[p.it] * p.K
		this.S  = p.S * p.K # total space
		this.Sk  = p.kshare .* this.S
		this.T = p.T
		return this
	end
end

country(; par = Dict(:K => 2, :kshare => [0.5,0.5], :factors => [1.0,1.0])) = Country(Param(par = par))

function show(io::IO, C::Country)
    # print(io,"Region: ϕ=$(round(m.ϕ,digits=3)), pop=$(pop(m)), area=$(round(area(m),digits=2))")
	for ik in 1:C.K
		m = C.R[ik]
		show(io,m)
    	# @printf(io,"Region %d: θu=%1.3f, θr=%1.3f, ϕ=%1.4f, area=%1.2f, Lu=%1.3f, Lr=%1.3f, pop=%1.3f, pr=%1.3f",ik,m.θu, m.θr, m.ϕ, area(m), m.Lu, m.Lr,pop(m))
    	if ik < C.K
			@printf(io,"\n ")
		end
	end
end

"Country Rural Market Clearing"
Rmk(C::Country) = sum(C.R[ik].icr + C.R[ik].Lr * cr(C.R[ik].ϕ,C.pp[ik],C.R[ik]) for ik in 1:C.K) -
                                   sum(Yr(C.R[ik],C.pp[ik]) for ik in 1:C.K)


"Obtain a Time Series from an array of Country as a DataFrame"
function dataframe(C::Vector{Country})
	K = length(C[1].R)
	# cols = setdiff(fieldnames(LandUse.Region),(:cr01,:cu01,:inodes,:iweights,:nodes))
	# region 1
	tt = C[1].T
	ir = 1
	p = C[1].pp[1]
	df = dataframe([C[it].R[1] for it in 1:length(tt)],p)
	# df = DataFrame(year = C[1].T, region = ir)
	df.region = [ir for i in 1:length(tt)]
	# for fi in cols
	# 	df[!,fi] = [getfield(C[it].R[ir],fi) for it in 1:length(C[1].T)]
	# end
	# other regions
	if K > 1
		for ir in 2:K
			df2 = dataframe([C[it].R[ir] for it in 1:length(tt)],p)
			df2.region = [ir for i in 1:length(tt)]


			# df2 = DataFrame(year = C[1].T, region = ir)
			# for fi in cols
			# 	df2[!,fi] = [getfield(C[it].R[ir],fi) for it in 1:length(C[1].T)]
			# end
			append!(df,df2)
		end
	end
	df
end


"""
	update!(c::Country,p::Vector{Param},x::Vector{Float64})

Update a `Country` object with a current vector of choice variables `x` supplied by the solver.
The ordering in `x` is:

1. b: ratio of labor to land in region 1
2. r: land rent
3. pr: relative price rural good
4. 4 : (K+3), Sr: amount of land use in rural production (complement of Srh)
5. (K+4) : (2K+3), Lu: urban pop in each k
5. (2K+4) : end, θu: urban prod in each k
"""
function update!(c::Country,x::Vector{Float64};estimateθ=false)
	K = c.K
	p = c.pp

	# update country wide components
	c.LS   = x[1]   # constant labor/land share in region 1. i called that b>0 in the doc.
	c.r    = x[2]   # land rent
	c.pr   = x[3]   # relative price rural good
	Srk    = x[4:(K+3)]  # Sr for each region k
	Lu     = x[(K+4):(2K+3)]  # Lu for each region k
	ϕs     = x[(2K+4):(3K+3)]  # Lu for each region k
	if estimateθ
		θus    = x[(3K+3+1):end]  # θu for each region k
	end


	@assert K == length(c.R)

	p1 = p[1]  # get first region's parameter to save typing

	# global implications
	σ1 = (p1.σ - 1) / p1.σ
	σ2 = 1 / (p1.σ - 1)
	# rural land price
	c.ρr = foc_Sr( c.LS , c.pr, p1)
	# rural wage
	c.wr = foc_Lr( c.LS , c.pr, p1)

	# update each region
	# 2. update other equations in each region
	for ik in 1:K
		cx = [c.ρr, 
		      ϕs[ik], 
			  c.r, 
			  c.LS * Srk[ik], 
			  c.pr, 
			  Srk[ik],
			  estimateθ ? θus[ik] : p[ik].θu,
			  p[ik].θr
			  ]
		update!(c.R[ik],p[ik],cx, Lu = Lu[ik])
	end
end


"""
Compute exponential decay model of density for each region in a country
"""
function expmodel(C::Country)
	p = C.pp[1]  # get country 1's param vector
	K = p.K

	coefs = Vector{Float64}[]
	bins  = Vector{Float64}[]
	# bins  = StepRangeLen[]
	densities  = Vector{Float64}[]

	# country 1 is assumed largest
	for ik in 1:K
		# need ratio of radii between each region
		# ϕratio = C.R[ik].ϕ / C.R[1].ϕ
		# ikbins = range(1.0 * ϕratio, p.int_bins * ϕratio, length = p.int_bins)
		# ikbins = range(1.0 , p.int_bins, length = p.int_bins)
		ikbins = C.R[ik].ϕmids
		ndensities = C.R[ik].iDensities ./ C.R[ik].iDensities[1]
		gradient,emod = expmodel(ikbins, ndensities)
		push!(coefs, gradient)
		push!(densities, ndensities)
		push!(bins, ikbins)
	end
	(coefs, densities, bins)
end




##########
# archive



"""
Computes the entries of the residual vector ``u``
"""
function EqSys!(F::Vector{Float64},C::Country)

	K = C.K
	p = C.pp

	fi = 1  # running index
	# 1. agg labor market clearing
	F[fi] = C.L - sum(pop(i) for i in C.R)

	# 2. land market clearing in each region K
	for ik in 1:K
		fi += 1
		F[fi] = C.Sk[ik] - C.R[ik].Sr - C.R[ik].ϕ^2 * π - C.R[ik].Srh
	end

	# # 3. rural labor/land ratio is constant across K (and equal to region 1)
	# for ik in 2:K
	# 	fi += 1
	# 	F[fi] = (C.R[ik].Lr / C.R[ik].Sr) - (C.R[1].Lr / C.R[1].Sr)
	# end

	# 3. Aggregate land rents
	fi += 1
	# F[fi] = C.r * p[1].L - sum(m.iq + m.ρr * (m.Sk - m.ϕ) for m in C.R)
	F[fi] = C.r * C.L - sum(C.R[ik].iq + C.R[ik].ρr * (C.Sk[ik] - C.R[ik].ϕ) for ik in 1:K)

	# 4. aggregate urban consumption good clears
	urban_prod = sum(Yu(C.R[ik],p[ik]) for ik in 1:K)
	fi += 1
	m = C.R
	F[fi] = urban_prod -
	         sum( m[ik].Lr * cur(p[ik],m[ik]) +
			      m[ik].icu +
				  m[ik].Srh * cu_input(m[ik].ϕ,p[ik],m[ik]) +
				  m[ik].icu_input +
				  m[ik].wu0 * m[ik].iτ for ik in 1:K)

	# K + 3 equations up to here

	# city size - Urban population relationship: equation (19)
	for ik in 1:K
		fi += 1
		F[fi] = m[ik].Lu - m[ik].iDensity
	end
	# push!(Ftrace,copy(F))

end

function solve!(F,x,C::Country)
	# push!(Xtrace,copy(x))
	# println(x)
	if any(x .< 0)
		F[:] .= 10.0 .+ x.^2
	else
		update!(C,x)
		EqSys!(F,C)
	end
end

function runk_archive(;par = Dict(:K => 2, :kshare => [0.5,0.5], :factors => [1.0,1.0]))

	C = LandUse.Country[]
   	sols = Vector{Float64}[]  # an empty array of vectors

	pp = LandUse.Param(par = par, use_estimatedθ = true)
	@assert pp.K > 1
	K = pp.K
	# pp = convert(cp,par = par)  # create Lk and Sk for each region. if par is not zero, then first level index (Int) is region.

	# 1. run a single region model
	x0,Mk,p0 = LandUse.runm()

	# 2. all regions in one country now. starting values for Sr from Mk.
	push!(C,LandUse.Country(pp))  # create that country

	# 3. make sure all countries start at θu[1]
	for ik in 2:K
		C[1].pp[ik].θu = C[1].pp[1].θu
	end

	# starting values.
	# 1. b: ratio of labor to land in region 1
	# 2. r: land rent
	# 3. pr: relative price rural good
	# 4. 4 - K+3, Sr: amount of land use in rural production (complement of Srh)
	# x0 = country_starts(Mk[1],K)
	x0 = zeros(2K + 3)
	x0[1] = Mk[1].Lr / Mk[1].Sr
	x0[2] = Mk[1].r
	x0[3] = Mk[1].pr
	for ik in 1:K
		x0[3 + ik] = Mk[ik].Sr
		x0[K + 3 + ik] = Mk[1].Lu
	end

	# scale_factors = ones(length(pp[1].T))
	# scale_factors[7:end] .= 0.5
	r = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,C[1]),x0,iterations = 100, show_trace=false)

	if converged(r)
		push!(sols, r.zero)
		update!(C[1],r.zero)
	else
		println(r)
		error("Country not converged")
	end
	for it in 2:5
	# for it in 2:length(pp.T)
		println(it)
		# reset tried global
		C_TRIED = [0]
		# display(hcat(sols...)')
		setperiod!(pp, it)   # set param to it - in particular θs
		# println(pp.θu)
		# if it > 12
		# 	C0,x = adapt_θ(cp,pp,sols[it-1],it,do_traceplot=true)
		# else
			C0,x = step_country(sols[it-1],pp,it,do_traceplot=pp.trace)
		# end
       	push!(sols,x)
       	push!(C,C0)
   end
   sols, C
end

function step_country(x0::Vector{Float64},pp::Param,it::Int; do_traceplot = true)
	C0 = LandUse.Country(pp)
	# println("period $it")
	# println("params $(C0.pp[1].θu)")
	# println("params $(C0.pp[2].θu)")
	r = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,C0),x0,iterations = 100, store_trace=do_traceplot, extended_trace = do_traceplot)

	# r = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,pp,C0),x0,iterations = 1000,
	# 						show_trace=false,
	# 						extended_trace=true,
	# 						factor = nlsolve_fac,
	# 						autoscale = true,ftol = 0.000004)
							# method = :newton,
							# linesearch = LineSearches.BackTracking(order=3))

	if converged(r)
		# push!(sols, r1.zero)
		update!(C0,r.zero)
		if do_traceplot
			countrytraceplot(r,it)
		end
		# reset traces
		# global Xtrace = Vector{Float64}[]
		# global Ftrace = Vector{Float64}[]
		return C0, r.zero
	else
		# if C_TRIED[1] < CTRY_MAXTRY
		# 	x0 = x0 .+ (randn(length(x0)) .* 0.01 .* x0)
		# 	println("starting at $x0")
		# 	C_TRIED[1] = C_TRIED[1] + 1
		# 	step_country(x0,pp,it, do_traceplot = do_traceplot)
		# 	if do_traceplot
		# 		countrytraceplot(r,it)
		# 	end
		# else
			println(r)
			error("Country not converged in period $it after $CTRY_MAXTRY trials")
		# end
	end
end

# for each period, start all countries at eps = 0 and
# step wise increase the slope until maxeps
function adapt_ϵ(cp::CParam,p::Vector{Param},x0::Vector{Float64},it::Int; do_traceplot = false)
	sols = Vector{Float64}[]  # an empty array of vectors
	C = Country[]  # an empty array of Countries
	push!(sols, x0)  # put solution previous period

	# range of elasticity slope values
	ϵs = range(0,stop = p[1].ϵsmax, length = p[1].ϵnsteps)[2:end]

	for (i,ϵ) in enumerate(ϵs)
		# @debug "adapting to ϵ=$ϵ in $it"
		setfields!(p, :ϵs, ϵ)  # set current value for elaticity function slope on all params
		c,x = step_country(sols[i],p,cp,it,do_traceplot=do_traceplot)
		push!(sols,x)
		push!(C,c)
	end
	return C[end], sols[end]
end

function adapt_θ(C::Country,x0::Vector{Float64},it::Int; do_traceplot = true,maxstep = 0.1)

	# how many steps to take by looking at coutnry 1
	step = C.pp[1].Θu[it] - C[1].pp.Θu[it-1]
	s = Int(fld(step,maxstep)) + 1

	# range of values to achieve next step
	θs = [range(C[k].pp.θu - step, stop = C[k].pp.θu, length = s+1)[2:end]  for k in 1:C.K] # first one is done already
	sols = Vector{Float64}[]
	cs = Country[]
	push!(sols,x0)

	for (i,θ) in enumerate(θs)
		@debug "θ adapting" step=i θ=θ
		LandUse.setfields!(p, :θu, θ)   # sets on each member of p
		LandUse.setfields!(p, :θr, θ)
		# C0,x = adapt_ϵ(cp,p,sols[i],it,do_traceplot=true)
		C0,x = step_country(sols[i],p,cp,it,do_traceplot=do_traceplot)
		push!(sols,x)
		push!(cs,C0)
   	end

	return (cs[end],sols[end])    # return final step

end

function NLopt_wrap(result::Vector, x::Vector, grad::Matrix,C::Country,p::Vector{Param})
	if length(grad) > 0
		# not implemented
	end
	solve!(result,x,p,C)
end

function nlopt_solveC(p::Vector{Param},C::Country,x0::Vector{Float64})
	K = length(C.R)
	nx = 3 + K
	opt = Opt(:LN_COBYLA,nx)
	# opt = Opt(:LN_BOBYQA,nx)
	# opt = Opt(:AUGLAG_EQ,nx)
	# opt = Opt(:GN_ISRES,nx)
	opt.lower_bounds = fill(0.001,nx)
	# opt.upper_bounds = fill(100.0,nx)
	# opt.upper_bounds = [1.0,1.0,1.0,1.0]
	f0(x::Vector, grad::Vector) = 1.0
	opt.min_objective = f0  # fix at a constant function
	equality_constraint!(opt,(r,x,g) -> NLopt_wrap(r,x,g,C,p), fill(1e-9,nx))
	# m.r    = x[1]   # land rent
	# m.Lr   = x[2]   # employment in rural sector
	# m.pr   = x[3]   # relative price rural good
	# m.Sr   = x[4]   # amount of land used in rural production
	# (optf,optx,ret) = optimize(opt, [0.1 * p.S, 0.5 * p.L, 0.3775, 0.545 * p.S])
	# if isnothing(x0)
	# 	x0 = [0.1 * p.S, 0.5 * p.L, 0.3775, 0.545 * p.S]   # an almost arbitrary point in the interior domain.
	# else
		@assert length(x0) == ndims(opt)
	# end
	(optf,optx,ret) = optimize(opt, x0)
end


function country_starts(M,K)
       x0 = zeros(K + 3)
       x0[1] = M.Lr / M.Sr
       x0[2] = M.r
       x0[3] = M.pr
       for ik in 1:K
               x0[3 + ik] = M.Sr
       end
       x0
end
