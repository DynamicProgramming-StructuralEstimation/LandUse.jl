
abstract type Model end


"""
Region type

Urban expansion model with flexible commuting cost and differential supply elasticity.
This represents a single region. A region contains one potential urban area, and a rural area.
A `Region` can be part of a wider `Country`.
"""
mutable struct Region <: Model
	ρr   :: Float64   # land price in rural sector
	qr   :: Float64     # housing price in rural sector
	ρ0   :: Float64     # land price in center of city
	ρ0_y :: Float64     # land price in center of city over income
	ρr_y :: Float64     # land price in center of city over income
	q0   :: Float64     # housing price in center of city
	qq1   :: Float64     # housing price first quintile
	qbar :: Float64     # average housing price in city
	Lr   :: Float64   # employment in rural sector
	Lu   :: Float64   # employment in urban sector
	wu0  :: Float64   # wage in urban sector (at center)
	wr   :: Float64   # wage in rural sector
	Sr   :: Float64   # Amount of land used in rural production
	Srh  :: Float64   # Amount of land used for rural housing
	r    :: Float64   # per capita land rental income
	pr   :: Float64   # relative price of rural good
	H0   :: Float64   # housing supply at center
	Hr   :: Float64   # housing supply at fringe
	h0   :: Float64   # housing demand at center
	hr   :: Float64   # housing demand at fringe
	dbar :: Float64   # total avg density
	d0   :: Float64   # density in central area from 0 to first period fringe
	d0_   :: Float64   # density at point zero
	dr   :: Float64   # density at fringe
	cr0  :: Float64   # cons of rural in center
	cr1  :: Float64   # cons of rural at fringe
	cu0  :: Float64   # cons of urban in center
	cu1  :: Float64   # cons of urban at fringe
	ϕ    :: Float64   # radius of the city
	ϕ10    :: Float64   # 0.1 of radius
	ϕ90    :: Float64   # 0.9 of radius
	cityarea    :: Float64   # area of the city
	rel_cityarea    :: Float64   # area of the city relative to rural area
	citydensity    :: Float64   # average density of the city
	ϕbreaks :: StepRangeLen  # bins of distance
	speedbreaks :: Vector{Float64}  # breaks of speed
	ϕmids :: Vector{Float64}  # 
	xsr  :: Float64  # excess subsistence rural worker
	U    :: Float64  # common utility level
	θu   :: Float64  # current productivity
	θr   :: Float64  # current productivity
	pcy  :: Float64  # per capita income
	GDP  :: Float64  # GDP
	y    :: Float64  # disposable income
	Yu    :: Float64  # urban production
	Yr    :: Float64  # rural production
	mode0 :: Float64  # mode at center
	modeϕ :: Float64  # mode at fringe
	ctime0 :: Float64  # commute time at center
	ctimeϕ :: Float64  # commute time at fringe

	Cu :: Float64  # aggregate urban consumption
	Cr :: Float64  # aggregate rural consumption
	Ch :: Float64  # aggregate housing consumption
	C  :: Float64  # aggregate total consumption
	# integration setup
	inodes   :: Vector{Float64}  # theoretical integration nodes
	iweights :: Vector{Float64}  # int weights
	nodes    :: Vector{Float64}  # points where to evaluate integrand (inodes scaled into [0,ϕ])
	nodes_center    :: Vector{Float64}  # points where to evaluate integrand (inodes scaled into [0,ϕ1])
	nodes_10    :: Vector{Float64}  
	nodes_90    :: Vector{Float64}  
	nodes_bins  :: Matrix{Float64}  # (nbins, int_nodes)
	nodes_speeds  :: Matrix{Float64}

	# resulting integrals from [0,ϕ]
	icu_input :: Float64   # ∫ cu_input(l) dl
	iDensity  :: Float64   # ∫ D(l) dl
	iDensity_center  :: Float64   # ∫ D(l) dl for l ∈ [0,ϕ1]
	icu       :: Float64   # ∫ cu(l) D(l) dl
	icr       :: Float64   # ∫ cr(l) D(l) dl
	iτ        :: Float64   # ∫ τ(l) D(l) dl
	iq        :: Float64   # ∫ ρ(l) dl
	iy        :: Float64   # ∫ w(l) D(l) dl
	ihexp       :: Float64   # ∫ q(l) h(l) D(l) dl
	imode       :: Float64   # ∫ mode(l) D(l) dl
	ictime       :: Float64   # ∫ ctime(l) D(l) dl
	iDensity_q10  :: Float64   # ∫ D(l) 2π dl for l ∈ [0,q10]
	iDensity_q90  :: Float64   # ∫ D(l) 2π dl for l ∈ [q10,q90]
	iDensities  :: Vector{Float64}   # ∫ D(l) 2π dl for l ∈ bins(ϕ)
	iDensitySpeeds  :: Vector{Float64}   # ∫ D(l) 2π dl for l ∈ {(0,walk),(walk,bus), (bus,ϕ)}

	function Region(p::Param)
		# creates a model fill with NaN
		m      = new()
		m.ρr   = NaN
		m.qr   = NaN
		m.ρ0   = NaN
		m.ρ0_y   = NaN
		m.q0   = NaN
		m.qq1   = NaN
		m.qbar = NaN
		m.Lr   = NaN
		m.Lu   = NaN
		m.wu0  = NaN
		m.wr   = NaN
		m.Sr   = NaN
		m.Srh  = NaN
		m.r    = NaN
		m.pr   = NaN
		m.Hr   = NaN
		m.H0   = NaN
		m.hr          = NaN
		m.h0          = NaN
		m.dr          = NaN
		m.dbar        = NaN
		m.d0_         = NaN
		m.d0          = NaN
		m.ϕ           = NaN
		m.ϕ10           = NaN
		m.ϕ90           = NaN
		m.cityarea    = NaN
		m.rel_cityarea    = NaN
		m.citydensity = NaN
		m.ϕbreaks       = range(0,stop = 1, length = p.int_bins+1)
		m.speedbreaks   = zeros(p.nspeeds + 1)  # [0, lofmode(0.25), lofmode(0.41), ϕ]
		m.ϕmids         = zeros(p.int_bins)
		m.xsr         = NaN
		m.θu          = p.θu
		m.θr          = p.θr
		m.U           = NaN
		m.pcy         = NaN
		m.GDP         = NaN
		m.y           = NaN
		m.inodes = p.inodes
		m.iweights = p.iweights
		m.nodes = zeros(p.int_nodes)
		m.nodes_center = zeros(p.int_nodes)
		m.nodes_10 = zeros(p.int_nodes)
		m.nodes_90 = zeros(p.int_nodes)
		m.nodes_90 = zeros(p.int_nodes)
		m.nodes_bins = zeros(p.int_bins,p.int_nodes)
		m.nodes_speeds = zeros(p.nspeeds,p.int_nodes)
		m.icu_input = NaN
		m.iDensity  = NaN
		m.iDensity_center = NaN
		m.iDensity_q10    = NaN
		m.iDensity_q90    = NaN
		m.iDensities      = fill(NaN,p.int_bins)
		m.iDensitySpeeds  = zeros(p.nspeeds)
		m.icu       = NaN
		m.iτ        = NaN
		m.iq        = NaN
		m.iy        = NaN
		m.ihexp       = NaN
		m.imode       = NaN
		m.mode0       = NaN
		m.modeϕ       = NaN
		m.ictime       = NaN
		m.ctime0       = NaN
		m.ctimeϕ       = NaN

		return m
	end
end



function show(io::IO, ::MIME"text/plain", m::Model)
    print(io,"Single Region:\n")
    print(io,"    pop   : $(pop(m)) \n")
    print(io,"    Lr   : $(m.Lr ) \n")
    print(io,"    Lu   : $(m.Lu ) \n")
    print(io,"    area : $(round(area(m),digits=2)) \n")
    print(io,"    θu    : $(m.θu  ) \n")
    print(io,"    θr    : $(m.θr  ) \n")
    print(io,"    ϕ    : $(m.ϕ  ) \n")
    print(io,"    cityarea    : $(m.cityarea  ) \n")
    print(io,"    rel_cityarea    : $(rel_cityarea(m))\n")
    print(io,"    Sr   : $(m.Sr ) \n")
    print(io,"    Srh  : $(m.Srh) \n")
    print(io,"    ρr   : $(m.ρr ) \n")
    print(io,"    qr   : $(m.qr ) \n")
    print(io,"    wu0  : $(m.wu0) \n")
    print(io,"    wr   : $(m.wr ) \n")
    print(io,"    r    : $(m.r  ) \n")
    print(io,"    pr   : $(m.pr ) \n")
    print(io,"    xsr  : $(m.xsr) \n")
    print(io,"    U    : $(m.U) \n")
end

function show(io::IO, m::Model)
    @printf(io,"Region: θu=%06.3f, θr=%06.3f, ϕ=%1.4f, cityarea=%1.4f, rel_cityarea=%1.2f, area=%1.2f, Lu=%1.3f, Lr=%1.3f, pop=%1.3f, pr=%1.3f",m.θu, m.θr, m.ϕ, cityarea(m),rel_cityarea(m), area(m), m.Lu, m.Lr,pop(m),m.pr)
end


pop(m::Model) = m.Lu + m.Lr
cityarea(m::Model) = m.ϕ^2 * π
area(m::Model) = cityarea(m) + m.Sr + m.Srh

"""
	rel_cityarea(m::Model)

Computes city area relative to total rural land use. This is the model moment related to the data moment *artificialized land relative to agricultural land in 2015*.
We measure this to be 18.3% in the data.
"""
rel_cityarea(m::Model) = cityarea(m) / (m.Sr + m.Srh)


"""
	update!(m::Region,p::Param,x::Vector{Float64})

update a single [`Region`](@ref) with a [`Param`](@ref) at choices `x`.
"""
function update!(m::Region,p::Param,x::Vector{Float64}; Lu...)
	# println(x)
	m.ρr    = x[1]
	m.ϕ    = x[2]
	m.r    = x[3]   # land rent
	m.Lr   = x[4]   # employment in rural sector
	m.pr   = x[5]   # relative price rural good
	m.Sr   = x[6]   # amount of land used in rural production
	p.θu   = x[7]   # keep theta as partof param vector for simplicity
	p.θr   = x[8]
	p.θrt[p.it] = p.θr  # fill out time series
	p.θut[p.it] = p.θu  # fill out time series
	m.θr = p.θr  # fill out time series
	m.θu = p.θu  # fill out time series

	m.ϕ10 = m.ϕ * 0.1
	m.ϕ90 = m.ϕ * 0.9
	m.ϕbreaks = range(0.0, stop = m.ϕ, length = p.int_bins+1)
	m.ϕmids = midpoints(collect(m.ϕbreaks))

	m.speedbreaks = [0.0]
	for i in p.speed_thresholds
		xspeed = lofmode(i,p,m.Lu)
		if xspeed < m.ϕ
			push!(m.speedbreaks, xspeed)
		else
			push!(m.speedbreaks, m.ϕ)
			break
		end
		if i == p.speed_thresholds[end]
			push!(m.speedbreaks, m.ϕ)
		end
	end



	# update equations
	if length(Lu) == 0
		m.Lu   = p.L - m.Lr   # employment in urban sector
	else
		m.Lu = Lu[1]
	end
	m.wu0  = wu0(m.Lu,p)   # wage rate urban sector at city center (distance = 0)
	m.wr   = m.wu0 - τ(m.ϕ,m.ϕ,p,m.Lu)
	# m.wr   = foc_Lr(m.Lr / m.Sr , m.pr, p)
	# m.ρr   = foc_Sr(m.Lr / m.Sr , m.pr, p)
	# m.ρr   = 0.059
	# m.ϕ    = getfringe(p.θu, m.wr ,p)

	m.cityarea = cityarea(m)
	m.citydensity = m.Lu / m.cityarea

	m.xsr  = xsr(p,m)
	m.Srh  = Srh(p,m)

	# now can update rel_cityarea
	m.rel_cityarea = rel_cityarea(m)

	m.qr   = qr(p,m)
	m.q0   = q(0.0,p,m)
	m.qq1   = q(m.ϕ / 5,p,m)
	m.ρ0   = ρ(0.0,p,m)
	m.Hr   = H(m.ϕ,p,m)
	m.hr   = h(m.ϕ,p,m)
	m.H0   = H(0.0,p,m)
	m.h0   = h(0.0,p,m)
	m.dbar   = NaN
	m.d0_   = D(0.0,p,m)
	m.dr   = D(m.ϕ,p,m)
	# compute consumption at locations 0 and 1 to check both positive in both sectors.
	m.cr0 = cr(0.0,p,m)
	m.cr1 = cr(m.ϕ,p,m)
	m.cu0 = cu(0.0,p,m)
	m.cu1 = cu(m.ϕ,p,m)
	# m.cr01 = (cr(0.0,p,m)-p.cbar, cr(1.0,p,m)-p.cbar)
	# m.cu01 = (cu(0.0,p,m)       , cu(1.0,p,m)       )
	m.U    = ((m.cr0 .>= 0.0) && (m.cr1 .>= 0.0) && (m.cu0 .>= 0.0) && (m.cu1 .>= 0.0) ) ? utility(0.0,p,m) : NaN
	# if !all((m.cr01 .>= 0.0) .* (m.cu01 .>= 0.0) )
	# println("neg cons")
	# end
	# display(m)
	m.nodes[:]        .= (m.ϕ   + 0.0) / 2 .+ (m.ϕ / 2)   .* m.inodes   # maps [-1,1] into [0,ϕ]
	m.nodes_center[:] .= (p.ϕ1  + 0.0) / 2 .+ (p.ϕ1 / 2)  .* m.inodes   # maps [-1,1] into [0,ϕ1]
	m.nodes_10[:]     .= (m.ϕ10 + 0.0) / 2 .+ (m.ϕ10 / 2) .* m.inodes   # maps [-1,1] into [0,ϕ/10]
	m.nodes_90[:]     .= (m.ϕ90 + m.ϕ) / 2 .+ ((m.ϕ - m.ϕ90) / 2 .* m.inodes)   # maps [-1,1] into [9ϕ/10, ϕ]

	for ib in 1:p.int_bins
		m.nodes_bins[ib,:] .= (m.ϕbreaks[ib] + m.ϕbreaks[ib+1]) / 2 .+ ((m.ϕbreaks[ib+1] - m.ϕbreaks[ib]) / 2 .* m.inodes)
	end

	for ib in 1:(length(m.speedbreaks)-1)
		m.nodes_speeds[ib,:] .= (m.speedbreaks[ib] + m.speedbreaks[ib+1]) / 2 .+ ((m.speedbreaks[ib+1] - m.speedbreaks[ib]) / 2 .* m.inodes)
	end
	
	integrate!(m,p)
	# m.d0   = D(p.ϕ1,p,m)
	m.d0   = m.iDensity_center

	m.mode0 = mode(0.01 * m.ϕ,p,m.Lu)
	m.ctime0 = 0.01 * m.ϕ / m.mode0
	m.modeϕ = mode(m.ϕ,p,m.Lu)
	m.ctimeϕ = m.ϕ / m.modeϕ

	# income measures
	m.pcy = pcy(m,p)
	m.GDP = GDP(m,p)
	m.y   = y(m,p)
	m.Yu  = Yu(m,p)
	m.Yr  = Yr(m,p)

	m.ρ0_y = m.ρ0 / m.y
	m.ρr_y = m.ρr / m.y

	# compute aggregate Consumption shares
	m.Cu = m.icu + m.Lr * cur(p,m)  # urban cons inside city plus total urban cons in rural
	m.Cr = m.pr * (m.icr + m.Lr * cr(m.ϕ,p,m))  # rural cons inside city plus total rural cons in rural
	m.Ch = m.ihexp + m.qr * h(m.ϕ,p,m) * m.Lr   # rural cons inside city plus total rural cons in rural
	m.C  = m.Cu + m.Cr + m.Ch
end


"""
	integrate!(m::Region,p::Param)

perform numerical integration using gauss-laguerre. 
"""
function integrate!(m::Region,p::Param)

	two_π_l       = 2π .* m.nodes
	two_π_lcenter = 2π .* m.nodes_center
	two_π_l10     = 2π .* m.nodes_10
	two_π_l90     = 2π .* m.nodes_90

	m.icu_input       = (m.ϕ/2) * sum(m.iweights[i] * two_π_l[i] * cu_input(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.iDensity        = (m.ϕ/2) * sum(m.iweights[i] * two_π_l[i]        * D(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.iDensity_center = (p.ϕ1/2) * sum(m.iweights[i] * two_π_lcenter[i] * D(m.nodes_center[i],p,m) for i in 1:p.int_nodes)[1]/(p.ϕ1^2 * π)
	m.iDensity_q10    = (m.ϕ10 / 2) * sum(m.iweights[i] * two_π_l10[i] * D(m.nodes_10[i],p,m) for i in 1:p.int_nodes)[1]/(m.ϕ10^2 * π)
	ring90 = (m.ϕ )^2 * π - (m.ϕ90)^2 * π
	m.iDensity_q90    = ((m.ϕ - m.ϕ90) / 2) * sum(m.iweights[i] * two_π_l90[i] * D(m.nodes_90[i],p,m) for i in 1:p.int_nodes)[1]/ring90
	m.icu       = (m.ϕ/2) * sum(m.iweights[i] * two_π_l[i] * cu(m.nodes[i],p,m) * D(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.icr       = (m.ϕ/2) * sum(m.iweights[i] * two_π_l[i] * cr(m.nodes[i],p,m) * D(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.iτ        = (m.ϕ/2) * sum(m.iweights[i] * two_π_l[i] * (m.wu0 - w(m.Lu,m.nodes[i],m.ϕ,p)) * D(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.iq        = (m.ϕ/2) * sum(m.iweights[i] * two_π_l[i] * ρ(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.iy        = (m.ϕ/2) * sum(m.iweights[i] * two_π_l[i] * w(m.Lu,m.nodes[i],m.ϕ,p) * D(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.ihexp     = (m.ϕ/2) * sum(m.iweights[i] * two_π_l[i] * (q(m.nodes[i],p,m) * h(m.nodes[i],p,m) * D(m.nodes[i],p,m)) for i in 1:p.int_nodes)[1]

	# integrate density by distance bin
	allrings = diff(m.ϕbreaks .^2 .* π)
	for ib in 1:p.int_bins
		@views nb = m.nodes_bins[ib,:]
		m.iDensities[ib] = ((m.ϕbreaks[ib+1] - m.ϕbreaks[ib]) / 2) * sum(m.iweights[i] * 2π * nb[i] * D(nb[i],p,m) for i in 1:p.int_nodes)[1] / allrings[ib]
	end

	# integrate density by speed bin
	speedrings = diff(m.speedbreaks .^2 .* π)
	for ib in 1:length(speedrings)
		@views nb = m.nodes_speeds[ib,:]
		m.iDensitySpeeds[ib] = ((m.speedbreaks[ib+1] - m.speedbreaks[ib]) / 2) * sum(m.iweights[i] * 2π * nb[i] * D(nb[i],p,m) for i in 1:p.int_nodes)[1] # / speedrings[ib]
	end

	# averages
	m.qbar      = (m.ϕ/(2 * m.Lu)) * sum(m.iweights[i] * two_π_l[i] * (q(m.nodes[i],p,m) * D(m.nodes[i],p,m)) for i in 1:p.int_nodes)[1]
	m.imode     = (m.ϕ/(2 * m.Lu)) * sum(m.iweights[i] * two_π_l[i] * (mode(m.nodes[i],p,m.Lu) * D(m.nodes[i],p,m)) for i in 1:p.int_nodes)[1]
	m.ictime    = (m.ϕ/(2 * m.Lu)) * sum(m.iweights[i] * two_π_l[i] * ((m.nodes[i] / mode(m.nodes[i],p,m.Lu)) * D(m.nodes[i],p,m)) for i in 1:p.int_nodes)[1]

end


"""
	Eqsys!(F::Vector{Float64},m::Region,p::Param)

compute system of equations for the general (with flexible ϵ).
"""
function Eqsys!(F::Vector{Float64},m::Region,p::Param)

	F[1] = m.wr - foc_Lr(m.Lr / m.Sr , m.pr, p)
	F[2] = m.ρr - foc_Sr(m.Lr / m.Sr , m.pr, p)

	# city size - Urban population relationship: equation (19)
	F[3] = m.Lu - m.iDensity

	#  total land rent per capita: equation (23)
	F[4] = m.iq + m.ρr * (m.Sr + m.Srh) - m.r * p.L

	# land market clearing: after equation (20)
	F[5] = p.S - p.λ - m.cityarea - m.Sr - m.Srh
	# F[1] = m.ρr - 0.059


	# urban goods market clearing. equation before (25) but not in per capita terms
	#      rural cu cons + urban cu cons + rural constr input + urban constr input + commuting - total urban production
	F[6] = m.Lr * cur(p,m) + m.icu + m.Srh * cu_input(m.ϕ,p,m) + m.icu_input + m.iτ - wu0(m.Lu, p)*m.Lu

end


function solve!(F,x,p::Param,m::Model)
	# println(x)
	if any( x .< 0 )
		# F[:] .= PEN
	else
		update!(m,p,x)
		# if isa(m,FModel)
		# 	# m.xtrace = hcat(m.xtrace,x)
		# end
		# # try
			Eqsys!(F,m,p)
		# catch
		# 	@warn "error in eqsys"
		# end
	end
end


#
# Model Component Functions


a(Lu::Float64,p::Param) = p.a * Lu^p.ηa
cτ(Lu::Float64,p::Param) = (0.5 * a(Lu,p))^2 / (2 * p.ζ)


mode(l::Float64,p::Param,Lu::Float64) = ((2*p.ζ)/cτ(Lu,p))^(1/(1+p.ηm)) * l^((1 - p.ηl)/(1+p.ηm)) * wu0(Lu, p)^((1 - p.ηw)/(1+p.ηm))
lofmode(m::Float64, p::Param,Lu::Float64) = ((2*p.ζ)/cτ(Lu,p))^(-1/(1-p.ηl)) * m^((1+p.ηm) / (1-p.ηl)) * wu0(Lu, p)^(-(1 - p.ηw)/(1-p.ηl))

γ(l::Float64,ϕ::Float64,p::Param) = p.γ / (1.0 + ϵ(l,ϕ,p))

"commuting distance"
d(l::Float64,ϕ::Float64, p::Param) = p.d1 * ϕ + l / (1 + p.d2 * ϕ)

"commuting cost: location x → cost"
τ(x::Float64,ϕ::Float64,p::Param,Lu::Float64) = a(Lu,p) * wu0(Lu, p)^(p.ξw) * d(x,ϕ,p)^(p.ξl)

"inverse commuting cost. cost x → location. Notice we don't consider that cost is zero beyond ϕ: we want to find ϕ here to start with."
invτ(Lu::Float64,wr::Float64,wu::Float64,p::Param) = ( (wu - wr) / ( a(Lu,p) * wu^(p.ξw)) )^(1.0/p.ξl)


"""
Get Fringe from indifference condition

At the fringe ``\\phi`` we have the condition

```math
w(0) - \\tau(\\phi) = w_r
```

which can be rearranged to obtain a map from ``w(0) - w_r = \\tau(\\phi)``.

The function takes ``w(0) - w_r`` as argument `x`. then we give ``\\tau(\\phi)``
	to its inverse function to get back ``\\phi``
"""
getfringe(Lu::Float64,w0::Float64,wr::Float64,p::Param) = w0 > wr ? invτ(Lu,wr,w0,p) : 0.0

"urban wage at location ``l``"
wu(Lu::Float64,ϕ::Float64,l::Float64,p::Param) = wu0(Lu,p) .- τ(l,ϕ,p,Lu)

"wage at location ``l``"
w(Lu::Float64,l::Float64,ϕ::Float64,p::Param) = l >= ϕ ? wr(Lu,ϕ,p) : wu(Lu,ϕ,l,p)

"urban wage at center"
wu0(Lu::Float64,p::Param) = p.Ψ * p.θu * Lu^p.η

"rural wage from indifference condition at ϕ. Eq (11)"
wr(Lu::Float64,ϕ::Float64,p::Param) = wu0(Lu,p) .- τ(ϕ,ϕ,p,Lu)

"FOC of rural firm wrt labor Lr"
foc_Lr(L_over_S::Float64,pr::Float64, p::Param) = p.α * pr * p.θr * (p.α + (1-p.α)*( 1.0/ L_over_S )^((p.σ-1)/p.σ))^(1.0 / (p.σ-1))

"FOC of rural firm wrt land Sr"
foc_Sr(L_over_S::Float64,pr::Float64, p::Param) = (1-p.α)* pr * p.θr * (p.α * (L_over_S)^((p.σ-1)/p.σ) + (1-p.α))^(1.0 / (p.σ-1))


"excess subsistence urban worker"
xsu(l::Float64,p::Param,m::Model) = w(m.Lu,l,m.ϕ,p) .+ m.r .- m.pr .* p.cbar .+ p.sbar

"excess subsistence rural worker"
xsr(p::Param,m::Model) = m.r + wr(m.Lu,m.ϕ,p) - m.pr * p.cbar + p.sbar

"Residential Land in Rural sector. equation (20)"
Srh(p::Param,m::Model) = m.Lr * (γ(m.ϕ,m.ϕ,p) * m.xsr) / m.ρr

"optimal urban good consumption at location ``l``. Equation (8)"
cu(l::Float64,p::Param,m::Model) = (1.0 - p.γ)*(1.0 - p.ν)*(w(m.Lu,l,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) - p.sbar

"optimal rural good consumption at location ``l``. Equation (7) divided by p"
cr(l::Float64,p::Param,m::Model) = (1.0 - p.γ) * p.ν * ((w(m.Lu,l,m.ϕ,p) .+ (m.r + p.sbar - m.pr * p.cbar)) ./ m.pr) + p.cbar

"urban good consumption in rural sector. Equation (8)"
cur(p::Param,m::Model) = cu(m.ϕ,p,m)

"cost of construction at location ``l``"
cfun(l::Float64,p::Param) = p.c0 + p.c1 * l + p.c2 * l^2
cost(l::Float64,ϕ::Float64,p::Param) = l >= ϕ ? cfun(ϕ,p) : cfun(l,p)

"housing supply shifter at ``l``"
χ(l::Float64,ϕ::Float64,p::Param) = (1.0 / cost(l,ϕ,p))^ϵ(l,ϕ,p)

"housing supply elasticity at ``l``"
ϵ(l::Float64,ϕ::Float64,p::Param) = p.ϵflat ? p.ϵr : p.ϵr * l / ϕ + p.ϵs * (ϕ - l)/ϕ

"house price function at ``l``. equation (12)"
q(l::Float64,p::Param,m::Model) = m.qr * (xsu(l,p,m) / m.xsr).^(1.0/p.γ)

"land price function at ``l``. equation (15)"
ρ(l::Float64,p::Param,m::Model) = (q(l,p,m).^(1.0 + ϵ(l,m.ϕ,p))) / (1.0 + ϵ(l,m.ϕ,p))

# "rural house price from land price"
# function qr(p::Param,m::Model)
# 	( (1+p.ϵr) * m.ρr * cfun(m.ϕ,p)^p.ϵr ) ^(1.0/(1+p.ϵr))
# end

"rural house price from land price"
function qr(p::Param,m::Model)
	( (1+p.ϵr) * m.ρr ) ^(1.0/(1+p.ϵr))
end


"housing demand at location ``l``"
# h(l::Float64,p::Param,m::Model) = (p.γ / m.qr) * xsu(l,p,m)^((p.γ-1)/p.γ) * m.xsr^(1/p.γ)
h(l::Float64,p::Param,m::Model) = p.γ * (w(m.Lu,l,m.ϕ,p) + m.r - m.pr * p.cbar + p.sbar) / q(l,p,m)

"housing supply at location ``l``"
H(l::Float64,p::Param,m::Model) = q(l,p,m).^ϵ(l,m.ϕ,p)

"Population Density at location ``l``"
D(l::Float64,p::Param,m::Model) = H(l,p,m) / h(l,p,m)

"Population Density at location ``l``. second version, independent of Lu"
function D2(l::Float64,p::Param,r::Region)
	γl = 1.0 / γ(l,r.ϕ,p)
	r.ρr * (γl * (r.wr + r.r + p.sbar - r.pr*p.cbar)^(-γl) * (w(r.Lu,l,r.ϕ,p) + r.r + p.sbar - r.pr*p.cbar)^(γl - 1))
end

"Amount of Urban Good (numeraire) required to build housing at ``l``"
cu_input(l::Float64,p::Param,m::Model) = ϵ(l,m.ϕ,p) / (1+ϵ(l,m.ϕ,p)) .* q(l,p,m) .* H(l,p,m)

"Utility function equation (5)"
utility(l::Float64,p::Param,m::Model) = (cr(l,p,m) - p.cbar)^(p.ν * (1-p.γ)) * (cu(l,p,m) + p.sbar)^((1-p.ν) * (1-p.γ)) * h(l,p,m)^p.γ


"aggregate per capita income"
pcy(m::Model,p::Param) = m.r + m.wr * m.Lr / p.L + m.iy / p.L


"GDP"
GDP(m::Model,p::Param) = (m.pr * Yr(m,p) + Yu(m,p)) / p.L

"disposable income: GDP net of commuting costs"
y(m::Model,p::Param) = GDP(m,p) - m.iτ / p.L


"Production of Rural Good"
function Yr(m::Model,p::Param)
	σ1 = (p.σ - 1) / p.σ
	σ2 = p.σ / (p.σ - 1)
	p.θr * (p.α * (m.Lr^σ1) + (1-p.α) * (m.Sr^σ1) )^σ2
end

"Production of Urban Good"
Yu(m::Model,p::Param) = p.θu * m.Lu

"Rural Market Clearing"
Rmk(m::Model,p::Param) = m.icr + m.Lr * cr(m.ϕ,p,m) - Yr(m,p)


"Obtain a Time Series for a single region as a DataFrame"
function dataframe(M::Vector{T},p::Param) where T <: Model
	tt = length(p.T)
	df = DataFrame(year = p.T)
	for fi in setdiff(fieldnames(eltype(M)),(:cr01,:cu01,:inodes,:iweights,:nodes))
		df[!,fi] = [getfield(M[it],fi) for it in 1:length(p.T)]
	end
	df.area = [area(M[it]) for it in 1:length(p.T)]
	df.pop  = [pop(M[it]) for it in 1:length(p.T)]

	# add rural land rents
	df.rr = df.ρr .* (df.Sr .+ df.Srh)

	# compute commuting cost at initial fringe in each period
	initϕ = df.ϕ[1]
	df.τ_ts = zeros(tt)
	df.rural_emp_model = zeros(tt)
	df.p_laspeyres = zeros(tt)
	df.p_paasche = zeros(tt)
	df.p_growth = zeros(tt)
	df.p_index = zeros(tt)
	df[1,:p_index] = M[1].pr
	for i in 1:tt
		setperiod!(p,i)
		df[i, :τ_ts] = τ(initϕ,initϕ, p, M[i].Lu)
		df[i, :rural_emp_model] = M[i].Lr / p.L
		if i > 1
			df[i, :p_laspeyres] = ( M[i].pr * M[i-1].Yr + M[i-1].Yu ) / ( M[i-1].pr *  M[i-1].Yr + M[i-1].Yu )
			df[i, :p_paasche]   = ( M[i].pr * M[i].Yr + M[i].Yu ) / ( M[i-1].pr *  M[i].Yr + M[i].Yu )
			df[i, :p_growth]   = sqrt(df[i, :p_paasche]) * sqrt(df[i, :p_laspeyres])
			df[i, :p_index]      = df[i-1, :p_index] * df[i, :p_growth]

		end
	end
	df[!,:r_real] = df[!,:r] .* df[!,:pop] ./ df[!, :p_index]
	df[!,:rr_real] = df[!,:rr] ./ df[!, :p_index]
	df[!,:ru_real] = df[!,:iq] ./ df[!, :p_index]
	df[!,:qr_real] = df[!,:qr] ./ df[!, :p_index]
	df[!,:qbar_real] = df[!,:qbar] ./ df[!, :p_index]
	df[!,:ρ0_real] = df[!,:ρ0] ./ df[!, :p_index]

	# normalizations / transforms
	# normalize densities at first period
	sort!(df, :year)
	dens = select(df,:d0,  :dr, :citydensity)
	mapcols!(x -> x ./ x[1],dens)
	rename!(dens, [:d0_n, :dr_n, :avgd_n] )
	df = [df dens]

	# normalize pr at 1950
	i1950 = argmin( abs.(p.moments.year .- 1950) )
	transform!(df, :pr => (x -> x ./ x[i1950]) => :pr_1950)

	# compute spending shares
	transform!(df, [:Ch, :C] => ((a,b) -> a ./b) => :hshare,
	               [:Cu, :C] => ((a,b) -> a ./b) => :ushare,
	               [:Cr, :C] => ((a,b) -> a ./b) => :rshare)

	transform!(df, :Lr => (x -> x ./ (p.Lt .* p.K)) => :Lr_n, :Lu => (x -> x ./ (p.Lt .* p.K)) => :Lu_n)

	df
end

function pdiff(M::Vector{Model}, v::Symbol; t1=1,t2=7)
	(getfield(M[t2],v) - getfield(M[t1],v)) / getfield(M[t1],v)
end

"normalize a vector wrt it's first element"
firstnorm(x) = x ./ x[1]

"normalize a vector wrt it's maximal element"
maxnorm(x) = x ./ maximum(x)