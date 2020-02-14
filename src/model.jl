
abstract type Model end


"""
Region type

Urban expansion dodel with flexible commuting cost and differential supply elasticity.
This represents a single region. A region contains one potential urban area, and a rural area.
A `Region` can be part of a wider `Country`.
"""
mutable struct Region <: Model
	ρr   :: Float64   # land price in rural sector
	qr   :: Float64     # housing price in rural sector
	Lr   :: Float64   # employment in rural sector
	Lu   :: Float64   # employment in urban sector
	wu0  :: Float64   # wage in urban sector (at center)
	wr   :: Float64   # wage in rural sector
	Sr   :: Float64   # Amount of land used in rural production
	Srh  :: Float64   # Amount of land used for rural housing
	r    :: Float64   # per capita land rental income
	pr   :: Float64   # relative price of rural good
	ϕ    :: Float64   # size of the city
	xsr  :: Float64  # excess subsistence rural worker
	U    :: Float64  # common utility level
	cr01 :: Tuple{Float64,Float64}  # consumption at locations 0 and 1. temps to check.
	cu01 :: Tuple{Float64,Float64} # consumption at locations 0 and 1. temps to check.
	# integration setup
	inodes   :: Vector{Float64}  # theoretical integration nodes
	iweights :: Vector{Float64}  # int weights
	nodes    :: Vector{Float64}  # points where to evaluate integrand (inodes scaled into [0,ϕ])

	# resulting integrals from [0,ϕ]
	icu_input :: Float64   # ∫ cu_input(l) dl
	iDensity  :: Float64   # ∫ D(l) dl
	icu       :: Float64   # ∫ cu(l) D(l) dl
	icr       :: Float64   # ∫ cr(l) D(l) dl
	iτ        :: Float64   # ∫ τ(l) D(l) dl
	iq        :: Float64   # ∫ ρ(l) dl
	iy        :: Float64   # ∫ w(l) D(l) dl
	function Region(p::Param)
		# creates a model fill with NaN
		m      = new()
		m.ρr   = NaN
		m.qr   = NaN
		m.Lr   = NaN
		m.Lu   = NaN
		m.wu0  = NaN
		m.wr   = NaN
		m.Sr   = NaN
		m.Srh  = NaN
		m.r    = NaN
		m.pr   = NaN
		m.ϕ    = NaN
		m.xsr  = NaN
		m.U    = NaN
		m.cr01 = (NaN,NaN)
		m.cu01 = (NaN,NaN)
		m.inodes, m.iweights = gausslegendre(p.int_nodes)
		m.nodes = zeros(p.int_nodes)
		icu_input = NaN
		iDensity  = NaN
		icu       = NaN
		iτ        = NaN
		iq        = NaN
		iy        = NaN
		return m
	end
end

function show(io::IO, ::MIME"text/plain", m::Region)
    print(io,"Single Region:\n")
    print(io,"    pop   : $(pop(m)) \n")
    print(io,"    Lr   : $(m.Lr ) \n")
    print(io,"    Lu   : $(m.Lu ) \n")
    print(io,"    area : $(round(area(m),digits=2)) \n")
    print(io,"    ϕ    : $(m.ϕ  ) \n")
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

function show(io::IO, m::Region)
    # print(io,"Region: ϕ=$(round(m.ϕ,digits=3)), pop=$(pop(m)), area=$(round(area(m),digits=2))")
    @printf(io,"Region: ϕ=%1.3f, pop=%1.3f, area=%1.2f",m.ϕ, pop(m), area(m))
end


pop(m::Model) = m.Lu + m.Lr
area(m::Model) = m.ϕ + m.Sr + m.Srh


"""
	update!(m::Region,p::Param,x::Vector{Float64})

update a single region a parameter vector at choices `x`.
"""
function update!(m::Region,p::Param,x::Vector{Float64})
	m.r    = x[1]   # land rent
	m.Lr   = x[2]   # employment in rural sector
	m.pr   = x[3]   # relative price rural good
	m.Sr   = x[4]   # amount of land used in rural production

	# update equations
	m.Lu   = p.L - m.Lr   # employment in urban sector
	m.wu0  = wu0(m.Lu,p)   # wage rate urban sector at city center (distance = 0)
	m.wr   = foc_Lr(m.Lr / m.Sr , m.pr, p)
	m.ρr   = foc_Sr(m.Lr / m.Sr , m.pr, p)
	m.ϕ    = invτ(m.wr / p.θu,p)

	m.xsr  = xsr(p,m)
	m.Srh  = Srh(p,m)
	m.qr   = qr(p,m)
	# compute consumption at locations 0 and 1 to check both positive in both sectors.
	m.cr01 = (cr(0.0,p,m)-p.cbar, cr(1.0,p,m)-p.cbar)
	m.cu01 = (cu(0.0,p,m)       , cu(1.0,p,m)       )
	m.U    = all( (m.cr01 .>= 0.0) .* (m.cu01 .>= 0.0) ) ? utility(0.0,p,m) : NaN
	# if !all((m.cr01 .>= 0.0) .* (m.cu01 .>= 0.0) )
	# println("neg cons")
	# end
	# display(m)
	m.nodes[:] .= m.ϕ / 2 .+ (m.ϕ / 2) .* m.inodes   # maps [-1,1] into [0,ϕ]
	integrate!(m,p)
end


"""
	integrate!(m::Region,p::Param)

could be made much faster by filling a matrix col-wise with integration nodes and
doing a matmul on it? have to allocate memory though.
"""
function integrate!(m::Region,p::Param)

	m.icu_input = (m.ϕ/2) * sum(m.iweights[i] * cu_input(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.iDensity  = (m.ϕ/2) * sum(m.iweights[i] * D(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.icu       = (m.ϕ/2) * sum(m.iweights[i] * cu(m.nodes[i],p,m) * D(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.icr       = (m.ϕ/2) * sum(m.iweights[i] * cr(m.nodes[i],p,m) * D(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.iτ        = (m.ϕ/2) * sum(m.iweights[i] * τ(m.nodes[i],m.ϕ,p) * D(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.iq        = (m.ϕ/2) * sum(m.iweights[i] * ρ(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.iy        = (m.ϕ/2) * sum(m.iweights[i] * w(m.Lu,m.nodes[i],m.ϕ,p) * D(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	# @debug "integrate!" icu_input=m.icu_input iDensity=m.iDensity icu=m.icu iτ=m.iτ iq=m.iq phi=m.ϕ
	# @assert m.iDensity > 0 "integral of density is negative"
end

"""
	Eqsys!(F::Vector{Float64},m::Region,p::Param)

compute system of equations for the general (with flexible ϵ).
"""
function Eqsys!(F::Vector{Float64},m::Region,p::Param)
	# land market clearing: after equation (20)
	F[1] = p.S - p.λ - m.ϕ - m.Sr - m.Srh

	# city size - Urban population relationship: equation (19)
	F[2] = m.Lu - m.iDensity

	#  total land rent per capita: equation (23)
	F[3] = m.r * p.L - m.iq - m.ρr * (p.S - m.ϕ)

	# urban goods market clearing. equation before (25) but not in per capita terms
	#      rural cu cons + urban cu cons + rural constr input + urban constr input + commuting - total urban production
	F[4] = m.Lr * cur(p,m) + m.icu + m.Srh * cu_input(m.ϕ,p,m) + m.icu_input + m.wu0 * m.iτ - wu0(m.Lu, p)*m.Lu
end


"""
Model with linear commuting cost and fixed supply elasticity.
This model admits closed form integrals.
"""
mutable struct FModel <: Model
	ρr      :: Float64   # land price in rural sector
	qr      :: Float64     # housing price in rural sector
	Lr      :: Float64   # employment in rural sector
	Lu      :: Float64   # employment in urban sector
	wu0     :: Float64   # wage in urban sector (at center)
	wr      :: Float64   # wage in rural sector
	Sr      :: Float64   # Amount of land used in rural production
	Srh     :: Float64   # Amount of land used for rural housing
	r       :: Float64   # per capita land rental income
	pr      :: Float64   # relative price of rural good
	ϕ       :: Float64   # size of the city
	xsr      :: Float64  # excess subsistence rural worker
	U        :: Float64  # common utility level
	cr01     :: Tuple{Float64,Float64}  # consumption at locations 0 and 1. temps to check.
	cu01     :: Tuple{Float64,Float64} # consumption at locations 0 and 1. temps to check.
	# Ftrace :: Matrix{Float64}
	# xtrace :: Matrix{Float64}
	function FModel(p::Param)
		# creates a model fill with NaN
		m     = new()
		m.ρr  = NaN
		m.qr  = NaN
		m.Lr  = NaN
		m.Lu  = NaN
		m.wu0 = NaN
		m.wr  = NaN
		m.Sr  = NaN
		m.Srh = NaN
		m.r   = NaN
		m.pr  = NaN
		m.ϕ   = NaN
		# m.xsr  = NaN
		# m.U    = NaN
		m.xsr    = NaN
		# m.Ftrace = zeros(6,1)
		# m.xtrace = zeros(6,1)
		return m
	end
end

"""
	update!(m::FModel,p::Param,x::Vector{Float64})

update a `Fmodel` from a choice vector x
"""
function update!(m::FModel,p::Param,x::Vector{Float64})
	m.r    = x[1]   # land rent
	m.Lr   = x[2]   # employment in rural sector
	m.pr   = x[3]   # relative price rural good
	m.Sr   = x[4]   # amount of land used in rural production

	# update params
	σ1 = (p.σ-1)/p.σ
	σ2 = 1.0 / (p.σ-1)
	γ2 = p.γ / (1 + p.ϵr)


	# update equations
	m.Lu   = p.L - m.Lr   # employment in urban sector
	m.wu0  = wu0(m.Lu,p)   # wage rate urban sector at city center (distance = 0)
	m.wr   = foc_Lr(m.Lr / m.Sr , m.pr, p)
	# m.wr   = p.α * m.pr * p.θr * (p.α + (1-p.α)*(m.Sr / m.Lr)^σ1)^σ2
	# m.ρr   = (1-p.α)*m.pr * p.θr * (p.α * (m.Lr / m.Sr)^σ1 + (1-p.α))^σ2
	m.ρr   = foc_Sr(m.Lr / m.Sr , m.pr, p)
	m.ϕ = invτ(m.wr / p.θu,p)

	m.xsr  = m.wr + m.r - m.pr * p.cbar - p.sbar
	# m.xsu = m.wu0 + m.r - m.pr * p.cbar - p.sbar
	m.Srh  = ( γ2 / m.ρr ) * m.xsr * m.Lr
	m.qr   = qr(p,m)

end

"""
	Eqsys!(F::Vector{Float64},m::FModel,p::Param)

compute system of equations for fixed elasticities. uses closed form solutions
	for integrals.
"""
function Eqsys!(F::Vector{Float64},m::FModel,p::Param)
	σ1 = (p.σ-1)/p.σ
	σ2 = 1.0 / (p.σ-1)
	γ2 = p.γ / (1 + p.ϵr)

	# land market clearing: after equation (18)
	F[1] = p.S - p.λ - m.ϕ - m.Sr - m.Srh


	# city size - Urban population relationship: analytic integral solution to equation (16)
	w2 = p.θu * m.Lu^p.η
	τϕ = 1-p.τ * m.ϕ
	w1 = p.Ψ * w2 * τϕ
	uu = w2+m.r - m.pr * p.cbar - p.sbar
	ur = w1+m.r - m.pr * p.cbar - p.sbar
	xx = m.r - m.pr * p.cbar - p.sbar


	F[2] = 1+p.τ * w2 * m.Lu/(m.ρr)-(uu/(w2*(1-p.τ*m.ϕ)+m.r-m.pr*p.cbar+p.sbar))^(1/γ2)

	#  total land rent per capita: closed form solution
	F[3] = m.ρr * m.Sr + m.ρr * m.Srh +
	       γ2 * m.ρr/(p.τ*(1+γ2)*w2)*(uu^(1/γ2+1)/((w2*τϕ+m.r-m.pr * p.cbar+p.sbar)^(1/γ2)) -
		   (w2*(1-p.τ * m.ϕ)+m.r-m.pr * p.cbar+p.sbar)) - m.r*p.L


	# urban goods market clearing.
	chi2 = 1.0

	F[4] = m.Lr * (1-p.γ)*(1-p.ν)* ur +
	chi2*(1-p.γ)*(1-p.ν)*m.ρr/((1+γ2)*w2*p.τ)*((w2+xx)^(1/γ2+1)/((w2*τϕ+xx)^(1/γ2))-(w2*τϕ+xx)) +
	chi2*γ2*m.ρr/((1+γ2)*p.τ*w2)*(((w2+xx)^(1/γ2+1))/((w2*τϕ+xx)^(1/γ2))-(w2*τϕ+xx)) -
	m.ϕ*chi2*m.ρr +
	p.ϵr *m.ρr*m.Srh +
	p.ϵr *chi2*γ2*m.ρr/(p.τ*(1+γ2)*w2)*((w2+xx)^(1/γ2+1)/((w2*τϕ+xx)^(1/γ2))-(w2*τϕ+xx)) -
	p.sbar*p.L -
	p.θu*m.Lu^(1+p.η)

	# record function values
	# m.Ftrace = hcat(m.Ftrace,F)

	# cr = (1.0 - p.γ) * p.ν * ur + m.pr * p.cbar
	# println("cr - p.cbar = $(cr - p.cbar * p.cbar)")
	# F[7] = minerr < 0.0 ? minerr^2 : 0.0
	# F[7] = cr - p.cbar * p.cbar > 0 ? 0.0 : 2*exp(cr - p.cbar * p.cbar)
	# F[7] = 0.0
	# println("penalty = $(F[7])")

	# F[7] = 0.0
end

function solve!(F,x,p::Param,m::Model)
	# println(x)
	if any( x .< 0 )
		F[:] .= PEN
	else
		update!(m,p,x)
		if isa(m,FModel)
			# m.xtrace = hcat(m.xtrace,x)
		end
		# try
			Eqsys!(F,m,p)
		# catch
		# 	@warn "error in eqsys"
		# end
	end
end








# Model Component Functions

γ(l::Float64,ϕ::Float64,p::Param) = p.γ / (1.0 + ϵ(l,ϕ,p))

"commuting cost"
τ(x::Float64,ϕ::Float64,p::Param) = (x > ϕ) ? 0.0 : p.τ * x

"inverse of commuting cost function"
invτ(x::Float64,p::Param) = (x > 1.0) ? 0.0 : (1.0 - x) / p.τ

"urban wage at location ``l``"
wu(Lu::Float64,l::Float64,ϕ::Float64,p::Param) = wu0(Lu,p) * (1.0 .- τ(l,ϕ,p))

"urban wage at location ``l``"
wu(l::Float64,ϕ::Float64,p::Param) = wu0(p) * (1.0 .- τ(l,ϕ,p))

"urban wage at center"
wu0(Lu::Float64,p::Param) = p.Ψ * p.θu * Lu^p.η

"urban wage at center"
wu0(p::Param) = p.Ψ * p.θu

"rural wage from indifference condition at ϕ. Eq (11)"
wr(Lu::Float64,ϕ::Float64,p::Param) = wu0(Lu,p)*(1.0 .- τ(ϕ,ϕ,p))

"FOC of rural firm wrt labor Lr"
foc_Lr(L_over_S::Float64,pr::Float64, p::Param) = p.α * pr * p.θr * (p.α + (1-p.α)*( 1.0/ L_over_S )^((p.σ-1)/p.σ))^(1.0 / (p.σ-1))

"FOC of rural firm wrt land Sr"
foc_Sr(L_over_S::Float64,pr::Float64, p::Param) = (1-p.α)* pr * p.θr * (p.α * (L_over_S)^((p.σ-1)/p.σ) + (1-p.α))^(1.0 / (p.σ-1))


"rural wage from indifference condition at ϕ. Eq (11)"
wr(ϕ::Float64,p::Param) = wu0(p)*(1.0 .- τ(ϕ,ϕ,p))

"wage at location ``l``"
w(Lu::Float64,l::Float64,ϕ::Float64,p::Param) = l >= ϕ ? wr(Lu,ϕ,p) : wu(Lu,l,ϕ,p)

"wage at location ``l`` indep of Lu"
w(l::Float64,ϕ::Float64,p::Param) = l >= ϕ ? wr(ϕ,p) : wu(l,ϕ,p)

"excess subsistence urban worker"
xsu(l::Float64,p::Param,m::Model) = w(m.Lu,l,m.ϕ,p) .+ m.r .- m.pr .* p.cbar .+ p.sbar

"excess subsistence rural worker"
xsr(p::Param,m::Model) = m.r + wr(m.Lu,m.ϕ,p) - m.pr * p.cbar + p.sbar

"Residential Land in Rural sector. equation (20)"
Srh(p::Param,m::Model) = m.Lr * (γ(m.ϕ,m.ϕ,p) * m.xsr) / m.ρr

"optimal urban good consumption at location ``l``. Equation (8)"
cu(l::Float64,p::Param,m::Model) = (1.0 - p.γ)*(1.0 - p.ν)*(w(m.Lu,l,m.ϕ,p) .+ (m.r - m.pr * p.cbar))

"optimal rural good consumption at location ``l``. Equation (7) divided by p"
cr(l::Float64,p::Param,m::Model) = (1.0 - p.γ) * p.ν * ((w(m.Lu,l,m.ϕ,p) .+ (m.r - m.pr * p.cbar)) ./ m.pr) + p.cbar

"urban good consumption in rural sector. Equation (8)"
cur(p::Param,m::Model) = cu(m.ϕ,p,m)

"cost of construction at location ``l``"
cfun(l::Float64,p::Param) = p.c0 + p.c1 * l + p.c2 * l^2
cost(l::Float64,ϕ::Float64,p::Param) = l >= ϕ ? cfun(ϕ,p) : cfun(l,p)

"housing supply shifter at ``l``"
χ(l::Float64,ϕ::Float64,p::Param) = (1.0 / cost(l,ϕ,p))^ϵ(l,ϕ,p)

"housing supply elasticity at ``l``"
ϵ(l::Float64,ϕ::Float64,p::Param) = p.ϵr * exp(-p.ϵs * max(ϕ-l,0.0))
# ϵ(l::Float64,ϕ::Float64,p::Param) = exp(log(p.ϵr)-p.ϵs*log(max(1 + ϕ - l,0.0)))

"house price function at ``l``. equation (12)"
q(l::Float64,p::Param,m::Model) = m.qr * (xsu(l,p,m) / m.xsr).^(1.0/p.γ)

"land price function at ``l``. equation (15)"
# ρ(l::Float64,p::Param,m::Model) = H(l,p,m) .* q(l,p,m) / (1.0 + p.ϵ)
ρ(l::Float64,p::Param,m::Model) = (χ(l,m.ϕ,p) .* q(l,p,m).^(1.0 + ϵ(l,m.ϕ,p))) / (1.0 + ϵ(l,m.ϕ,p))

# "rural house price from land price"
# function qr(p::Param,m::Model)
# 	( (1+p.ϵr) * m.ρr * cfun(m.ϕ,p)^p.ϵr ) ^(1.0/(1+p.ϵr))
# end

"rural house price from land price"
function qr(p::Param,m::Model)
	( (1+p.ϵr) * m.ρr ) ^(1.0/(1+p.ϵr))
end


"housing demand at location ``l``"
h(l::Float64,p::Param,m::Model) = (p.γ / m.qr) * xsu(l,p,m)^((p.γ-1)/p.γ) * m.xsr^(1/p.γ)

"housing supply at location ``l``"
H(l::Float64,p::Param,m::Model) = χ(l,m.ϕ,p) * q(l,p,m).^ϵ(l,m.ϕ,p)

"Population Density at location ``l``"
D(l::Float64,p::Param,m::Model) = H(l,p,m) / h(l,p,m)

"Population Density at location ``l``. second version, independent of Lu"
function D2(l::Float64,p::Param,r::Region)
	γl = 1.0 / γ(l,r.ϕ,p)
	r.ρr * (γl * (r.wr + r.r - r.pr*p.cbar)^(-γl) * (w(r.Lu,l,r.ϕ,p) + r.r - r.pr*p.cbar)^(γl - 1))
end

"Amount of Urban Good (numeraire) required to build housing at ``l``"
cu_input(l::Float64,p::Param,m::Model) = ϵ(l,m.ϕ,p) / (1+ϵ(l,m.ϕ,p)) .* q(l,p,m) .* H(l,p,m)

"Utility function equation (5)"
utility(l::Float64,p::Param,m::Model) = (cr(l,p,m) - p.cbar)^(p.ν * (1-p.γ)) * (cu(l,p,m))^((1-p.ν) * (1-p.γ)) * h(l,p,m)^p.γ


"aggregate per capita income"
pcy(m::Model,p::Param) = m.r + m.wr * m.Lr / p.L + m.iy / p.L

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


"Obtain a Time Series for the Model as a DataFrame"
function dataframe(M::Vector{Region},p::Param)
	df = DataFrame(year = p.T)
	for fi in setdiff(fieldnames(LandUse.Region),(:cr01,:cu01,:inodes,:iweights,:nodes))
		df[!,fi] = [getfield(M[it],fi) for it in 1:length(p.T)]
	end
	df
end
