"""
Model with general commuting cost.
"""
mutable struct Model
	ρr      :: Float64   # land price in rural sector
	qr      :: Float64     # housing price in rural sector
	ϵr      :: Float64     # housing supply elasticity in rural sector
	χr      :: Float64     # housing supply shifter in rural sector
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
	# integration setup
	inodes   :: Vector{Float64}  # theoretical integration nodes
	iweights :: Vector{Float64}  # int weights
	nodes   :: Vector{Float64}  # points where to evaluate integrand (inodes scaled into [0,ϕ])

	# resulting integrals from [0,ϕ]
	icu_input :: Float64   # ∫ cu_input(l) D(l) dl
	iDensity  :: Float64   # ∫ D(l) dl
	icu       :: Float64   # ∫ cu(l) D(l) dl
	iτ        :: Float64   # ∫ τ(l) D(l) dl
	ir        :: Float64   # ∫ ρ(l) D(l) dl
	function Model(p::Param)
		# creates a model fill with NaN
		m     = new()
		m.ρr  = NaN
		m.qr  = NaN
		m.ϵr  = NaN
		m.χr  = NaN
		m.Lr  = NaN
		m.Lu  = NaN
		m.wu0 = NaN
		m.wr  = NaN
		m.Sr  = NaN
		m.Srh = NaN
		m.r   = NaN
		m.pr  = NaN
		m.ϕ   = NaN
		m.xsr  = NaN
		m.U    = NaN
		m.inodes, m.iweights = gausslegendre(p.int_nodes)
		m.nodes = zeros(p.int_nodes)
		icu_input = NaN
		iDensity  = NaN
		icu       = NaN
		iτ        = NaN
		ir        = NaN
		return m
	end
end


function show(io::IO, ::MIME"text/plain", m::Model) 
    print(io,"LandUse Model:\n")
    print(io,"    m.ρr   : $(m.ρr ) \n")
    print(io,"    m.qr : $(m.qr ) \n")
    print(io,"    m.ϵr : $(m.ϵr ) \n")
    print(io,"    m.χr : $(m.χr ) \n")
    print(io,"    m.Lr   : $(m.Lr ) \n")
    print(io,"    m.Lu   : $(m.Lu ) \n")
    print(io,"    m.wu0  : $(m.wu0) \n")
    print(io,"    m.wr   : $(m.wr ) \n")
    print(io,"    m.Sr   : $(m.Sr ) \n")
    print(io,"    m.Srh  : $(m.Srh) \n")
    print(io,"    m.r    : $(m.r  ) \n")
    print(io,"    m.pr   : $(m.pr ) \n")
    print(io,"    m.ϕ    : $(m.ϕ  ) \n")
    print(io,"    m.xsr  : $(m.xsr) \n")
    print(io,"    m.U    : $(m.U) \n")
end


"""
	integrate!(m::Model,p::Param)

could be made much faster by filling a matrix col-wise with integration nodes and 
doing a matmul on it? have to allocate memory though.
"""
function integrate!(m::Model,p::Param)
	
	m.icu_input = (m.ϕ/2) * sum(m.iweights[i] * cu_input(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.iDensity  = (m.ϕ/2) * sum(m.iweights[i] * D(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.icu       = (m.ϕ/2) * sum(m.iweights[i] * cu(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.iτ        = (m.ϕ/2) * sum(m.iweights[i] * τ(m.nodes[i],m.ϕ,p) for i in 1:p.int_nodes)[1]
	m.ir        = (m.ϕ/2) * sum(m.iweights[i] * ρ(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	@debug "integrate!" icu_input=m.icu_input iDensity=m.iDensity icu=m.icu iτ=m.iτ ir=m.ir phi=m.ϕ
	if m.iDensity < 0
		println("D = $([D(m.nodes[i],p,m) for i in 1:p.int_nodes])")
		println("H = $([H(m.nodes[i],p,m) for i in 1:p.int_nodes])")
		println("h = $([h(m.nodes[i],p,m) for i in 1:p.int_nodes])")
	end
end

"""
	update!(m::Model,m0::CD0Model, p::Param)

update the general model from a CD0Model
"""
function update!(m::Model,m0::CD0Model,p::Param)
	m.ρr   = m0.ρr
	m.ϕ    = m0.ϕ 
	m.r    = m0.r 
	m.Lr   = m0.Lr
	m.pr   = m0.pr
	m.Sr   = m0.Sr

	# update equations
	m.Lu   = p.L - m.Lr   # employment in urban sector
	m.wu0  = wu0(m.Lu,p)   # wage rate urban sector at city center (distance = 0)
	m.wr   = wr(m.Lu,m.ϕ,p) # wage rate rural sector
	m.xsr  = xsr(p,m)
	m.qr   = qr(p,m)
	m.Srh  = Srh(p,m)
	m.U    = utility(0.0,p,m)
	m.nodes[:] .= m.ϕ / 2 .+ (m.ϕ / 2) .* m.inodes   # maps [-1,1] into [0,ϕ]

end

"""
	update!(m::Model,p::Param,x::Vector{Float64})

update the general model from a parameter vector
"""
function update!(m::Model,p::Param,x::Vector{Float64})
	m.ρr   = x[1]   # land price in rural sector
	m.ϕ    = x[2]   # city size
	m.r    = x[3]   # land rent
	m.Lr   = x[4]   # employment in rural sector
	m.pr   = x[5]   # relative price rural good
	m.Sr   = x[6]   # amount of land used in rural production

	# update equations
	m.Lu   = p.L - m.Lr   # employment in urban sector
	m.wu0  = wu0(m.Lu,p)   # wage rate urban sector at city center (distance = 0)
	m.wr   = wr(m.Lu,m.ϕ,p) # wage rate rural sector
	m.xsr  = xsr(p,m)
	m.Srh  = Srh(p,m)
	m.qr   = qr(p,m)
	m.U    = utility(0.0,p,m)
	m.ϵr   = ϵ(m.ϕ,m.ϕ,p)
	m.nodes[:] .= m.ϕ / 2 .+ (m.ϕ / 2) .* m.inodes   # maps [-1,1] into [0,ϕ]

end


γ(l::Float64,ϕ::Float64,p::Param) = p.γ / (1.0 + ϵ(l,ϕ,p))

"commuting cost"
τ(x::Float64,ϕ::Float64,p::Param) = x > ϕ ? 0.0 : p.τ * x

"urban wage at location ``l``"
wu(Lu::Float64,l::Float64,ϕ::Float64,p::Param) = p.Ψ * p.θu * Lu^p.η * (1.0 .- τ(l,ϕ,p))

"urban wage at center"
wu0(Lu::Float64,p::Param) = p.Ψ * p.θu * Lu^p.η

"rural wage from indifference condition at ϕ. Eq (11)"
wr(Lu::Float64,ϕ::Float64,p::Param) = wu(Lu,ϕ,ϕ,p)

"wage at location ``l``"
w(Lu::Float64,l::Float64,ϕ::Float64,p::Param) = l >= ϕ ? wr(Lu,ϕ,p) : wu(Lu,l,ϕ,p)

"excess subsistence urban worker"
xsu(l::Float64,p::Param,m::Model) = w(m.Lu,l,m.ϕ,p) .+ m.r .- m.pr .* p.cbar .+ p.sbar

"excess subsistence rural worker"
xsr(p::Param,m::Model) = m.r + wr(m.Lu,m.ϕ,p) - m.pr * p.cbar + p.sbar

"Residential Land in Rural sector. equation (20)"
Srh(p::Param,m::Model) = m.Lr * (γ(m.ϕ,m.ϕ,p) * m.xsr) / ((1.0+ϵ(m.ϕ,m.ϕ,p)) * m.ρr)

"optimal urban good consumption at location ``l``. Equation (8)"
cu(l::Float64,p::Param,m::Model) = (1.0 - p.γ)*(1.0 - p.ν)*(w(m.Lu,l,m.ϕ,p) .+ (m.r - m.pr * p.cbar))

"optimal rural good consumption at location ``l``. Equation (7)"
cr(l::Float64,p::Param,m::Model) = (1.0 - p.γ) * p.ν * (w(m.Lu,l,m.ϕ,p) .+ (m.r - m.pr * p.cbar)) + m.pr * p.cbar

"urban good consumption in rural sector. Equation (8)"
cur(p::Param,m::Model) = cu(m.ϕ,p,m)

"cost of construction at location ``l``"
cfun(l::Float64,p::Param) = p.c0 + p.c1 * l + p.c2 * l^2
cost(l::Float64,ϕ::Float64,p::Param) = l >= ϕ ? cfun(ϕ,p) : cfun(l,p)

"housing supply shifter at ``l``"
χ(l::Float64,ϕ::Float64,p::Param) = (1.0 / cost(l,ϕ,p))^ϵ(l,ϕ,p)

"housing supply elasticity at ``l``"
ϵ(l::Float64,ϕ::Float64,p::Param) = l <= ϕ ? p.ϵ0 + p.ϵ1 * l + p.ϵ2 * l^2 : p.ϵ0 + p.ϵ1 * ϕ + p.ϵ2 * ϕ^2

"house price function at ``l``. equation (12)"
q(l::Float64,p::Param,m::Model) = m.qr * (xsu(l,p,m) / m.xsr).^(1.0/p.γ)

"land price function at ``l``. equation (15)"
# ρ(l::Float64,p::Param,m::Model) = H(l,p,m) .* q(l,p,m) / (1.0 + p.ϵ)
ρ(l::Float64,p::Param,m::Model) = (χ(l,p,m) .* q(l,p,m).^(1.0 + ϵ(l,ϕ,p))) / (1.0 + ϵ(l,ϕ,p))

"rural house price from land price"
qr(p::Param,m::Model) = ( (1+p.ϵ) * m.ρr * cfun(m.ϕ,p)^p.ϵ ) ^(1.0/(1+p.ϵ))



"housing demand at location ``l``"
h(l::Float64,p::Param,m::Model) = (p.γ / m.qr) * xsu(l,p,m)^((p.γ-1)/p.γ) * m.xsr^(1/p.γ)

"housing supply at location ``l``"
H(l::Float64,p::Param,m::Model) = χ(l,p,m) * q(l,p,m).^ϵ(l,ϕ,p)

"Population Density at location ``l``"
D(l::Float64,p::Param,m::Model) = H(l,p,m) / h(l,p,m)

"Amount of Urban Good (numeraire) required to build housing at ``l``"
cu_input(l::Float64,p::Param,m::Model) = ϵ(l,ϕ,p) / (1+ϵ(l,ϕ,p)) .* q(l,p,m) .* H(l,p,m)

"Utility function equation (5)"
utility(l::Float64,p::Param,m::Model) = (cr(l,p,m) - p.cbar)^(p.ν * (1-p.γ)) * (cu(l,p,m))^((1-p.ν) * (1-p.γ)) * h(l,p,m)^p.γ

"""
	Solves the General Case Model 
"""
function solve!(F,x,p::Param,m::Model)

	if any(x .< 0)
		F .= PEN
	else
		update!(m,p,x)
		if (xsu(0.0,p,m) < 0.0) 
		# if (xsu(0.0,p,m) < 0.0) || (xsr(p,m) < 0.0)
			# F .= PEN
			println("true")
		else
			integrate!(m,p)
			Eqsys!(F,m,p)
		end
		# @debug "StructChange! values:" qr=qr Lr=Lr wu wr Sr r F1=F[1] F2=F[2]
	end
	
end

"""
	Eqsys!(F::Vector{Float64},m::Model,p::Param)

compute system of equations for general commuting cost case.
"""
function Eqsys!(F::Vector{Float64},m::Model,p::Param)
	σ1 = (p.σ-1)/p.σ
	σ2 = 1.0 / (p.σ-1)

	# rural firm w-FOC: equation (2)
	F[1] = m.wr - p.α * m.pr * p.θr * (p.α + (1-p.α)*(m.Sr / m.Lr)^σ1)^σ2  #

	# rural firm q-FOC: equation (3)
	F[2] = m.ρr - (1-p.α)*m.pr * p.θr * (p.α * (m.Lr / m.Sr)^σ1 + (1-p.α))^σ2 

	# city size - Urban population relationship: equation (16)
	F[3] = m.Lu - m.iDensity

	#  total land rent per capita: equation before (19)
	F[4] = m.r * p.L - m.ir - m.ρr * (1 - m.ϕ)

	# land market clearing: after equation (18)  
	F[5] = 1.0 - p.λ - m.ϕ - m.Sr - m.Srh 

	# urban goods market clearing. equation before (22) but not in per capita terms
	#      rural cu cons + urban cu cons + rural constr input + urban constr input + commuting - total urban production
	F[6] = m.Lr * cur(p,m) + m.icu + m.Srh * cu_input(m.ϕ,p,m) + m.icu_input + m.iτ - wu0(m.Lu, p)*m.Lu
	
end

