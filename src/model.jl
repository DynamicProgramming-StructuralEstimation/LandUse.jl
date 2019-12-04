"""
Model with general commuting cost.
"""
mutable struct Model
	qr      :: Float64   # land price in rural sector
	Lr      :: Float64   # employment in rural sector
	Lu      :: Float64   # employment in urban sector
	wu0     :: Float64   # wage in urban sector (at center)
	wr      :: Float64   # wage in rural sector
	Sr      :: Float64   # Amount of land used in rural production
	Srh     :: Float64   # Amount of land used for rural housing
	r       :: Float64   # per capita land rental income
	pr      :: Float64   # relative price of rural good
	ϕ       :: Float64   # size of the city
	ur      :: Float64  # net utility rural worker
	# integration setup
	inodes   :: Vector{Float64}  # theoretical integration nodes
	iweights :: Vector{Float64}  # int weights
	nodes   :: Vector{Float64}  # points where to evaluate integrand (inodes scaled into [0,ϕ])

	# resulting integrals from [0,ϕ]
	icu_input :: Float64   # ∫ cu_input(l) D(l) dl
	iDensity  :: Float64   # ∫ D(l) dl
	icu       :: Float64   # ∫ cu(l) D(l) dl
	iτ        :: Float64   # ∫ τ(l) D(l) dl
	ir        :: Float64   # ∫ qr(l) D(l) dl
	function Model(p::Param)
		# creates a model fill with NaN
		m     = new()
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
		m.ur  = NaN
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


"""
	integrate!(m::Model,p::Param)

could be made much faster by filling a matrix col-wise with integration nodes and 
doing a matmul on it? have to allocate memory though.
"""
function integrate!(m::Model,p::Param)
	
	m.icu_input = sum(m.iweights[i] * cu_input(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.iDensity  = sum(m.iweights[i] * D(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.icu       = sum(m.iweights[i] * cu(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
	m.iτ        = sum(m.iweights[i] * τ(m.nodes[i],m.ϕ,p) for i in 1:p.int_nodes)[1]
	m.ir        = sum(m.iweights[i] * qr(m.nodes[i],p,m) for i in 1:p.int_nodes)[1]
end

"""
	update!(m::Model,m0::CD0Model, p::Param)

update the general model from a CD0Model
"""
function update!(m::Model,m0::CD0Model,p::Param)
	m.qr   = m0.qr
	m.ϕ    = m0.ϕ 
	m.r    = m0.r 
	m.Lr   = m0.Lr
	m.pr   = m0.pr
	m.Sr   = m0.Sr

	# update equations
	m.Lu   = p.L - m.Lr   # employment in urban sector
	m.wu0  = wu0(m.Lu,p)   # wage rate urban sector at city center (distance = 0)
	m.wr   = wr(m.Lu,m.ϕ,p) # wage rate rural sector
	m.ur   = ur(p,m)
	m.Srh  = Srh(p,m)
	m.nodes[:] .= m.ϕ / 2 .+ (m.ϕ / 2) .* m.inodes   # maps [-1,1] into [0,ϕ]

end

"""
	update!(m::Model,p::Param,x::Vector{Float64})

update the general model from a parameter vector
"""
function update!(m::Model,p::Param,x::Vector{Float64})
	m.qr   = x[1]   # land price in rural sector
	m.ϕ    = x[2]   # city size
	m.r    = x[3]   # land rent
	m.Lr   = x[4]   # employment in rural sector
	m.pr   = x[5]   # relative price rural good
	m.Sr   = x[6]   # amount of land used in rural production

	# update equations
	m.Lu   = p.L - m.Lr   # employment in urban sector
	m.wu0  = wu0(m.Lu,p)   # wage rate urban sector at city center (distance = 0)
	m.wr   = wr(m.Lu,m.ϕ,p) # wage rate rural sector
	m.ur   = ur(p,m)
	m.Srh  = Srh(p,m)
	m.nodes[:] .= m.ϕ / 2 .+ (m.ϕ / 2) .* m.inodes   # maps [-1,1] into [0,ϕ]

end


"commuting cost"
τ(x::Float64,ϕ::Float64,p::Param) = x >= ϕ ? 0.0 : p.τ * x

"urban wage at location ``l``"
wu(Lu::Float64,l::Float64,ϕ::Float64,p::Param) = p.Ψ * p.θu * Lu^p.η * (1.0 .- τ(l,ϕ,p))

"urban wage at center"
wu0(Lu::Float64,p::Param) = p.Ψ * p.θu * Lu^p.η

"rural wage from indifference condition at ϕ. Eq (11)"
wr(Lu::Float64,ϕ::Float64,p::Param) = wu(Lu,ϕ,ϕ,p)

"net utility urban worker"
uu(l::Float64,p::Param,m::Model) = m.r - m.pr * p.cbar + p.sbar .+ wu(m.Lu,l,m.ϕ,p) 

"net utility rural worker"
ur(p::Param,m::Model) = m.r + wr(m.Lu,m.ϕ,p) - m.pr * p.cbar + p.sbar

"Residential Land in Rural sector"
Srh(p::Param,m::Model) = m.Lr * (p.γ * m.ur) / ((1.0+p.ϵ) * m.qr)

"urban good consumption at location ``l``. Equation (8)"
cu(l::Float64,p::Param,m::Model) = (1.0 - p.γ)*(1.0 - p.ν)*(wu(m.Lu,l,m.ϕ,p) .+ (m.r - m.pr * p.cbar))

"cost of construction at location ``l``"
cost(l::Float64,p::Param) = p.c0 + p.c1 * l + p.c2 * l^2


"house price function at ``l``"
q(l::Float64,p::Param,m::Model) = m.qr * (uu(l,p,m) / m.ur).^(1.0/p.γ)

"land price function at ``l``"
qr(l::Float64,p::Param,m::Model) = H(l,p,m) .* q(l,p,m) / (1.0 + p.ϵ)


"housing demand at location ``l``"
h(l::Float64,p::Param,m::Model) = l >= m.ϕ ? (p.γ / m.qr) * m.ur : (p.γ / m.qr) * m.ur^((p.γ-1)/p.γ) * m.ur^(1/p.γ)

"housing supply at location ``l``"
H(l::Float64,p::Param,m::Model) = cost(l,p).^(-p.ϵ) * q(l,p,m).^p.ϵ 

"Population Density at location ``l``"
D(l::Float64,p::Param,m::Model) = H(l,p,m) / h(l,p,m)

"Amount of Urban Good (numeraire) required to build housing at ``l``"
cu_input(l::Float64,p::Param,m::Model) = p.ϵ / (1+p.ϵ) .* q(l,p,m) .* H(l,p,m)

"""
	Solves the General Case Model 
"""
function solve!(F,x,p::Param,m::Model)

	if any(x .< 0)
		F .= PEN
	else
		update!(m,p,x)
		integrate!(m,p)
		# Eqsys!(F,m,p)
		# @debug "StructChange! values:" qr=qr Lr=Lr wu wr Sr r F1=F[1] F2=F[2]
	end
	
end

"""
	Eqsys!(F::Vector{Float64},m::Model,p::Param)

compute system of equations for general commuting cost case.
"""
function Eqsys!(F::Vector{Float64},m::Model,p::Param)
	σ1 = (p.σ-1)/p.σ
	σ2 = p.σ / (p.σ-1)
	σ_ = -1.0/p.σ
	y  = m.qr / p.L  + m.wu * (1.0 - m.τ * m.ϕ * m.Lr / p.L)

	F[1] = p.α * m.pr * p.θr * (p.α + (1 - p.α) * (m.Sr / m.Lr)^σ1)^σ2 - m.wr  # w-FOC of rural firm
	F[2] = (1 - p.α) * m.pr * p.θr * (p.α * (m.Sr / m.Lr)^σ1 + (1 - p.α))^σ2 - m.qr  # q-FOC of rural firm
	# Urban workers = density from 0 to ϕ. this is a closed form expression for this case
	F[3] = 1.0 + p.τ*m.wu*m.Lu/(p.χu*m.qr) - ((m.wu + m.r- m.pr*p.cbar+p.sbar)/(m.wu*(1-p.τ*m.ϕ)+m.r-m.p1*p.cbar+p.sbar))^(1/p.γ) 
	# definition of land rentes (equation 20). this is another closed form. 
	F[4] = m.qr * m.Sr + m.qr * p.χr * m.Srh + p.χu * p.γ * p.qr / (p.τ * (1 + p.γ) * m.wu) * 
		   (m.uu^(1/p.γ+1)/((m.wu*(1-p.τ*m.ϕ)+m.r-m.pr*p.cbar+p.sbar)^(1/p.γ))-(m.wu*(1-p.τ*m.ϕ)+m.r-m.p1*p.cbar+p.sbar)) - m.r * p.L
	F[5] = 1.0 - p.λ - m.ϕ - m.Srh - m.Sr   # land market clearing equation (18)
	# F[6] = 
end