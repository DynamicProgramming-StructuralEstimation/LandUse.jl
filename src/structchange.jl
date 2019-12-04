
# Cobb-Douglas Production Function in Rural Sector and no commuting cost

"""
Cobb-Douglas Production model without commuting cost.
"""
mutable struct CD0Model
	qr :: Float64   # land price in rural sector
	Lr :: Float64   # employment in rural sector
	Lu :: Float64   # employment in urban sector
	wu :: Float64   # wage in urban sector
	wr :: Float64   # wage in rural sector
	Sr :: Float64   # Amount of land used in rural production
	Srh:: Float64   # Amount of land used for rural housing
	r  :: Float64   # per capita land rental income
	pr :: Float64   # relative price of rural good
	ϕ  :: Float64   # size of the city
	uu :: Float64   # net utility of urban worker
	ur :: Float64   # net utility of rural worker
	function CD0Model(p::Param)
		m = new()
		m.qr = 0.3
		m.Lr = 0.2
		m.Lu   = p.L-m.Lr   # employment in urban sector
		m.wu   = p.θu     # wage rate urban sector with no commuting costs
		m.wr   = m.wu       # wage rate rural sector: equation (11)

		# amount of land used for r prod
		# equation (4) with σ=1
		m.Sr  = (1-p.α) / p.α * ( (m.wr * m.Lr) / m.qr )
		m.r   = 1/p.L*m.qr*(1-p.λ)  # per capita land rental income
		m.pr  = m.wr / (p.α * p.θr) * (m.Lr/m.Sr)^(1-p.α)  # rel price of rural conumption good. FOC of rural firm for Lr.
		m.ur = m.wr + m.r - m.pr * p.cbar
		m.uu = m.wu + m.r - m.pr * p.cbar
		m.ϕ    = (m.Lu/p.χu)*(p.γ * m.uu / m.qr )
		m.Srh  = (m.Lr/p.χr)*(p.γ * m.ur / m.qr )
		return m
	end
end

"""
	update!(m::CD0Model,p::Param,qr::Float64,Lr::Float64)

update a CD0 model
"""
function update!(m::CD0Model,p::Param,qr::Float64,Lr::Float64)
	m.qr = qr
	m.Lr = Lr
	m.Lu   = p.L-m.Lr   # employment in urban sector
	m.wu   = p.θu     # wage rate urban sector with no commuting costs
	m.wr   = m.wu       # wage rate rural sector: equation (11)

	# amount of land used for r prod
	# equation (4) with σ=1
	m.Sr  = (1-p.α) / p.α * ( (m.wr * m.Lr) / m.qr )
	m.r   = 1/p.L*m.qr*(1-p.λ)  # per capita land rental income
	m.pr  = m.wr / (p.α * p.θr) * (m.Lr/m.Sr)^(1-p.α)  # rel price of rural conumption good. FOC of rural firm for Lr.
	m.ur = m.wr + m.r - m.pr * p.cbar
	m.uu = m.wu + m.r - m.pr * p.cbar
	m.ϕ    = (m.Lu/p.χu)*(p.γ * m.uu / m.qr )
	m.Srh  = (m.Lr/p.χr)*(p.γ * m.ur / m.qr )

end

"""
	StructChange!(F,x,p)

solves a simple structural change model - no city structure.
in particular, no commuting cost, hence urban wage is ``θ_u`` everywhere.

This solves the model for ``σ = 1`` i.e. the cobb douglas case.
"""
function StructChange!(F,x,p::Param,m::CD0Model)

	qr   = x[1]   # land price in rural sector
	Lr   = x[2]   # employment in rural sector

	if (qr < 0) || (Lr < 0)
		F[1] = PEN
		F[2] = PEN
	else
		update!(m,p,qr,Lr)
		Eqsys!(F,m,p)

		# @debug "StructChange! values:" qr=qr Lr=Lr wu wr Sr r F1=F[1] F2=F[2]
	end
	
end

"""
	Eqsys!(F::Vector{Float64},m::CD0Model,p::Param)

compute system of equations for CD0 case.
"""
function Eqsys!(F::Vector{Float64},m::CD0Model,p::Param)
	F[1] = (1-p.γ)*(1-p.ν)* m.ur - p.θu*m.Lu
	F[2] = m.Sr + m.Srh + m.ϕ - (1-p.λ)
end



