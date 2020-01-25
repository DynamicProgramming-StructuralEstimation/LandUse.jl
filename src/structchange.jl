
# Cobb-Douglas preferences models to generate starting value for full model

"""
Cobb-Douglas Production model without commuting cost.
"""
mutable struct CD0Model
	ρr :: Float64   # land price in rural sector
	χr :: Float64   # house supply shifter in rural sector
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
	U :: Float64   # Utility
	function CD0Model(p::Param)
		m = new()
		m.ρr = 0.3
		m.Lr = 0.2
		m.Lu   = p.L-m.Lr   # employment in urban sector
		m.wu   = p.θu     # wage rate urban sector with no commuting costs
		m.wr   = m.wu       # wage rate rural sector: equation (11)

		# amount of land used for r prod
		# equation (4) with σ=1
		m.Sr  = (1-p.α) / p.α * ( (m.wr * m.Lr) / m.ρr )
		m.r   = 1/p.L*m.ρr*(1-p.λ)  # per capita land rental income
		m.pr  = m.wr / (p.α * p.θr) * (m.Lr/m.Sr)^(1-p.α)  # rel price of rural conumption good. FOC of rural firm for Lr.
		m.ur = m.wr + m.r - m.pr * p.cbar
		m.uu = m.wu + m.r - m.pr * p.cbar
		m.χr  = χ(1.0,m.ϕ,p)
		m.ϕ    = (m.Lu/χ(0.0,1.0,p))*(p.γ * m.uu / m.ρr )
		m.Srh  = (m.Lr/m.χr)*(p.γ * m.ur / m.ρr )
		m.U = NaN
		return m
	end
end

"""
	update!(m::CD0Model,p::Param,ρr::Float64,Lr::Float64)

update a CD0 model
"""
function update!(m::CD0Model,p::Param,ρr::Float64,Lr::Float64)
	m.ρr = ρr
	m.Lr = Lr
	m.Lu   = p.L-m.Lr   # employment in urban sector
	m.wu   = p.θu     # wage rate urban sector with no commuting costs
	m.wr   = m.wu       # wage rate rural sector: equation (11)

	# amount of land used for r prod
	# equation (4) with σ=1
	m.Sr  = (1-p.α) / p.α * ( (m.wr * m.Lr) / m.ρr )
	m.r   = 1/p.L*m.ρr*(1-p.λ)  # per capita land rental income
	m.pr  = m.wr / (p.α * p.θr) * (m.Lr/m.Sr)^(1-p.α)  # rel price of rural conumption good. FOC of rural firm for Lr.
	m.ur = m.wr + m.r - m.pr * p.cbar
	m.uu = m.wu + m.r - m.pr * p.cbar
	m.χr = χ(1.0,m.ϕ,p)
	m.ϕ    = (m.Lu/χ(0.0,1.0,p))*(p.γ * m.uu / m.ρr )
	m.Srh  = (m.Lr/m.χr)*(p.γ * m.ur / m.ρr )

	# compute utility
	cr = (1.0 - p.γ) * p.ν * m.uu + m.pr * p.cbar
	cu = (1.0 - p.γ)*(1.0 - p.ν)*m.uu
	qr = ( (1+p.ϵr) * m.ρr * 1.0^p.ϵr ) ^(1.0/(1+p.ϵr))

	h = p.γ * (m.uu) / qr
	m.U = cr - p.cbar > 0.0 ? (cr - p.cbar)^(p.ν * (1-p.γ)) * (cu)^((1-p.ν) * (1-p.γ)) * h^p.γ : NaN

end

"""
	solve!(F,x,p::Param,m::CD0Model)

solves a simple structural change model - no city structure.
in particular, no commuting cost, hence urban wage is ``θ_u`` everywhere.

This solves the model for ``σ = 1`` i.e. the cobb douglas case.
"""
function solve!(F,x,p::Param,m::CD0Model)

	ρr   = x[1]   # land price in rural sector
	Lr   = x[2]   # employment in rural sector
	uconstraint   = x[3]   # well defined utility

	if (ρr < 0) || (Lr < 0)
		F[1] = PEN
		F[2] = PEN
	else
		update!(m,p,ρr,Lr)
		Eqsys!(F,m,p)

		# @debug "StructChange! values:" ρr=ρr Lr=Lr wu wr Sr r F1=F[1] F2=F[2]
	end

end

"""
	Eqsys!(F::Vector{Float64},m::CD0Model,p::Param)

compute system of equations for CD0 case.
"""
function Eqsys!(F::Vector{Float64},m::CD0Model,p::Param)
	F[1] = (1-p.γ)*(1-p.ν)* m.ur - p.θu*m.Lu
	F[2] = m.Sr + m.Srh + m.ϕ - (1-p.λ)
	cr = (1.0 - p.γ) * p.ν * m.uu + m.pr * p.cbar
	F[3] = cr > p.cbar ? 0.0 : (cr - p.cbar)^2
end


function CD()
	p = LandUse.Param()
	m0 = LandUse.CD0Model(p)
	CD0_closure(F,x) = LandUse.solve!(F,x,p,m0)
	r0 = LandUse.nlsolve(CD0_closure, [1; 0.5; 0.0])
	return (m0,p,r0)
end


"""
	update!(m::GModel,m0::CD0Model, p::Param)

update the general model from a CD0Model
"""
function update!(m::GModel,m0::CD0Model,p::Param)
	m.ρr   = m0.ρr
	m.ϕ    = m0.ϕ
	m.r    = m0.r
	m.Lr   = m0.Lr
	m.pr   = m0.pr
	m.Sr   = m0.Sr

	# update params


	# update equations
	m.Lu   = p.L - m.Lr   # employment in urban sector
	m.wu0  = wu0(m.Lu,p)   # wage rate urban sector at city center (distance = 0)
	m.wr   = wr(m.Lu,m.ϕ,p) # wage rate rural sector
	m.xsr  = xsr(p,m)
	m.qr   = qr(p,m)
	m.cr01 = (cr(0.0,p,m)-p.cbar, cr(1.0,p,m)-p.cbar)
	m.cu01 = (cu(0.0,p,m)       , cu(1.0,p,m)       )
	# if !all((m.cr01 .>= 0.0) .* (m.cu01 .>= 0.0) )
	# 	@warn "current model produces negative consumption"
	# end
	m.Srh  = Srh(p,m)
	m.U    = all((m.cr01 .>= 0.0) .* (m.cu01 .>= 0.0) ) ? utility(0.0,p,m) : NaN
	m.nodes[:] .= m.ϕ / 2 .+ (m.ϕ / 2) .* m.inodes   # maps [-1,1] into [0,ϕ]

end



"""
Model with linear commuting cost and fixed elasticities.
"""
mutable struct FModel <: Model
	ρr :: Float64   # land price in rural sector
	qr :: Float64   # house price in rural sector
	Lr :: Float64   # employment in rural sector
	Lu :: Float64   # employment in urban sector
	wu0 :: Float64   # wage in urban sector at 0 distance
	wr :: Float64   # wage in rural sector
	Sr :: Float64   # Amount of land used in rural production
	Srh:: Float64   # Amount of land used for rural housing
	r  :: Float64   # per capita land rental income
	pr :: Float64   # relative price of rural good
	ϕ  :: Float64   # size of the city
	xsu :: Float64   # net utility of urban worker
	xsr :: Float64   # net utility of rural worker
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
		m.xsu    = NaN
		m.xsr    = NaN
		return m
	end
end

"""
	update!(m::FModel,p::Param,x::Vector{Float64})

update a `Fmodel` from a choice vector x
"""
function update!(m::FModel,p::Param,x::Vector{Float64})
	m.ρr   = x[1]   # land price in rural sector
	m.ϕ    = x[2]   # city size
	m.r    = x[3]   # land rent
	m.Lr   = x[4]   # employment in rural sector
	m.pr   = x[5]   # relative price rural good
	m.Sr   = x[6]   # amount of land used in rural production

	# update params
	γ2 = p.γ / (1+p.ϵr)

	# update equations
	m.Lu   = p.L - m.Lr   # employment in urban sector
	m.wu0  = wu0(m.Lu,p)   # wage rate urban sector at city center (distance = 0)
	m.wr   = wr(m.Lu,m.ϕ,p) # wage rate rural sector
	m.xsr = m.wr + m.r - m.pr * p.cbar - p.sbar
	m.xsu = m.wu0 + m.r - m.pr * p.cbar - p.sbar
	m.Srh  = ( γ2 / m.ρr ) * m.xsr * m.Lr
	m.qr   = qr(p,m)

end

"""
	Solves the `FModel`
"""
function solve!(F,x,p::Param,m::FModel)
# println(x)
	# if any(x .< 0)
	# 	F .= PEN
	# 	println("penalizing neg inputs")
	# # elseif m.xsr < 0
	# # 	# m.r + wr(m.Lu,m.ϕ,p) - m.pr * p.cbar + p.sbar < 0
	# # 	# that means that the price of rural good is chosen too high.
	# # 	# imply that this way overshoots the first two equations
	# # 	# F[1] = -PEN
	# # 	# F[2] = -PEN
	# # 	F .= PEN
	# # 	println("penalizing neg cons")
	# else
		# @debug "l=0" cr=cr(0.0,p,m)-p.cbar cu=cu(0.0,p,m) h=h(0.0,p,m)
		# @debug "l=1" cr=cr(1.0,p,m)-p.cbar cu=cu(1.0,p,m) h=h(1.0,p,m)
		update!(m,p,x)

		# if negative consumption in either sector
		# if isnan(m.U)
		# 	F .= PEN
		# else
			Eqsys!(F,m,p)
		# end
	# end

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

	# rural firm w-FOC: equation (2)
	F[1] = p.α * m.pr * p.θr * (p.α + (1-p.α)*(m.Sr / m.Lr)^σ1)^σ2 - m.wr
	# F[1] = ((m.Sr > 0) && (m.Lr > 0)) ? p.α * m.pr * p.θr * (p.α + (1-p.α)*(m.Sr / m.Lr)^σ1)^σ2 - m.wr : 2*exp(min(m.Sr,m.Lr))
	# if !((m.Sr > 0) && (m.Lr > 0))
	# 	println(min(m.Sr,m.Lr))
	# end
	# rural firm q-FOC: equation (3)
	F[2] = (1-p.α)*m.pr * p.θr * (p.α * (m.Lr / m.Sr)^σ1 + (1-p.α))^σ2 - m.ρr

	# city size - Urban population relationship: analytic integral solution to equation (16)
	w2 = p.θu * m.Lu^p.η
	τϕ = 1-p.τ * m.ϕ
	w1 = p.Ψ * w2 * τϕ
	uu = w2+m.r - m.pr * p.cbar - p.sbar
	ur = w1+m.r - m.pr * p.cbar - p.sbar
	xx = m.r - m.pr * p.cbar - p.sbar


	F[3] = 1+p.τ * w2 * m.Lu/(m.ρr)-(uu/(w2*(1-p.τ*m.ϕ)+m.r-m.pr*p.cbar+p.sbar))^(1/γ2)

	#  total land rent per capita: closed form solution
	F[4] = m.ρr * m.Sr + m.ρr * m.Srh +
	       γ2 * m.ρr/(p.τ*(1+γ2)*w2)*(uu^(1/γ2+1)/((w2*τϕ+m.r-m.pr * p.cbar+p.sbar)^(1/γ2)) -
		   (w2*(1-p.τ * m.ϕ)+m.r-m.pr * p.cbar+p.sbar)) - m.r*p.L

	# land market clearing: after equation (18)
	F[5] = 1.0 - p.λ - m.ϕ - m.Sr - m.Srh

	# urban goods market clearing.
	chi2 = 1.0

	F[6] = m.Lr * (1-p.γ)*(1-p.ν)* ur +
	chi2*(1-p.γ)*(1-p.ν)*m.ρr/((1+γ2)*w2*p.τ)*((w2+xx)^(1/γ2+1)/((w2*τϕ+xx)^(1/γ2))-(w2*τϕ+xx)) +
	chi2*γ2*m.ρr/((1+γ2)*p.τ*w2)*(((w2+xx)^(1/γ2+1))/((w2*τϕ+xx)^(1/γ2))-(w2*τϕ+xx)) -
	m.ϕ*chi2*m.ρr +
	p.ϵr *m.ρr*m.Srh +
	p.ϵr *chi2*γ2*m.ρr/(p.τ*(1+γ2)*w2)*((w2+xx)^(1/γ2+1)/((w2*τϕ+xx)^(1/γ2))-(w2*τϕ+xx)) -
	p.sbar*p.L -
	p.θu*m.Lu^(1+p.η)

	cr = (1.0 - p.γ) * p.ν * ur + m.pr * p.cbar
	# println("cr - p.cbar = $(cr - p.cbar * p.cbar)")
	# F[7] = minerr < 0.0 ? minerr^2 : 0.0
	F[7] = cr - p.cbar * p.cbar > 0 ? 0.0 : 2*exp(cr - p.cbar * p.cbar)
	# F[7] = 0.0
	# println("penalty = $(F[7])")

	# F[7] = 0.0
end
