
# models to generate starting value for full model

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
	get_starts()

Generate starting values for all years from the fixed ϵ model `FModel`.

In year 1:

1. Solve a model without commuting cost and 2 sectors. Choose land price and number of rural workers.
2. Construct a `FModel` from this solution.
3. Start with a very small city, ϕ = 0.00005, and solve with starting values from 2.

In subsequent years, take solution from 3. of previous year as starting value

"""
function get_starts()
	# 1. initialize parameter
	# p = Param(par=Dict(:ϵr => 0.0, :ϵs => 0.0))
	p = Param()
	fm = LandUse.FModel(p)  # create a fixed elasticity model
	startvals = Vector{Float64}[]  # an empty array of vectors

	# 2. For each time period
	for it in 1:length(p.T)
		setperiod!(p, it)   # set period on param to it

		# if first year
		if it == 1
			# 3. solve structural change model without space
			# --> use solution to that system, but replace ϕ with something very small.
			# --> produces starting value x00
			m0 = LandUse.CD0Model(p)
			r0 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,m0), [1; 0.5; 0.0])
			LandUse.update!(m0,p,r0.zero[1],r0.zero[2])
			# x00 = [m0.ρr; 0.00005; m0.r; m0.Lr; m0.pr; m0.Sr; m0.U]   # set very small city!
			x00 = [m0.ρr; 0.00005; m0.r; m0.Lr; m0.pr; m0.Sr]   # set very small city!

			# 4. solve general model with fixed elasticity starting from x00
			# --> closed form solutions for integrals
			# --> produces starting value x0
			# r1 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,fm),x00,iterations = 1000, autodiff = :forward)
			# r1 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,fm),x00,iterations = 100, method = :trust_region,show_trace = true, extended_trace = true)
			r1 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,fm),x00,iterations = 100)
			if converged(r1)
				push!(startvals, r1.zero)
				# println("rural market clears with $(Rmk(fm,p))")
			else
				error("first FModel not converged")
			end

		else  # in other years just start at previous solution
			lb = zeros(6)
			ub = fill(Inf,6)
			r1 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,fm),startvals[it-1],iterations = 100)
			# r1 = LandUse.mcpsolve((F,x) -> LandUse.solve!(F,x,p,fm),lb,ub, startvals[it-1],iterations = 10000, factor = 0.1)
			if converged(r1)
				push!(startvals, r1.zero)
				# println("rural market clears with $(Rmk(fm,p))")
			else
				print(r1)
				error("FModel not converged")
			end
		end
	end
	return startvals

end


"""
	adapt_ϵ(x0::Vector{Float64})

Adaptively increase slope coefficient ``s`` in elasticity of housing supply function [`ϵ`](@ref).
Starts from the first period solution of [`FModel`](@ref).
"""
function adapt_ϵ(x0::Vector{Float64})

	p = Param()
	m = GModel(p)

	startvals = Vector{Float64}[]  # an empty array of vectors
	push!(startvals, x0)  # put 1860 solution for flat epsilon function

	# range of elasticity slope values
	ϵs = range(0,stop = p.ϵsmax, length = p.ϵnsteps)[2:end]

	for (i,ϵ) in enumerate(ϵs)
		setfield!(p, :ϵs, ϵ)  # set current value for elaticity function slope

		r1 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,p,m),startvals[i],iterations = 1000)
		if converged(r1)
			push!(startvals, r1.zero)
		else
			error("adaptive search not converged for ϵ = $ϵ")
		end
	end
	return (startvals,p)
end
