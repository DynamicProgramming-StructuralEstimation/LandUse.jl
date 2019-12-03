module LandUse
	using JSON
	using NLsolve

	const PEN = 10000.0

	using NLsolve

	mutable struct Param
		γ     :: Float64 # housing weight
		ϵ     :: Float64 # housing supply elasticity
		η     :: Float64 # agglomeration forces
		ν     :: Float64 # weight of rural good consumption on consumption composite
		cbar  :: Float64 # agr cons subsistence level (0.14 for thailand)
		sbar  :: Float64 # neg urba subsistence level
		Θr    :: Vector{Float64} # set of rural sector TFP
		Θu    :: Vector{Float64} # set of urban sector TFP
		θr    :: Float64 # rural sector TFP
		θu    :: Float64 # urban sector TFP
		α     :: Float64 # labor weight on farm sector production function
		λ     :: Float64 # useless land: non-farm, non-urban land (forests, national parks...)
		τ     :: Vector{Float64} # commuting cost function
		χr    :: Float64 # "easiness" to convert land into housing in rural sector
		χu    :: Float64 # "easiness" to convert land into housing in urban sector
		L     :: Float64 # total population (note: total land normalized to 1)
		ψ     :: Float64 # urban ammenities relative to rural ammenities (note: creates wedge in urban-rural wage)
		T     :: Vector{Float64}
		σ     :: Float64 # land-labor elasticity of substitution in farm production function

		function Param(;par=Dict())
	        f = open(joinpath(dirname(@__FILE__),"params.json"))
	        j = JSON.parse(f)
	        close(f)
	        this = new()
	        for (k,v) in j
	            if v["value"] isa Vector{Any}
	                setfield!(this,Symbol(k),convert(Vector{Float64},v["value"]))
	            else
	                setfield!(this,Symbol(k),v["value"])
	            end
	        end

	        if length(par) > 0
	            # override parameters from dict par
	            for (k,v) in par
	                setfield!(this,k,v)
	            end
	        end
        	return this
    	end
	end

	mutable struct Model
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
		function Model(p::Param)
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
		update!(m::Model,p::Param,qr::Float64,Lr::Float64)

	update a CD model
	"""
	function update!(m::Model,p::Param,qr::Float64,Lr::Float64)
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
	function StructChange!(F,x,p::Param,m::Model)

		qr   = x[1]   # land price in rural sector
		Lr   = x[2]   # employment in rural sector


		if (qr < 0) || (Lr < 0)
			F[1] = PEN
			F[2] = PEN
		else
			update!(m,p,qr,Lr)
			# Lu   = p.L-Lr   # employment in urban sector
			# wu   = p.θu     # wage rate urban sector with no commuting costs
			# wr   = wu       # wage rate rural sector: equation (11)

			# # amount of land used for r prod
			# # equation (4) with σ=1
			# Sr  = (1-p.α) / p.α * ( (wr * Lr) / qr )
			# r   = 1/p.L*qr*(1-p.λ)  # per capita land rental income
			# pr  = wr / (p.α * p.θr) * (Lr/Sr)^(1-p.α)  # rel price of rural conumption good. FOC of rural firm for Lr.
			# ur = wr + r - pr * p.cbar
			# uu = wu + r - pr * p.cbar
			# ϕ    = (Lu/p.χu)*(p.γ * uu / qr )
			# Srh  = (Lr/p.χr)*(p.γ * ur / qr )
			
			# F[1] = (1-p.γ)*(1-p.ν)*ur - p.θu*Lu
			# F[2] = Sr + Srh + ϕ - (1-p.λ)
			Fsys_CD!(F,m,p)

			# @debug "StructChange! values:" qr=qr Lr=Lr wu wr Sr r F1=F[1] F2=F[2]
		end
		
	end

	"""
		Fsys_CD!(F::Vector{Float64},m::Model,p::Param)

	compute system of equations for CD case.
	"""
	function Fsys_CD!(F::Vector{Float64},m::Model,p::Param)
		F[1] = (1-p.γ)*(1-p.ν)* m.ur - p.θu*m.Lu
		F[2] = m.Sr + m.Srh + m.ϕ - (1-p.λ)
	end


	function main()
		p = Param()
		m = Model(p)
		StructChange_closure(F,x) = StructChange!(F,x,p,m)
		r = nlsolve(StructChange_closure, [1; 0.5])
	end

	"""
		update(p::Param,qr::Float64,Lr::Float64)

	update a param from rural land price and rural employment
	"""
	function updateCD(p::Param,qr::Float64,Lr::Float64)
		Lu   = p.L-Lr   # employment in urban sector
		wu   = p.θu     # wage rate urban sector with no commuting costs
		wr   = wu       # wage rate rural sector: equation (11)

		# amount of land used for r prod
		# equation (4) with σ=1
		Sr  = (1-p.α) / p.α * ( (wr * Lr) / qr )
		r   = 1/p.L*qr*(1-p.λ)  # per capita land rental income
		pr  = wr / (p.α * p.θr) * (Lr/Sr)^(1-p.α)  # rel price of rural conumption good. FOC of rural firm for Lr.
		ur = wr + r - pr * p.cbar
		uu = wu + r - pr * p.cbar
		ϕ    = (Lu/p.χu)*(p.γ * uu / qr )

	end








end # module
