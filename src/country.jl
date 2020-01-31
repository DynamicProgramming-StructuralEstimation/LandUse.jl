

mutable struct Country
	R  :: Vector{Region}   # a set of regions
	# p  :: Vector{Param}    # a set of region-specific parameters
	wr :: Float64          # a global rural wage
	pr :: Float64          # a global price of rural good
	ρr :: Float64          # a global land price
	r  :: Float64          # country-wide per capita rental income
	Sk :: Vector{Float64}  # share of total space for each region

	function Country(p::Vector{Param})
		this = new()
		this.R = [Region(pp) for pp in p]
		this.pr = NaN
		this.ρr = NaN
		this.r  = NaN
		this.Sk  = p[1].Sk
		return this
	end
end


function update!(c::Country,p::Vector{Param},x::Vector{Float64})
	# update country wide components
	c.ρr   = x[1]   # land price in rural sector
	c.wr   = x[2]   # rural wage
	c.r    = x[3]   # land rent
	c.pr   = x[4]   # relative price rural good


	K = length(p)
	@assert K == length(c.R)

	Lr = x[5:(K+4)]  # Lr choices
	Sr = x[(K+5):end]  # Lr choices

	# update each region
	# 2. update other equations in each region
	for ik in 1:K
		# update country-wide stuff in each region
		c.R[ik].pr = c.pr 
		c.R[ik].r  = c.r 
		c.R[ik].ρr = c.ρr 
		c.R[ik].wr = c.wr # wage rate rural sector

		# 1. set ϕ for each region
		c.R[ik].ϕ  = invτ(c.wr / p[ik].θu,p[ik])
		# 2. compute city size equation
		c.R[ik].nodes[:] .= c.R[ik].ϕ / 2 .+ (c.R[ik].ϕ / 2) .* c.R[ik].inodes   
		c.R[ik].Lu   = (c.R[ik].ϕ/2) * sum(c.R[ik].iweights[i] * D2(c.R[ik].nodes[i],p[ik],c.R[ik]) for i in 1:p[ik].int_nodes)[1]

		# 3. stick in rural workers and rural land
		c.R[ik].Lr   = Lr[ik]  
		c.R[ik].Sr   = Sr[ik]  
		c.R[ik].wu0  = wu0(c.R[ik].Lu,p[ik])   # wage rate urban sector at city center (distance = 0)
		# c.R[ik].wr   = wr(c.R[ik].Lu,c.R[ik].ϕ,p) # TEST
		c.R[ik].xsr  = xsr(p[ik],c.R[ik])
		c.R[ik].Srh  = Srh(p[ik],c.R[ik])
		c.R[ik].qr   = qr(p[ik],c.R[ik])
		integrate!(c.R[ik],p[ik])
	end
end

"""
Computes the entries of the residual vector ``u``
""" 
function EqSys!(F::Vector{Float64},C::Country,p::Vector{Param})

	K = length(C.R)  # num of regions
	fi = 1  # running index
	# 1. agg labor market clearing
	F[fi] = p[1].L - sum(pop(i) for i in C.R)

	# 2. land market clearing in each region K
	for ik in 1:K
		fi += 1
		F[fi] = C.Sk[ik] - C.R[ik].Sr - C.R[ik].ϕ - C.R[ik].Srh
	end

	# 3. rural labor/land ratio is constant across K
	for ik in 2:K
		fi += 1
		F[fi] = (C.R[ik].Lr / C.R[ik].Sr) - (C.R[1].Lr / C.R[1].Sr)
	end

	# 4. Aggregate land rents
	fi += 1
	F[fi] = C.r * p[1].L - sum(m.iq + m.ρr * (1 - m.ϕ) for m in C.R)

	# 5. aggregate urban consumption good clears
	urban_prod = sum(Yu(C.R[ik],p[ik]) for ik in 1:K)
	fi += 1
	m = C.R
	F[fi] = urban_prod - sum( m[ik].Lr * cur(p[ik],m[ik]) + m[ik].icu + m[ik].Srh * cu_input(m[ik].ϕ,p[ik],m[ik]) + m[ik].icu_input + m[ik].wu0 * m[ik].iτ for ik in 1:K)
	println(F)
end

function solve!(F,x,p::Vector{Param},C::Country)
	update!(C,p,x)
	EqSys!(F,C,p)
end

function runk()

    p = LandUse.Param(par = Dict(:ϵs => 0.0))
	LandUse.setperiod!(p,1)  # make sure we are in year 1
	C = LandUse.Country([p;p])
	x0 = 0.5*ones(2*2 + 4)
	r = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,[p;p],C),x0,iterations = 100, show_trace=true)

end