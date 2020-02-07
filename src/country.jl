

mutable struct Country
	R  :: Vector{Region}   # a set of regions
	wr :: Float64          # a global rural wage
	pr :: Float64          # a global price of rural good
	ρr :: Float64          # a global land price
	r  :: Float64          # country-wide per capita rental income
	LS  :: Float64          # Lr/Sr for region 1
	Sk :: Vector{Float64}  # total space for each region

	function Country(p::Vector{Param})
		this = new()
		this.R = [Region(pp) for pp in p]
		this.pr = NaN
		this.ρr = NaN
		this.r  = NaN
		this.LS  = NaN
		this.Sk  = p[1].kshare .* p[1].S
		return this
	end
end

"""
	update!(c::Country,p::Vector{Param},x::Vector{Float64})

Update a `Country` object with a current vector of choice variables `x` supplied by the solver.
The ordering in `x` is:

1. b: ratio of labor to land in region 1
2. r: land rent
3. pr: relative price rural good
4. 4 - K+3, Sr: amount of land use in rural production (complement of Srh)
"""
function update!(c::Country,p::Vector{Param},x::Vector{Float64})
	# update country wide components
	c.LS   = x[1]   # constant labor/land share in region 1. i called that b>0 in the doc.
	c.r    = x[2]   # land rent
	c.pr   = x[3]   # relative price rural good
	Srk    = x[4:end]  # Sr for each region k

	K = length(p)
	@assert K == length(c.R)

	p1 = p[1]  # get first region's parameter to save typing

	# global implications
	σ1 = (p1.σ - 1) / p1.σ
	σ2 = 1 / (p1.σ - 1)
	# rural land price
	c.ρr = (1-p1.α) * c.pr * p1.θr * (p1.α * c.LS^σ1 + (1-p1.α) )^σ2
	# rural wage
	c.wr = c.ρr * p1.α / (1-p1.α) * c.LS^(-1/p1.σ)

	# 

	# update each region
	# 2. update other equations in each region
	for ik in 1:K
		# update country-wide stuff in each region
		c.R[ik].pr = c.pr 
		c.R[ik].r  = c.r 
		c.R[ik].ρr = c.ρr 
		c.R[ik].wr = c.wr # wage rate rural sector

		# we chose Sr in each region:
		c.R[ik].Sr = Srk[ik]
		# now we know Lr for each region:
		c.R[ik].Lr = c.LS * Srk[ik]

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

	# 3. rural labor/land ratio is constant across K (and equal to region 1)
	for ik in 2:K
		fi += 1
		F[fi] = (C.R[ik].Lr / C.R[ik].Sr) - (C.R[1].Lr / C.R[1].Sr)
	end

	# 4. Aggregate land rents
	fi += 1
	# F[fi] = C.r * p[1].L - sum(m.iq + m.ρr * (m.Sk - m.ϕ) for m in C.R)
	F[fi] = C.r * p[1].L - sum(C.R[ik].iq + C.R[ik].ρr * (C.Sk[ik] - C.R[ik].ϕ) for ik in 1:K)

	# 5. aggregate urban consumption good clears
	urban_prod = sum(Yu(C.R[ik],p[ik]) for ik in 1:K)
	fi += 1
	m = C.R
	F[fi] = urban_prod - sum( m[ik].Lr * cur(p[ik],m[ik]) + m[ik].icu + m[ik].Srh * cu_input(m[ik].ϕ,p[ik],m[ik]) + m[ik].icu_input + m[ik].wu0 * m[ik].iτ for ik in 1:K)

	# 2K + 2 equations up to here

	# add check of density in each city
	# for ik in 1:K
	# 	fi += 1
	# 	F[fi] = C.R[ik].Lu - C.R[ik].iDensity
	# end

	# # now 3K + 2 equations

end

function solve!(F,x,p::Vector{Param},C::Country)
	if any(x .< 0)
		F[:] .= NaN
	else
		update!(C,p,x)
		EqSys!(F,C,p)
	end
end

function runk()


	# 1. run a single region with pop = 1 and area = 1
	x,M,p = LandUse.run()
	# 2. run a country with 2 regions, total pop 2 and total area =2. starting from solution in period 1 of 1.

    p = LandUse.Param(par = Dict(:S => 2.0, :L => 2.0, :kshare => [0.5,0.5], :K => 2))  # double space and pop for 2 equally sized regions.
	C = LandUse.Country([p;p])  # create that country

	# starting values.
	x0 = zeros(2*p.K + 4)
	x0[1] = 


	x0 = 0.5*ones(2*2 + 4)
	r = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,[p;p],C),x0,iterations = 100, show_trace=false,extended_trace=true)
	# r = LandUse.mcpsolve((F,x) -> LandUse.solve!(F,x,[p;p],C),zeros(length(x0)),fill(Inf,length(x0)),x0,iterations = 100, show_trace=true,reformulation = :minmax)

end