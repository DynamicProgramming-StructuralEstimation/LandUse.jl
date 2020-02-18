

mutable struct Country
	R  :: Vector{Region}   # a set of regions
	wr :: Float64          # a global rural wage
	pr :: Float64          # a global price of rural good
	ρr :: Float64          # a global land price
	r  :: Float64          # country-wide per capita rental income
	LS  :: Float64          # Lr/Sr for region 1
	Sk :: Vector{Float64}  # total space for each region
	Xtrace ::Matrix{Float64}
	Ftrace ::Matrix{Float64}

	function Country(p::Vector{Param})
		this = new()
		this.R = [Region(pp) for pp in p]
		this.pr = NaN
		this.ρr = NaN
		this.r  = NaN
		this.LS  = NaN
		this.Sk  = p[1].kshare .* p[1].S
		this.Xtrace = zeros(3 + length(p),1)
		this.Ftrace = zeros(3 + length(p),1)
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
	c.ρr = foc_Sr( c.LS , c.pr, p1)
	# rural wage
	c.wr = foc_Lr( c.LS , c.pr, p1)

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

		# 3. update remaining fields
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

	# # 3. rural labor/land ratio is constant across K (and equal to region 1)
	# for ik in 2:K
	# 	fi += 1
	# 	F[fi] = (C.R[ik].Lr / C.R[ik].Sr) - (C.R[1].Lr / C.R[1].Sr)
	# end

	# 3. Aggregate land rents
	fi += 1
	# F[fi] = C.r * p[1].L - sum(m.iq + m.ρr * (m.Sk - m.ϕ) for m in C.R)
	F[fi] = C.r * p[1].L - sum(C.R[ik].iq + C.R[ik].ρr * (C.Sk[ik] - C.R[ik].ϕ) for ik in 1:K)

	# 4. aggregate urban consumption good clears
	urban_prod = sum(Yu(C.R[ik],p[ik]) for ik in 1:K)
	fi += 1
	m = C.R
	F[fi] = urban_prod - sum( m[ik].Lr * cur(p[ik],m[ik]) + m[ik].icu + m[ik].Srh * cu_input(m[ik].ϕ,p[ik],m[ik]) + m[ik].icu_input + m[ik].wu0 * m[ik].iτ for ik in 1:K)

	# K + 3 equations up to here
	C.Ftrace = hcat(C.Ftrace,F)

end

function solve!(F,x,p::Vector{Param},C::Country)
	# println("current x = $x")
		C.Xtrace = hcat(C.Xtrace,x)
	# push!(C.Xtrace,x)
	# println("last element in xtrace = $(C.Xtrace[end])")
	global count += 1
	# println("iteration $count")
	# println(x)
	if any(x .< 0)
		F[:] .= NaN
	else
		update!(C,p,x)
		EqSys!(F,C,p)
	end
end

function runk(;par = Dict(:S => 1.0, :L => 1.0, :kshare => [0.6,0.4], :K => 2),maxstep = 0.01)


	# 1. run a single region with pop = 1 and area = 1
	x,M,p = LandUse.run(par=par)
	K = p.K
	# 2. run a country with 2 regions, total pop 2 and total area =2. starting from solution in period 1 of 1.

    # p1 = LandUse.Param(par = par)  # double space and pop for 2 equally sized regions.
    # p2 = LandUse.Param(par = par)
	pp = [LandUse.Param(par = par) for ik in 1:K]
	C = Country[]
	sols = Vector{Float64}[]  # an empty array of vectors

	push!(C,LandUse.Country(pp))

	# starting values.
	# 1. b: ratio of labor to land in region 1
	# 2. r: land rent
	# 3. pr: relative price rural good
	# 4. 4 - K+3, Sr: amount of land use in rural production (complement of Srh)
	x0 = zeros(K + 3)
	x0[1] = M[1].Lr / M[1].Sr
	x0[2] = M[1].r
	x0[3] = M[1].pr
	for ik in 1:K
		x0[3 + ik] = M[1].Sr
	end

	r = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,pp,C[1]),x0,iterations = 100, show_trace=false,extended_trace=true)
	# r = LandUse.mcpsolve((F,x) -> LandUse.solve!(F,x,[p;p],C),zeros(length(x0)),fill(Inf,length(x0)),x0,iterations = 100, show_trace=true,reformulation = :minmax)



	if !converged(r)
		error("not converged")
	else
		LandUse.update!(C[1],pp,r.zero)
		push!(sols, r.zero)
	end
	Δθ = diff(pp[1].Θu)   # notice that is the set of all θus (it's a capital Θ !)
	for it in 2:length(pp[1].T)
	# for it in 2:7
		setperiod!(pp, it)   # set all params to it
		if Δθ[it-1] > maxstep
			# adaptive step
			C0,x = adapt_θ(sols[it-1],pp,Δθ[it-1],it,maxstep=maxstep)
			push!(sols,x)
			push!(C,C0)
		else
			# direct
			C0,x = step_country(sols[it-1],pp,it)
			push!(sols,x)
			push!(C,C0)
		end
	end
	sols, C

end

function adapt_θ(x0::Vector{Float64},p::Vector{Param},step::Float64,it::Int;maxstep = 0.04)

	# how many steps to take
	s = Int(fld(step,maxstep)) + 1

	# range of values to achieve next step
	θs = range(p[1].θu - step, stop = p[1].θu, length = s+1)[2:end]   # first one is done already
	sols = Vector{Float64}[]
	cs = Country[]
	push!(sols,x0)

	for (i,θ) in enumerate(θs)
		@debug "θ adapting" step=i θ=θ
		setfields!(p, :θu, θ)   # sets on each member of p
		setfields!(p, :θr, θ)
		c,x = step_country(sols[i],p,it)
		push!(sols,x)
		push!(cs,c)
	end

	return (cs[end],sols[end])    # return final step

end

function step_country(x0::Vector{Float64},pp::Vector{Param},it)
	C0 = LandUse.Country(pp)

	r1 = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,pp,C0),x0,iterations = 100, show_trace=false,extended_trace=true)
	if converged(r1)
		update!(C0,pp,r1.zero)
		traceplot(C0,it)
		return C0,r1.zero
	else
		traceplot(C0,it)
		error("Country not converged in period $it")
	end
end

function run7()
	par = Dict(:S => 1.0, :L => 2.0, :kshare => [0.6,0.1,0.1,0.1,0.1], :K => 5)
	x,C = LandUse.runk(par = par)
end
function run4()
	par = Dict(:S => 1.0, :L => 1.0, :kshare => [0.5,0.3,0.1,0.1], :K => 4)
	x,C = LandUse.runk(par = par)
end
