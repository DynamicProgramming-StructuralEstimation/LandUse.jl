

mutable struct Country
	R  :: Vector{Region}   # a set of regions
	wr :: Float64          # a global rural wage
	pr :: Float64          # a global price of rural good
	ρr :: Float64          # a global land price
	r  :: Float64          # country-wide per capita rental income
	LS  :: Float64          # Lr/Sr for region 1
	Sk :: Vector{Float64}  # total space for each region
	L  :: Float64  # total population
	S  :: Float64 # total space

	function Country(cp::CParam,p::Vector{Param})
		this = new()
		this.R = [Region(pp) for pp in p]
		this.pr = NaN
		this.ρr = NaN
		this.r  = NaN
		this.LS  = NaN
		this.L  = cp.L
		this.S  = cp.S
		this.Sk  = cp.kshare .* cp.S
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

		# we know Sr in each region:
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
	F[fi] = C.L - sum(pop(i) for i in C.R)

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
	F[fi] = C.r * C.L - sum(C.R[ik].iq + C.R[ik].ρr * (C.Sk[ik] - C.R[ik].ϕ) for ik in 1:K)

	# 4. aggregate urban consumption good clears
	urban_prod = sum(Yu(C.R[ik],p[ik]) for ik in 1:K)
	fi += 1
	m = C.R
	F[fi] = urban_prod - sum( m[ik].Lr * cur(p[ik],m[ik]) + m[ik].icu + m[ik].Srh * cu_input(m[ik].ϕ,p[ik],m[ik]) + m[ik].icu_input + m[ik].wu0 * m[ik].iτ for ik in 1:K)

	# K + 3 equations up to here
	push!(Ftrace,F)

end

function solve!(F,x,p::Vector{Param},C::Country)
	push!(Xtrace,x)
	# println(x)
	if any(x .< 0)
		F[:] .= 10.0 .+ x.^2
	else
		update!(C,p,x)
		EqSys!(F,C,p)
	end
end

function NLopt_wrap(result::Vector, x::Vector, grad::Matrix,C::Country,p::Vector{Param})
	if length(grad) > 0
		# not implemented
	end
	solve!(result,x,p,C)
end

function nlopt_solveC(p::Vector{Param},C::Country,x0::Vector{Float64})
	K = length(C.R)
	nx = 3 + K
	opt = Opt(:LN_COBYLA,nx)
	# opt = Opt(:LN_BOBYQA,nx)
	# opt = Opt(:AUGLAG_EQ,nx)
	# opt = Opt(:GN_ISRES,nx)
	opt.lower_bounds = fill(0.001,nx)
	# opt.upper_bounds = fill(100.0,nx)
	# opt.upper_bounds = [1.0,1.0,1.0,1.0]
	f0(x::Vector, grad::Vector) = 1.0
	opt.min_objective = f0  # fix at a constant function
	equality_constraint!(opt,(r,x,g) -> NLopt_wrap(r,x,g,C,p), fill(1e-9,nx))
	# m.r    = x[1]   # land rent
	# m.Lr   = x[2]   # employment in rural sector
	# m.pr   = x[3]   # relative price rural good
	# m.Sr   = x[4]   # amount of land used in rural production
	# (optf,optx,ret) = optimize(opt, [0.1 * p.S, 0.5 * p.L, 0.3775, 0.545 * p.S])
	# if isnothing(x0)
	# 	x0 = [0.1 * p.S, 0.5 * p.L, 0.3775, 0.545 * p.S]   # an almost arbitrary point in the interior domain.
	# else
		@assert length(x0) == ndims(opt)
	# end
	(optf,optx,ret) = optimize(opt, x0)
end

function runk(;cpar = Dict(:S => 1.0, :L => 1.0, :kshare => [0.6,0.4], :K => 2),
				par = Dict())


	cp = LandUse.CParam(par = cpar)
	K = cp.K
	pp = convert(cp,par = par)  # create Lk and Sk for each region

	# 1. run single region model for each region
	Mk = [LandUse.run(pp[ik])[2][1] for ik in 1:K]   # [2] is model, [1] is first period

	# reset time
	setperiod!(pp,1)

	# 2. all regions in one country now. starting values for Sr from Mk.
	C = LandUse.Country(cp,pp)  # create that country

	# starting values.
	# 1. b: ratio of labor to land in region 1
	# 2. r: land rent
	# 3. pr: relative price rural good
	# 4. 4 - K+3, Sr: amount of land use in rural production (complement of Srh)
	x0 = zeros(K + 3)
	x0[1] = Mk[1].Lr / Mk[1].Sr  # take first countrys Lr/Sr ratio
	x0[2] = Mk[1].r
	x0[3] = Mk[1].pr
	for ik in 1:K
		x0[3 + ik] = Mk[1].Sr
	end

	r = LandUse.nlsolve((F,x) -> LandUse.solve!(F,x,pp,C),x0,iterations = 100, show_trace=false)
	# r = LandUse.nlopt_solveC(pp,C,x0)
	# r = LandUse.mcpsolve((F,x) -> LandUse.solve!(F,x,[p;p],C),zeros(length(x0)),fill(Inf,length(x0)),x0,iterations = 100, show_trace=true,reformulation = :minmax)

	# if (r[3] == :ROUNDOFF_LIMITED) | (r[3] == :SUCCESS)
	# 	LandUse.update!(C,pp,r[2])
	# 	println("eq system:")
	# 	LandUse.EqSys!(x0,C,pp)
	# 	println(x0)
	# else
	# 	error("not converged")
	# end
	#
	# return C,pp
	if converged(r)
		# push!(sols, r1.zero)
		update!(C,pp,r.zero)
		# @assert abs(Rmk(m0,p)) < 1e-7   # walras' law
		# push!(m,m0)
	else
		println(r)
		error("Country not converged")
	end
	C,r.zero

end
