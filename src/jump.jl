


"""
solve model at current paramter p and starting at point x0
"""
function jm(p::LandUse.Param,mo::LandUse.Region,x0::NamedTuple; estimateθ = false)

	# setup Model object
	m = JuMP.Model(Ipopt.Optimizer)
	set_optimizer_attribute(m, MOI.Silent(), true)
	# lbs = [x0...] .* 0.3

	# variables
	@variable(m, ρr >=         0.01 * x0.ρr  , start = x0.ρr)
	@variable(m, 0.01 * x0.ϕ <= ϕ <= p.S     , start = x0.ϕ )
	@variable(m, r >=      0.1 * x0.r       , start = x0.r )
	@variable(m, 0.1 * x0.Lr <= Lr <= p.L   , start = x0.Lr)
	@variable(m, pr >=        0.05 * x0.pr  , start = x0.pr)
	@variable(m, 0.05 * x0.Sr <= Sr <= p.S   , start = x0.Sr)
	#
	# @variable(m, ρr       , start = x0.ρr)
	# @variable(m,  ϕ   , start = x0.ϕ )
	# @variable(m, r       , start = x0.r )
	# @variable(m, Lr  , start = x0.Lr)
	# @variable(m, pr          , start = x0.pr)
	# @variable(m,  Sr   , start = x0.Sr)
	if estimateθ
		@variable(m, θu , start = x0.θu)
		@variable(m, θr , start = x0.θr)
	else
		θu = p.θu
		θr = p.θr
	end

	# nonlinear expressions

	@NLexpression(m, wu0, θu * (p.L - Lr)^p.η)
	@NLexpression(m, wr , p.Ψ * (wu0 - p.a * (wu0^(p.ξw)) * (ϕ^(p.ξl))) )
	@NLexpression(m, qr , ((1+p.ϵr) * ρr)^(1.0/(1+p.ϵr)) )
	@NLexpression(m, r_pr_csbar, r - pr * p.cbar + p.sbar )
	@NLexpression(m, xsr, wr + r_pr_csbar )
	@NLexpression(m, hr, p.γ * xsr / qr )
	@NLexpression(m, Hr, qr^p.ϵr )
	# @NLexpression(m, Srh, Lr * (xsr * p.γ / (1.0 + p.ϵr)) / ρr )
	@NLexpression(m, Srh, Lr * hr / Hr )
	@NLexpression(m, cur, (1.0 - p.γ)*(1.0 - p.ν)*(wr + r_pr_csbar) - p.sbar)
	@NLexpression(m, cu_inputr,  (qr^(1 + p.ϵr)) * p.ϵr / (1.0+p.ϵr) )

	# expressions indexed at location l
	@NLexpression(m, nodes[i = 1:p.int_nodes], ϕ / 2 + ϕ / 2 * mo.inodes[i] )
	if p.ϵflat
		@NLexpression(m, ϵ[i = 1:p.int_nodes], p.ϵr)
	else
		@NLexpression(m, ϵ[i = 1:p.int_nodes], p.ϵr * nodes[i] / ϕ + p.ϵs * (ϕ - nodes[i])/ϕ)
	end
	@NLexpression(m, τ[i = 1:p.int_nodes], p.a * wu0^(p.ξw) * nodes[i]^(p.ξl) )
	@NLexpression(m, w[i = 1:p.int_nodes], wu0 - τ[i] )
	# @warn "hard coding abs() for q function" maxlog=1
	# @NLexpression(m, q[i = 1:p.int_nodes], qr * (abs((w[i] + r_pr_csbar) / xsr))^(1.0/p.γ))
	@NLexpression(m, q[i = 1:p.int_nodes], qr * ((w[i] + r_pr_csbar) / xsr)^(1.0/p.γ))
	@NLexpression(m, H[i = 1:p.int_nodes], q[i]^ϵ[i])
	@NLexpression(m, h[i = 1:p.int_nodes], p.γ * (w[i] + r_pr_csbar) / q[i])
	@NLexpression(m, ρ[i = 1:p.int_nodes], (q[i]^(1.0 + ϵ[i])) / (1.0 + ϵ[i]) )
	@NLexpression(m, cu[i = 1:p.int_nodes], (1.0 - p.γ)*(1.0 - p.ν)*(w[i] + r_pr_csbar) - p.sbar)
	@NLexpression(m, D[i = 1:p.int_nodes] , H[i] / h[i])
	@NLexpression(m, cu_input[i = 1:p.int_nodes], q[i] * H[i] * ϵ[i] / (1.0+ϵ[i]) )



	# integrals
	@NLexpression(m, iDensity,  (ϕ/2) * sum(mo.iweights[i] * 2π * nodes[i] * D[i] for i in 1:p.int_nodes))
	@NLexpression(m, iρ,        (ϕ/2) * sum(mo.iweights[i] * 2π * nodes[i] * ρ[i] for i in 1:p.int_nodes))
	@NLexpression(m, icu,       (ϕ/2) * sum(mo.iweights[i] * 2π * nodes[i] * cu[i] * D[i] for i in 1:p.int_nodes))
	@NLexpression(m, icu_input, (ϕ/2) * sum(mo.iweights[i] * 2π * nodes[i] * cu_input[i] for i in 1:p.int_nodes))
	@NLexpression(m, iτ,        (ϕ/2) * sum(mo.iweights[i] * 2π * nodes[i] * (wu0 - w[i]) * D[i] for i in 1:p.int_nodes))

	# objective
	# if estimateθ
		# if p.it == 1
			# @objective(m, Min, (Lr / p.Lt[p.it] - p.moments[p.it,:Employment_rural])^2)
		# else
			@objective(m, Min, (pr - p.moments[p.it,:P_rural])^2)
			# @objective(m, Min, (Lr / p.Lt[p.it] - p.moments[p.it,:Employment_rural])^2 + (pr - p.moments[p.it,:P_rural])^2)
		# end
	# end

	# nonlinear constraints (they are actually linear but contain nonlinear expressions - which means we need the nonlinear setup)
	# F[1] = m.wr - foc_Lr(m.Lr / m.Sr , m.pr, p)
	@NLconstraint(m, wr - p.α * pr * θr * (p.α + (1-p.α)*( Sr / Lr )^((p.σ-1)/p.σ))^(1.0 / (p.σ-1)) == 0.0)

	# F[2] = m.ρr - foc_Sr(m.Lr / m.Sr , m.pr, p)
	@NLconstraint(m, ρr - (1-p.α)* pr * θr * (p.α * ( Lr / Sr )^((p.σ-1)/p.σ) + (1-p.α))^(1.0 / (p.σ-1)) == 0.0)

	# F[3] = m.Lu - m.iDensity
	@NLconstraint(m, p.L - Lr == iDensity)

	# F[4] = m.iρ + m.ρr * (m.Sr + m.Srh) - m.r * p.L
	@NLconstraint(m, iρ + ρr * (Sr + Srh) == r * p.L)

	# F[5] = p.S - p.λ - m.ϕ - m.Sr - m.Srh
	@NLconstraint(m, p.S - p.λ == ϕ^2 * π + Sr + Srh)

	# F[6] = m.Lr * cur(p,m) + m.icu + m.Srh * cu_input(m.ϕ,p,m) + m.icu_input + m.wu0 * m.iτ - wu0(m.Lu, p)*m.Lu
	@NLconstraint(m, Lr * cur + icu + Srh * cu_inputr + icu_input + iτ ==  wu0 * (p.L - Lr))

	@NLconstraint(m, cur >= 0)
	# @NLconstraint(m, cui[i = 1:p.int_nodes], cu[i] >= 0)
	# @NLconstraint(m, xsr >= 0)


	if p.it == 1 && estimateθ
		@constraint(m, θu == θr)
	end

	JuMP.optimize!(m)
		# println("qr = $(value(qr))")
		# println("q = $(value(q[1]))")
		# println("wr         = $(value(wr        ))")
		# println("qr         = $(value(qr        ))")
		# println("r_pr_csbar = $(value(r_pr_csbar))")
		# println("xsr        = $(value(xsr       ))")
		# println("hr         = $(value(hr        ))")

	# println("rhor = $(value(ρr))")
	# println("ϕ    = $(value(ϕ ))")
	# println("r    = $(value(r ))")
	# println("Lr   = $(value(Lr))")
	# println("pr   = $(value(pr))")
	# println("Sr   = $(value(Sr))")
	# println("θr   = $(value(θr))")
	# println("θu   = $(value(θu))")

	# check termination status
	if termination_status(m) != MOI.LOCALLY_SOLVED
		# println("error in period $(p.it)")
		# println("Termination status: $(termination_status(m))")
		# println("rhor = $(value(ρr))")
		# println("ϕ    = $(value(ϕ ))")
		# println("r    = $(value(r ))")
		# println("Lr   = $(value(Lr))")
		# println("pr   = $(value(pr))")
		# println("Sr   = $(value(Sr))")
		# println("θr   = $(value(θr))")
		# println("θu   = $(value(θu))")
		error("model not locally solved")
	else
		# if p.it == 1
			# println("period $(p.it) solution")
			#
			# println("rhor = $(value(ρr))")
			# println("ϕ    = $(value(ϕ ))")
			# println("r    = $(value(r ))")
			# println("Lr   = $(value(Lr))")
			# println("pr   = $(value(pr))")
			# println("Sr   = $(value(Sr))")
			# println("θr   = $(value(θr))")
			# println("θu   = $(value(θu))")
		# end


		(ρr = value(ρr), ϕ = value(ϕ), r = value(r), Lr = value(Lr), pr = value(pr), Sr = value(Sr), θu = estimateθ ? value(θu) : θu, θr = estimateθ ? value(θr) : θr)
	end
end


"""
solve a `Country` at point x0
"""
function jc(C::Country,x0::Vector)

	pp = C.pp   # vector of params
	p  = pp[1]   # param of region 1
	K = C.K

	# setup Model object
	m = JuMP.Model(Ipopt.Optimizer)
	set_optimizer_attribute(m, MOI.Silent(), true)
	# lbs = [x0...] .* 0.3

	# variables
	@variable(m, LS >=  0.1 * x0[1]  , start = x0[1])
	@variable(m, r  >=   0.1 * x0[2]  , start = x0[2])
	@variable(m, pr >=  0.05 * x0[3]  , start = x0[3])
	@variable(m, 0.1 * x0[3 + ik]   <= Sr[ik = 1:K] <= C.Sk[ik], start = x0[3 + ik])
	@variable(m, 0.05 * x0[3+K + ik] <= Lu[ik = 1:K] <= C.L     , start = x0[3+K + ik])
	# @variable(m, minimum([pp[ik].θu for ik in 1:K]) >= wr >= 0.0)

	# all θs are in p[ik].θ
	σ1 = (p.σ - 1) / p.σ
	σ2 = 1 / (p.σ - 1)

	# most rural sector expressions are identical in all regions
	# rural land price

	# rural pop


	# constant expressions in each country's rural part
	@NLexpression(m, ρr , (1-p.α)* pr * p.θr * (p.α * (LS)^(σ1) + (1-p.α))^σ2)
	@NLexpression(m, qr , ((1+p.ϵr) * ρr)^(1.0/(1+p.ϵr)) )
	@NLexpression(m, wr , p.α * pr * p.θr * (p.α + (1-p.α)*( 1.0 / LS )^(σ1))^(σ2) )
	@NLexpression(m, r_pr_csbar, r - pr * p.cbar + p.sbar )
	@NLexpression(m, xsr, wr + r_pr_csbar )
	@NLexpression(m, hr, p.γ * xsr / qr )
	@NLexpression(m, Hr, qr^p.ϵr )
	@NLexpression(m, cur, (1.0 - p.γ)*(1.0 - p.ν)*(wr + r_pr_csbar) - p.sbar)
	@NLexpression(m, cu_inputr,  (qr^(1 + p.ϵr)) * p.ϵr / (1.0+p.ϵr) )

	# indexed by only k
	@NLexpression(m, Lr[ik = 1:K], LS * Sr[ik])  # rural pop from labor to land share LS
	@NLexpression(m, Srh[ik = 1:K], Lr[ik] * hr / Hr )   # housing space for rural pop
	@NLexpression(m, wu0[ik = 1:K], pp[ik].θu * Lu[ik]^p.η)  # urban wage in each city

	@NLexpression(m, ϕ[ik = 1:K], ( (wu0[ik] - wr) / (p.a * wu0[ik]^(p.ξw)) )^(1.0/p.ξl))  # fringe for each region from inverse moving cost function


	# expressions indexed at location l in each k
	@NLexpression(m, nodes[i = 1:p.int_nodes, ik = 1:K], ϕ[ik] / 2 + ϕ[ik] / 2 * p.inodes[i] )
	@NLexpression(m, ϵ[i = 1:p.int_nodes, ik = 1:K], p.ϵr * exp(-p.ϵs * (ϕ[ik]-nodes[i,ik])))
	@NLexpression(m, τ[i = 1:p.int_nodes,ik = 1:K], p.a * pp[ik].θu^(p.ξw) * nodes[i,ik]^(p.ξl) )
	@NLexpression(m, w[i = 1:p.int_nodes,ik = 1:K], pp[ik].θu - τ[i,ik] )
	# @warn "hard coding abs() for q function" maxlog=1
	# @NLexpression(m, q[i = 1:p.int_nodes], qr * (abs((w[i] + r_pr_csbar) / xsr))^(1.0/p.γ))
	@NLexpression(m,        q[i = 1:p.int_nodes,ik = 1:K], qr * ((w[i,ik] + r_pr_csbar) / xsr)^(1.0/p.γ))
	@NLexpression(m,        H[i = 1:p.int_nodes,ik = 1:K], q[i,ik]^ϵ[i,ik])
	@NLexpression(m,        h[i = 1:p.int_nodes,ik = 1:K], p.γ * (w[i,ik] + r_pr_csbar) / q[i,ik])
	@NLexpression(m,        ρ[i = 1:p.int_nodes,ik = 1:K], (q[i,ik]^(1.0 + ϵ[i,ik])) / (1.0 + ϵ[i,ik]) )
	@NLexpression(m,       cu[i = 1:p.int_nodes,ik = 1:K], (1.0 - p.γ)*(1.0 - p.ν)*(w[i,ik] + r_pr_csbar) - p.sbar)
	@NLexpression(m,        D[i = 1:p.int_nodes,ik = 1:K] , H[i,ik] / h[i,ik])
	@NLexpression(m, cu_input[i = 1:p.int_nodes,ik = 1:K], q[i,ik] * H[i,ik] * ϵ[i,ik] / (1.0+ϵ[i,ik]) )

	# integrals for each region ik
	@NLexpression(m, iDensity[ik = 1:K],  (ϕ[ik]/2) * sum(p.iweights[i] * 2π * nodes[i,ik] * D[i,ik] for i in 1:p.int_nodes))
	@NLexpression(m, iρ[ik = 1:K],        (ϕ[ik]/2) * sum(p.iweights[i] * 2π * nodes[i,ik] * ρ[i,ik] for i in 1:p.int_nodes))
	@NLexpression(m, icu[ik = 1:K],       (ϕ[ik]/2) * sum(p.iweights[i] * 2π * nodes[i,ik] * cu[i,ik] * D[i,ik] for i in 1:p.int_nodes))
	@NLexpression(m, icu_input[ik = 1:K], (ϕ[ik]/2) * sum(p.iweights[i] * 2π * nodes[i,ik] * cu_input[i,ik] for i in 1:p.int_nodes))
	@NLexpression(m, iτ[ik = 1:K],        (ϕ[ik]/2) * sum(p.iweights[i] * 2π * nodes[i,ik] * (pp[ik].θu - w[i,ik]) * D[i,ik] for i in 1:p.int_nodes))

	# constraints
	# @NLconstraint(m, aux_con_wr, wr == p.α * pr * p.θr * (p.α + (1-p.α)*( 1.0 / LS )^(σ1))^(σ2) )
	# @NLconstraint(m, aux_con_Srh[ik = 1:K], Srh[ik] >= 0.001 )
	@NLconstraint(m, C.L == sum(Lr[ik] + Lu[ik] for ik in 1:K))   # agg labor market clearing
	@NLconstraint(m, land_clearing[ik = 1:K], C.Sk[ik] == Sr[ik] + ϕ[ik]^2 * π + Srh[ik])   # land market clearing in each region
	@NLconstraint(m, r * C.L == sum(iρ[ik] + ρr * (Sr[ik] + Srh[ik]) for ik in 1:K))   # agg land rent definition
	# input of urban good equal urban good production
	@NLconstraint(m, sum(Lr[ik] * cur + icu[ik] + Srh[ik] * cu_inputr + icu_input[ik] + iτ[ik] for ik in 1:K) == sum(pp[ik].θu * Lu[ik] for ik in 1:K))
	@NLconstraint(m, city_size[ik = 1:K], Lu[ik] == iDensity[ik])

	# objective function
	@objective(m, Min, 1.0)  # constant function
	# @NLexpression(m, emp_share, )
	# @NLobjective(m, Min, ((sum(Lr[ik] for ik in 1:K) / C.L)  - p.moments[1,:Employment_rural])^2)



	JuMP.optimize!(m)

	if termination_status(m) == MOI.LOCALLY_SOLVED
		out = zeros(3 + 2K)
		out[1] = value(LS)
		out[2] = value(r)
		out[3] = value(pr)
		for ik in 1:K
			out[3 + ik] = value(Sr[ik])
			out[3 + K + ik] = value(Lu[ik])
		end
		return (out,value.(ϕ))

		# 		names = [:LS,:r,:pr, 
		# 		[Symbol("Sr_$ik") for ik in 1:K]..., 
		# 		[Symbol("Lu_$ik") for ik in 1:K]...,
		# 		[Symbol("ϕ_$ik") for ik in 1:K]...]
		# values = [value(LS),value(r), value(pr),
		# 		 [value(Sr[ik])  for ik in 1:K]..., 
		# 		 [value(Lu[ik])  for ik in 1:K]...,
		# 		 [value(ϕ[ik]) for ik in 1:K]...]

		# (; zip(names, values)...)  # return as named tuple
	else
		println(termination_status(m))
		for ik in 1:K
			println("k = $ik")
			println("pp[ik].θu=$(pp[ik].θu),xx=$(pp[ik].θu - value(wr)),wr=$(value(wr)),  ϕ = $(value(ϕ[ik])),iρ = $(value(iρ[ik])),ρr=$(value(ρr)),Sr=$(value(Sr[ik])),Srh = $(value(Srh[ik]))")
		end
	    error("The model was not solved correctly.")
	end
end

function jjc()
	m = runm()
	x0 = m[2][1]
	C = country()
	x = Float64[]
	push!(x, x0.Lr / x0.Sr)
	push!(x, x0.r)
	push!(x, x0.pr)
	for ik in 1:C.K
		push!(x,x0.Sr)
	end
	for ik in 1:C.K
		push!(x,x0.Lu)
	end
	jc(C,x)
end


function stmodel2(p::LandUse.Param)
	x00 = nocommute(p)
	# x00 = nlsolve((F,x) -> nocommute!(F,x,p),x0)
	gamma2 = p.γ / (1 + p.ϵr)

	ρ  = x00.ρ
	Lr = x00.Lr
	Lu = p.L - Lr
	wu = p.θu*Lu^p.η
	wr = wu
	Sr = (((1 - p.α)/ p.α) * wr / ρ)^p.σ * Lr # farm land input
	r  = ρ * (p.S - p.λ) / p.L
	pr = wr / (p.α * p.θr) * (p.α + (1-p.α) * (Lr / Sr)^((1-p.σ)/p.σ))^(1/(1-p.σ))
	ϕ  = gamma2 * (wu + r - pr * p.cbar + p.sbar ) * Lu / (ρ)
	Srh= gamma2 * (wr + r - pr * p.cbar + p.sbar ) * Lr / (ρ)

	(ρr = ρ, ϕ = ϕ/10, r = r, Lr = Lr, pr = pr, Sr = Sr)

end

function nocommute(p::LandUse.Param)

    model = Model(Ipopt.Optimizer)

    # user-defined functions
    # Sr(Lr,ρ,wr,α) = (((1 - α)/ α) * wr / ρ)^σ * Lr)
    # register(model, :Sr, 4, Sr, autodiff=true)

    @variable(model, 0.001 <= ρ ,  start = 1.0)
    @variable(model, 0.001 <= Lr <= p.L, start = 0.5)

    # nonlinear expressions
    @NLexpression(model, wu, p.θu * (p.L - Lr)^p.η)
    @NLexpression(model, Sr, (((1 - p.α)/ p.α) * wu / ρ)^p.σ * Lr)
    @NLexpression(model, pr, wu / (p.α * p.θr) * (p.α + (1-p.α) * (Lr / Sr)^((1-p.σ)/p.σ))^(1/(1-p.σ)))

    @NLexpression(model, ϕ,   p.γ / (1 + p.ϵr) * (wu + (ρ * (p.S - p.λ) / p.L) - pr * p.cbar + p.sbar) * (p.L - Lr) / ρ)
    @NLexpression(model, Srh, p.γ / (1 + p.ϵr) * (wu + (ρ * (p.S - p.λ) / p.L) - pr * p.cbar + p.sbar) * Lr / ρ)


    # objective
    @objective(model, Max, 1.0)

    # constraints

    @NLconstraint(model, (1 - p.γ) * (1 - p.ν) * (wu + (ρ * (p.S - p.λ) / p.L) - pr * p.cbar + p.sbar ) + p.ϵr * ρ * (Srh + ϕ) - p.sbar * p.L - (p.L - Lr) * p.θu == 0.0)
    @NLconstraint(model, Sr + Srh + ϕ - (p.S - p.λ) == 0.0)

    optimize!(model)

    if termination_status(model) != MOI.LOCALLY_SOLVED
        error("non-optimal exit")
    else
        return (ρ = value(ρ), Lr = value(Lr))
    end
end
