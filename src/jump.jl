


"""
solve model at current paramter p and starting at point x0
"""
function jm(p::LandUse.Param,mo::LandUse.Region,x0::NamedTuple; estimateθ = true)

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
	@NLexpression(m, wr , p.Ψ * (wu0 - p.a * (wu0^(p.tauw)) * (ϕ^(p.taul))) )
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
	@NLexpression(m, τ[i = 1:p.int_nodes], p.a * wu0^(p.tauw) * nodes[i]^(p.taul) )
	@NLexpression(m, w[i = 1:p.int_nodes], wu0 - τ[i] )
	# @warn "hard coding abs() for q function" maxlog=1
	# @NLexpression(m, q[i = 1:p.int_nodes], qr * (abs((w[i] + r_pr_csbar) / xsr))^(1.0/p.γ))
	@NLexpression(m, q[i = 1:p.int_nodes], qr * ((w[i] + r_pr_csbar) / xsr)^(1.0/p.γ))
	@NLexpression(m, H[i = 1:p.int_nodes], q[i]^p.ϵr)
	@NLexpression(m, h[i = 1:p.int_nodes], p.γ * (w[i] + r_pr_csbar) / q[i])
	@NLexpression(m, ρ[i = 1:p.int_nodes], (q[i]^(1.0 + p.ϵr)) / (1.0 + p.ϵr) )
	@NLexpression(m, cu[i = 1:p.int_nodes], (1.0 - p.γ)*(1.0 - p.ν)*(w[i] + r_pr_csbar) - p.sbar)
	@NLexpression(m, D[i = 1:p.int_nodes] , H[i] / h[i])
	@NLexpression(m, cu_input[i = 1:p.int_nodes], q[i] * H[i] * p.ϵr / (1.0+p.ϵr) )



	# integrals
	@NLexpression(m, iDensity,  (ϕ/2) * sum(mo.iweights[i] * 2π * nodes[i] * D[i] for i in 1:p.int_nodes))
	@NLexpression(m, iρ,        (ϕ/2) * sum(mo.iweights[i] * 2π * nodes[i] * ρ[i] for i in 1:p.int_nodes))
	@NLexpression(m, icu,       (ϕ/2) * sum(mo.iweights[i] * 2π * nodes[i] * cu[i] * D[i] for i in 1:p.int_nodes))
	@NLexpression(m, icu_input, (ϕ/2) * sum(mo.iweights[i] * 2π * nodes[i] * cu_input[i] for i in 1:p.int_nodes))
	@NLexpression(m, iτ,        (ϕ/2) * sum(mo.iweights[i] * 2π * nodes[i] * (wu0 - w[i]) * D[i] for i in 1:p.int_nodes))

	# objective
	# if estimateθ
		if p.it == 1
			@objective(m, Min, (Lr / p.Lt[p.it] - p.moments[p.it,:Employment_rural])^2)
		else
			@objective(m, Min, (Lr / p.Lt[p.it] - p.moments[p.it,:Employment_rural])^2 + (pr - p.moments[p.it,:P_rural])^2)
		end
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
		println("error in period $(p.it)")
		println("Termination status: $(termination_status(m))")
		println("rhor = $(value(ρr))")
		println("ϕ    = $(value(ϕ ))")
		println("r    = $(value(r ))")
		println("Lr   = $(value(Lr))")
		println("pr   = $(value(pr))")
		println("Sr   = $(value(Sr))")
		println("θr   = $(value(θr))")
		println("θu   = $(value(θu))")
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
