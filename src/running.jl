

xtrace(x::NLsolve.SolverResults) = vcat([(x.trace[i].metadata["x"])' for i in 1:x.iterations]...)
ftrace(x::NLsolve.SolverResults) = vcat([(x.trace[i].metadata["f(x)"])' for i in 1:x.iterations]...)


function solve_once(p::Param,m0::Model,x0::Vector{Float64})
	nlsolve((F,x) -> solve!(F,x,p,m0),x0,iterations = p.iters,store_trace = p.trace, extended_trace = p.trace)
end


"""
run Single region model for all time periods
"""
function run(p::Param; estimateθ = false)

	setperiod!(p,1)
	# x0 = startval(p)
	x0 = nearstart(p)

	# m = Region(p)
	# x = jm(p,m,x0)
	# px = x0.pr / p.moments[1,:P_rural]
	# transform!(p.moments,:P_rural => (x -> x .* px) => :P_rural)

	sols = NamedTuple[]
	push!(sols,x0)
	M = Region[]

	for it in 1:length(p.T)
		# println(it)
		setperiod!(p,it)
		m = Region(p)
		x = jm(p,m,sols[it], estimateθ = estimateθ)
		push!(sols,x)
		if it == 1
			p.ϕ1 = x.ϕ * p.ϕ1x
		end
		update!(m,p,[x...])

		push!(M,m)
		if it == 1
			# adjust relative price in data to first period solution
			px = x.pr / p.moments[1,:P_rural]
			# transform!(p.moments,:P_rural => (x -> x .* px) => :P_rural)
		end

	end

	(sols[2:end],M,p) # solutions, models, and parameter

end


"""
run Single region model for all time periods
"""
function run(x0::NamedTuple, p::Param)

	sols = NamedTuple[]
	push!(sols,x0)
	M = Region[]

	for it in 1:length(p.T)
		# println(it)
		setperiod!(p,it)
		m = Region(p)
		x = jm(p,m,sols[it])
		push!(sols,x)
		if it == 1
			p.ϕ1 = x.ϕ * p.ϕ1x
		end
		update!(m,p,[x...])

		push!(M,m)
		if it == 1
			# adjust relative price in data to first period solution
			# px = x.pr / p.moments[1,:P_rural]
			# transform!(p.moments,:P_rural => (x -> x .* px) => :P_rural)
		end

	end

	(sols[2:end],M,p) # solutions, models, and parameter

end

"""
run Multi-region model for all time periods starting from 
the single city starting value. Works only for not too different θu values.
"""
function runk(;par = Dict(:K => 2,:kshare => [0.5,0.5], :factors => [1.0,1.0], :gs => zeros(2)), estimateθ = false)

	# get single city solution in first period
	p = LandUse.Param(par = par, use_estimatedθ = false)
	@assert p.K > 1

	setperiod!(p,1)
	x0 = nearstart(p)
	m = Region(p)
	x0 = jm(p,m,x0, estimateθ = false)
	update!(m,p,[x0...])

	x = Float64[]
	push!(x, m.Lr / m.Sr)
	push!(x, m.r)
	push!(x, m.pr)
	for ik in 1:p.K
		push!(x,m.Sr)
	end
	for ik in 1:p.K
		push!(x,m.Lu)
	end
	for ik in 1:p.K
		push!(x,p.θu)
	end
	runk_impl(x,p, estimateθ = estimateθ)
end




function check(it; d1 = 0.04, d2= 1.0)
	par = Dict(:d1 => d1, :d2 => d2,:K => 2,:kshare => [0.5,0.5], :factors => [1.0,1.0], :gs => zeros(2))

	x0,M,p = run(Param(par = par))
	m = M[it]  # period 
	setperiod!(p,it)
	sols =Vector{Float64}[]
	C = Country[]  # an emtpy array of countries


	x = Float64[]
	push!(x, m.Lr / m.Sr)
	push!(x, m.r)
	push!(x, m.pr)
	for ik in 1:p.K
		push!(x,m.Sr)
	end
	for ik in 1:p.K
		push!(x,m.Lu)
	end
	for ik in 1:p.K
		push!(x,p.θu)
	end
	push!(sols, x)
	# for it in 1:length(p.T)
		# println(it)
		# setperiod!(p,it)
		c = Country(p)  # performs scaling of productivity upon creation
		xmod = jc(c,sols[1],estimateθ = false) # returns a JuMP model as last element
		if termination_status(xmod[end]) == MOI.LOCALLY_SOLVED
			# clean up results and save
			x,ϕs = xmod[1], xmod[2]
		else		
			println("period = $it")
			println(termination_status(xmod[end]))  # error
			println(JuMP.all_variables(xmod[end]))
			return JuMP.primal_feasibility_report(xmod[end])
		end
		update!(c,x,estimateθ = false)

		push!(C,c)
		# push!(ϕvs,ϕs)
		# push!(dϕvs,dϕs)
	
	# (sols,C,p) # solutions, models, and parameter

	m.ϕ - c.R[1].ϕ
end


function check2()
	par = Dict(:d1 => 0.04, :d2 => 1.0,:K => 2,:kshare => [0.5,0.5], :factors => [1.0,1.0], :gs => zeros(2))

	x0,M,p = run(Param(par = par))
	m = M[1]  # period 2
	setperiod!(p,1)
	sols =Vector{Float64}[]
	C = Country[]  # an emtpy array of countries


	x = Float64[]
	push!(x, m.Lr / m.Sr)
	push!(x, m.r)
	push!(x, m.pr)
	for ik in 1:p.K
		push!(x,m.Sr)
	end
	for ik in 1:p.K
		push!(x,m.Lu)
	end
	for ik in 1:p.K
		push!(x,p.θu)
	end
	push!(sols, x)
	for it in 1:length(p.T)
		# println(it)
		setperiod!(p,it)
		c = Country(p)  # performs scaling of productivity upon creation
		xmod = jc(c,sols[it],estimateθ = false) # returns a JuMP model as last element
		if termination_status(xmod[end]) == MOI.LOCALLY_SOLVED
			# clean up results and save
			x,ϕs = xmod[1], xmod[2]
		else		
			println("period = $it")
			println(termination_status(xmod[end]))  # error
			println(JuMP.all_variables(xmod[end]))
			return JuMP.primal_feasibility_report(xmod[end])
		end

		# x,ϕs = jc(c,sols[it],estimateθ = estimateθ)
		# x,ϕs,dϕs = jc(c,sols[it],estimateθ = estimateθ)
		push!(sols,x)
		if it == 1
			for ik in 1:p.K
				c.pp[ik].ϕ1 = ϕs[ik] * c.pp[ik].ϕ1x
			end
		else
			for ik in 1:p.K
				c.pp[ik].ϕ1 = C[1].R[ik].ϕ * c.pp[ik].ϕ1x
			end
		end
		# overwrite θu if estimated
		# if estimateθ
		# 	for ik in 1:p.K
		# 		c.pp[ik].θu = x[3 + 2p.K + ik]
		# 	end
		# end

		update!(c,x,estimateθ = false)

		push!(C,c)
		# push!(ϕvs,ϕs)
		# push!(dϕvs,dϕs)
	end
	
	(sols,C,p) # solutions, models, and parameter
end


function runk_impl(x0::Vector,p::Param; estimateθ = false)
	sols =Vector{Float64}[]
	# ϕvs =Vector{Float64}[]
	# dϕvs =Vector{Float64}[]
	push!(sols,x0)

	C = Country[]  # an emtpy array of countries

	for it in 1:length(p.T)
		# println(it)
		setperiod!(p,it)
		c = Country(p)  # performs scaling of productivity upon creation
		xmod = jc(c,sols[it],estimateθ = estimateθ) # returns a JuMP model as last element
		if termination_status(xmod[end]) == MOI.LOCALLY_SOLVED
			# clean up results and save
			x,ϕs = xmod[1], xmod[2]
		else		
			println("period = $it")
			println(termination_status(xmod[end]))  # error
			println(JuMP.all_variables(xmod[end]))
			return JuMP.primal_feasibility_report(xmod[end])
		end

		# x,ϕs = jc(c,sols[it],estimateθ = estimateθ)
		# x,ϕs,dϕs = jc(c,sols[it],estimateθ = estimateθ)
		push!(sols,x)
		if it == 1
			for ik in 1:p.K
				c.pp[ik].ϕ1 = ϕs[ik] * c.pp[ik].ϕ1x
			end
		else
			for ik in 1:p.K
				c.pp[ik].ϕ1 = C[1].R[ik].ϕ * c.pp[ik].ϕ1x
			end
		end
		# overwrite θu if estimated
		if estimateθ
			for ik in 1:p.K
				c.pp[ik].θu = x[3 + 2p.K + ik]
			end
		end

		update!(c,x,estimateθ = estimateθ)

		push!(C,c)
		# push!(ϕvs,ϕs)
		# push!(dϕvs,dϕs)
	end
	# (sols,C,p, ϕvs, dϕvs) # solutions, models, and parameter
	(sols,C,p) # solutions, models, and parameter
end


"helper function to prepare country param"
function startvals_par(K::Int; θus = 1.2:0.01:1.35)
	if K == 2
		θus = 1.2:0.01:1.32
		facs = [1.0, θus[1]]
	elseif K == 3	
		facs = [1.0, 1.07, θus[1]]
	elseif K == 4	
		facs = [1.0, 1.07, 1.1, θus[1]]
	elseif K == 5	
		facs = [1.0, 1.07, 1.08, 1.09, θus[1]]
	end
		
	(par = Dict(:K => K,:kshare => [1/K for i in 1:K], :factors => facs, :gs => zeros(K)), θus = θus)
end

"""
collect valid starting values for the k country case by reusing the first period solution
	of each preceding step in θu along the sequence
"""
function startvals_impl(K::Int; θu = 1.2:0.01:1.35)
	
	par, θus = startvals_par(K,θus = θu)

	x,C,p = runk(par = par)

	x0 = Vector{Float64}[] 
	push!(x0, x[2])
	for (i,θu) in enumerate(θus)
		# println("θu = $θu")
		par[:factors][K] = θu
		p = LandUse.Param(par = par)
		x1,C1,p1 = runk_impl(x0[i],p)
		push!(x0,x1[2])
	end
	(par = par, x0 = x0[end])
end


"""
	startvals_k()

Find feasible starting values for multi city case (2 to 5 cities) and write to disk. By default return the saved dict with start param and initial guess vector `x0`.
"""
function startvals_k(K::Int; overwrite = false)
	if overwrite
		# recompute all starting values and write to disk
		d = Dict()
		for k in 2:K
			@info "doing case k=$k"
			d[k] = startvals_impl(k)
		end
		bson(joinpath(@__DIR__,"..","out","multistarts.bson"), d)

	else
		# read from disk
		d = BSON.load(joinpath(@__DIR__,"..","out","multistarts.bson"))
	end
	d[K]
end

"5 country case"
function k5()
	(par, x0) = startvals_k(5)
	par[:gs] = zeros(5)
	p = Param(par = par)
	runk_impl(x0,p)
end


function k3()

	(par, x0) = startvals_impl(3, θu = 1.2:0.01:1.24)
	par[:gs] = [0.0,0.001,0.01]
	
	p = Param(par = par)
	x,C,p = runk_impl(x0,p)
	d = dataframe(C)
	gg = groupby(select(filter(x -> x.region .∈ Ref([1,3]), d), :region, :year, :Lu), :year)
    combine(gg, :Lu => (x -> maximum(x) / minimum(x)) => :rel_Lu)
end

"""
return model relative population in each year to largest city
largest city is city 1
"""
function relpop(C::Vector{Country})
	d = dataframe(C)
	cla = select(C[1].pp[1].citylist, :rank, :group, :ngroup)
	transform!(cla, [:group, :ngroup] => ((a,b) -> string.(a,"(n=",b,")")) => :grouplabel)

	d2 = leftjoin(select(d, :year, :region, :Lu), C[1].pp[1].citylist, on = [:region => :rank, :year])
	gg = groupby(d2, :year)
    g2 = combine(gg, [:Lu, :region] => ((a,b) -> a ./ a[b .== 1]) => :rel_Lu, :region, :Lu, :grouplabel)
	combine(groupby(g2, [:grouplabel, :year]),  :rel_Lu => mean, :region) # mean amongst groups
end

function k(K;pars = Dict(),estimateθ = true)
	LandUse.runk(par = merge(Dict(:K => K,:kshare => [1/K for i in 1:K], :factors => ones(K), :gs => zeros(K)), pars),estimateθ = estimateθ)
end

function k20(;overwrite = false)
	if overwrite
		x,C,p = LandUse.runk(par = Dict(:K => 20,:kshare => [1/20 for i in 1:20], :factors => ones(20), :gs => zeros(20)),estimateθ = true)
		d = dataframe(C)

		FileIO.save(joinpath(LandUse.dboutdata, "k20.jld2"), Dict("df" => d))
		
		x,C,p,d
	else
		df = FileIO.load(joinpath(LandUse.dboutdata, "k20.jld2"))
		df["df"]
	end

	# K = 20
	# # par = Dict(:K => K, :kshare => [1/K for i in 1:K], :factors => [1.09, 1.02, 1.02, 1.01,1.01, [1.005 for i in 1:7]...,ones(8)...], :gs => [0.003,zeros(K-1)...])
	# par = Dict(:K => K, :kshare => [1/K for i in 1:K], :factors => ones(K), :gs => zeros(K))
	# # par = Dict(:K => K, :kshare => [1/K for i in 1:K], :factors => [1.1, 1.01, 1.01, 1.005,1.005, [1.0025 for i in 1:7]...,ones(8)...], :gs => [0.0025,zeros(K-1)...])
	# # par = Dict(:K => K, :kshare => [1/K for i in 1:K], :factors => [1.05, 1.01, 1.01, 1.005,1.005, [1.0025 for i in 1:7]...,ones(8)...], :gs => [0.002,zeros(K-1)...])
	# # par = Dict(:K => K, :kshare => [1/K for i in 1:K], :factors => [1.02, 1.001, 1.001, 1.01,1.01, [1.005 for i in 1:7]...,ones(8)...], :gs => [0.001,zeros(K-1)...])
	# # x,C,p = runk(par = par)
	# runk(par = par)
	# # relpop(C)
	# # dd = relpop(C)
	# @df subset(dd, :region => x-> x.> 1, :year => x-> x.> 1870) plot(:year, :rel_Lu_mean, group = :grouplabel, title = "model relative to paris")
	
end

function runm()
	run(Param())
end
function plot1()
	x,M,p = run(Param())
	ts_plots(M,p)
end
function plot1cs(it)
	x,M,p = run(Param())
	cs_plots(M[it],p,it)
end
function dash(it)
	x,M,p = run(Param())
	dashboard(M,p,it)
end
function cdash(it)
	x,M,p = k3()
	dashboard(M,it)
end
function export_params()
	x,M,p = runm()
	latex_param()
	d = DataFrame(year = collect(p.T), thetau = p.θut, thetar = p.θrt, pr = [M[it].pr for it in 1:length(M)], Lt = p.Lt)
	CSV.write(joinpath(dbtables,"export_theta_pr.csv"),d)

	x0 = LandUse.nearstart(Param())
	CSV.write(joinpath(dbtables,"export_x0.csv"),DataFrame([x0]))
end
