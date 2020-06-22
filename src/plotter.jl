
function setup_color(l::Int)
    # setup the colors
    grad   = ColorSchemes.inferno
    colors = [grad[k] for k in range(0,stop=1,length=l)]
    # colors = repeat(colors, 1, length(colors))
    return colors
end

function ts_plots(M,p::Param;fixy = false)
	d = dataframe(M,p.T)
	dd = Dict()
	df = @linq d |>
		 select(:year,:Ch,:Cu,:Cr,:C ) |>
		 transform(h = :Ch ./ :C,u = :Cu ./ :C,r = :Cr ./ :C) |>
		 select(:year, :h, :u , :r)
	ds = stack(df, Not(:year))

	# marker
	mmark = (:circle, 4)

	dd[:spending] = @df ds plot(:year,:value, group = :variable,
			   linewidth = 2, title = "Spending Shares",
			   ylims = (0.0,0.85), marker = mmark, legend = :right)

	ds2 = stack(select(d,:year,:Lu, :Lr), Not(:year))
	dd[:pop] = @df ds2 plot(:year, :value, group = :variable,
					 linewidth = 2, title = "Population",
					 ylims = (0,1), marker = mmark, legend = :topleft)

	dd[:ρ0] = @df d plot(:year, :ρ0, linewidth = 2, color = "black", leg =false, title = "Central Land Values", marker = mmark, ylims = fixy ? (0.0,29.0) : false )
 	dd[:ρr] = @df d plot(:year, :ρr, linewidth = 2, color = "black", leg =false, title = "Fringe Land Values", marker = mmark, ylims = fixy ? (0.19,0.36) : false )
 	dd[:q0] = @df d plot(:year, :q0, linewidth = 2, color = "black", leg =false, title = "Central House Prices", marker = mmark, ylims = fixy ? (0.0,5.0) : false )
 	dd[:qr] = @df d plot(:year, :qr, linewidth = 2, color = "black", leg =false, title = "Fringe House Prices", marker = mmark, ylims = fixy ? (0.5,1.5) : false )
	dd[:hr] = @df d plot(:year, :hr, linewidth = 2, color = "black", leg =false, title = "Fringe Housing Demand", marker = mmark, ylims = fixy ? (0.0,5.0) : false )

	df = @linq d |>
		 select(:year,:r, :y ,:iq) |>
		 transform(r_y = 100*(:r ./ :y), ru_y = 100*(:iq ./ :y))
	dd[:r_y] = @df df plot(:year, [:r_y, :ru_y], labels = ["Total" "Urban Only"],
	            linewidth = 2, title = "Land Rents over Income (%)", marker = mmark, ylims = fixy ? (0,50) : false)


	dtemp = transform(select(d, :year, :r, :ρ0, :ρr, :q0, :qr, :hr, :Hr),
	                  :r => ((x) -> 100*(x ./ x[7])) => :r,
	                  :ρ0 => ((x) -> 100*(x ./ x[7])) => :ρ0,
	                  :ρr => ((x) -> 100*(x ./ x[7])) => :ρr,
	                  :q0 => ((x) -> 100*(x ./ x[7])) => :q0,
	                  :hr => ((x) -> 100*(x ./ x[7])) => :hr,
	                  :Hr => ((x) -> 100*(x ./ x[7])) => :Hr,
	                  :qr => ((x) -> 100*(x ./ x[7])) => :qr
					  )
	dd[:r] = @df dtemp plot(:year, :r, linewidth = 2, color = "black", leg =false, title = "Land Rents (1945=100)", marker = mmark, ylims = fixy ? (95,500) : false)
	dd[:r_rho] = @df dtemp plot(:year, [:r, :ρ0],linewidth = 2, lab = [L"r" L"\rho_0" ],
	                           title = "Land Rents and Central Land Values (1945=100)",
							   leg = :topleft,
							   marker = mmark, ylims = fixy ? (50,500) : false)

    dd[:q0_qr] = @df dtemp plot(:year, [:q0, :qr],linewidth = 2, lab = [L"q_0" L"q_r"],
	                           title = "Central and Fringe House Prices (1945=100)",
							   leg = :topleft,
							   marker = mmark, ylims = fixy ? (85,150) : false)
	# is missing density!
    # dd[:hr100] = @df dtemp plot(:year, :hr, linewidth = 2, color = "black", leg =false, title = "Fringe Housing Demand (1945=100)", marker = mmark)
    # dd[:Hr100] = @df dtemp plot(:year, :Hr, linewidth = 2, color = "black", leg =false, title = "Fringe Housing Supply (1945=100)", marker = mmark)



	ds2 = stack(select(d,:year,:qr, :ρr), Not(:year))

	d3  = @linq d |>
			select(:year,:θu, :θr) |>
			transform(theta_u = :θu, theta_r = :θr) |>
			select(:year,:theta_u, :theta_r)

	# ds3 = stack(d3, Not(:year))
	# pl3 = @df ds3 plot(:year, :value, group = :variable,
	# 			      linewidth = 2, title = "Consumption",
	# 				  ylims = (0.0,3.5))

	ds3 = stack(d3, Not(:year))
	dd[:productivity] = @df ds3 plot(:year, :value, group = :variable,
					  linewidth = 2, title = "Productivity", ylims = fixy ? (0,20) : false)
	ds4 = stack(select(d,:year,:ϕ), Not(:year))
	dd[:phi] = @df ds4 plot(:year, :value, group = :variable,
					 linewidth = 2, title = "City Size", color = "black",
					 leg = false, marker = mmark, ylims = fixy ? (0.0,0.15) : false)
    dd[:Sr] = @df d plot(:year, :Sr,
				  linewidth = 2, color = "black",title = "Agricultural Land",
				  leg = false, marker = mmark, ylims = fixy ? (0.6,1.0) : false)
	# ds4 = stack(select(d,:year,:pr), Not(:year))
	df4 = @linq d |>
		# transform(avgd = :Lu ./ :ϕ, tauphi = (1 .- map(x -> τ(x,x,p),:ϕ) .* :ϕ) ./ :pr)
		transform(avgd = :Lu ./ :ϕ, tauphi = (1 .- map(x -> τ(x,x,p),:ϕ) .* :ϕ) )
		# transform(avgd = :Lu ./ :ϕ, tauphi = (1 .- p.τ .* :ϕ))

	ds4 = stack(select(df4,:year,:d0,:dq1, :dq2, :dq3, :dq4, :dr, :avgd), Not(:year))
	ds5 = stack(select(df4,:year,:avgd), Not(:year))
	ds6 = stack(select(df4,:year,:tauphi), Not(:year))
	# ds4 = stack(select(df4,:year, :avgd), Not(:year))
	dd[:densities] = @df ds4 plot(:year, :value, group = :variable,
					 linewidth = 2, title = "Densities", leg = :topright, ylims = fixy ? (0,130) : false)
	dd[:avdensity] = @df ds5 plot(:year, :value, group = :variable,
					 linewidth = 2, title = "Avg Density", marker = mmark, legend = false,color = "black", ylims = fixy ? (0,50) : false)
	dd[:tauphi] = @df ds6 plot(:year, :value, group = :variable,
					  linewidth = 2, title = "1 - tau(phi)",leg = false)
	# plot(pl,pl2,pl3,pl4,pl5, layout = (1,5),size = (700,250))
	# plot(pl2,pl6,pl4,pl7 ,layout = (1,4),size = (900,250))
	dd
end


"make separate time series plot for each region"
function plot_ts(C::Vector{Country})
	df = dataframe(C)
	K = length(C[1].R)
	# vars = (:ρr, :qr, :Lr, :Lu, :wu0, :wr, :Sr, :Srh, :r, :pr, :ϕ, :icu_input, :iDensity, :icu, :icr, :iτ, :iq, :iy)
	# @df df plot(:year, cols(2:size(df,2)))
	s = stack(df, Not([:year, :region]))
	# sims = [(:Lr,:Lu); (:Sr, :ϕ, :Srh); (:qr, :r); (:wr , :wu0)]
	sims = [[:Lr,:Lu], [:Sr, :ϕ, :Srh], [:qr, :r], [:wr , :wu0]]
	titles = ["Labor"; "Land"; "Rents"; "Wages"]
	nms = [[L"L_r" L"L_u"], [L"S_r" L"\phi" L"S_{rh}"] , [L"q" L"r"], [L"w_r" L"w_u"]]

	plts = Any[]
	for ir in 1:K
		plt = Any[]
		for i in 1:length(sims)
			x = @linq s |>
				where((:variable .∈ Ref(sims[i])) .& (:region .== ir))
			px = @df x plot(:year, :value, group = :variable,
							title = titles[i],
							titlefontsize=10,
							label=nms[i],
							legend = :topleft,
							linewidth=2,marker = (:circle,3))
			push!(plt, px)
		end
		push!(plts, plot(plt...,layout = (2,2)))
	end
	plts
end

function impl_plot_slopes(C::Vector{Country})
	d = dataframe(C)
	d.lϕ = log.(d.ϕ)
	d.lu = log.(d.Lu)
	gd = groupby(d,:year)
	gd = combine(gd, AsTable([:lϕ, :lu]) => (x -> round.(diff(vcat(extrema(x.lϕ)...)) ./ diff(vcat(extrema(x.lu)...)),digits = 1)) => :slope)
	d  = innerjoin(d,gd,on = :year)
	transform!(d, AsTable([:year, :slope]) => (x -> string.(x.year) .* ": slope=" .* string.(x.slope) ) => :year_s)

	cols = range(colorant"red",colorant"blue",length = length(unique(d.year)))
	dd = select(d, :year_s, :region, :lϕ, :lu)
	pl2 = @df dd plot(:lu,:lϕ,group = :year_s,
					ylab = L"\log \phi",
					xlab = L"\log L_u",
					marker = (:circle, 4),
					colour = cols',
					legend = :outertopright,
					# ylims = (-3.8,-2.7),
					# xlims = K == 3 ? (-1.5,-0.8) : (-1.0,-0.5),
					)
	pl2
end

"time series showing all regions together"
function impl_plot_ts_all(C::Vector{Country})
	df = dataframe(C)
	K = length(C[1].R)
	# vars = (:ρr, :qr, :Lr, :Lu, :wu0, :wr, :Sr, :Srh, :r, :pr, :ϕ, :icu_input, :iDensity, :icu, :icr, :iτ, :iq, :iy)
	# @df df plot(:year, cols(2:size(df,2)))
	s = stack(df, Not([:year, :region]))
	# sims = [(:Lr,:Lu); (:Sr, :ϕ, :Srh); (:qr, :r); (:wr , :wu0)]
	# sims = [[:Lu], [:ϕ], [:pop], [:ρ0], [:q0], [:Srh]]
	# sims = ["Lu", "ϕ", "pop", "ρ0", "q0","Srh", "Sr", "Hr"]
	sims = ["Lu", "ϕ", "pop","dbar"]
	# titles = ["Urban Labor"; "Urban Size"; "Total pop"; "Central Land values"; "Central House prices"; "Rural Housing (Srh)"; "Rural Land (Sr)"; "r Housing supply"]
	titles = ["Urban Labor"; "Urban Size"; "Total pop"; "Avg Density"]

	plts = Any[]
	for i in 1:length(sims)
			x = @linq s |>
				where((:variable .== sims[i]))
			px = @df x plot(:year, :value, group = :region,
							title = titles[i],
							titlefontsize=10,
							legend = i == 1 ? true : false,
							linewidth=2,marker = (:circle,3),
							legendfontsize = 5)
			push!(plts, px)

	end
	# push!(plts,plot())  # fill up with empty
	# push!(plts, scatter((1:3)', xlim = (4,5), legend = true, framestyle = :none, labels = ["1" "2" "3"]))
	# plot(plts...,layout = (1,4))
	plts


end


"single spatial cross sections region"
function plot_space(m::Region,p::Param,it::Int; fixy = false)

	lvec = collect(range(0.,area(m),length=100))  # space for region ik
	setperiod!(p,it)  # set time on parameter!
	plts = Any[]
	if fixy
		push!(plts, plot(lvec, [ϵ(i,m.ϕ,p) for i in lvec],title = L"\epsilon(l)", ylims = (2,4.1), linewidth = 2, leg = false))
		push!(plts, plot(lvec, [D(i,p,m) for i in lvec],title = L"D(l)", ylims = (-3,60), linewidth = 2, leg = false))
		push!(plts, plot(lvec, [H(i,p,m) for i in lvec],title = L"H(l)", ylims = (-0.1,15), linewidth = 2, leg = false))
		push!(plts, plot(lvec, [ρ(i,p,m) for i in lvec],title = L"\rho(l)", ylims = (-0.1,10), linewidth = 2, leg = false))
		push!(plts, plot(lvec, [w(m.Lu,i,m.ϕ,p) for i in lvec],title = L"w(l)", ylims = (-0.1,5.5), linewidth = 2, leg = false))
		# push!(plts, plot(lvec, [h(i,p,m) for i in lvec],title = L"h(l)", ylims = (-0.1,1.5), linewidth = 2, leg = false))
		push!(plts, plot(lvec, [q(i,p,m) for i in lvec],title = L"q(l)", ylims = (0.5,5.5), linewidth = 2, leg = false))
	else
		push!(plts, plot(lvec, [ϵ(i,m.ϕ,p) for i in lvec],title = L"\epsilon(l)",linewidth = 2, leg = false))
		push!(plts, plot(lvec, [D(i,p,m) for i in lvec],title = L"D(l)",linewidth = 2, leg = false))
		push!(plts, plot(lvec, [H(i,p,m) for i in lvec],title = L"H(l)",linewidth = 2, leg = false))
		push!(plts, plot(lvec, [ρ(i,p,m) for i in lvec],title = L"\rho(l)",linewidth = 2, leg = false))
		push!(plts, plot(lvec, [w(m.Lu,i,m.ϕ,p) for i in lvec],title = L"w(l)",linewidth = 2, leg = false))
		# push!(plts, plot(lvec, [h(i,p,m) for i in lvec],title = L"h(l)",linewidth = 2, leg = false))
		push!(plts, plot(lvec, [q(i,p,m) for i in lvec],title = L"q(l)", linewidth = 2, leg = false))

	end

	ti = plot(title = "Spatial Cross-Section in $(p.T[it])", grid = false, showaxis = false, bottom_margin = -30Plots.px)
	pl = plot(ti,plot(plts...,layout = (2,3), link = :x), layout = @layout([A{0.05h}; B]))
	pl
end

"plot repeated spatial cross sections for each region"
function anim_space(C::Vector{Country},pp::Vector{Param})
	K = length(C[1].R)

	anims = Plots.Animation[]
	for ik in 1:K
		plt = Any[]
		p = pp[ik]
		# pl0 = plot(leg = false)
		anim = Animation()

		for (jt,it) in enumerate(C[1].T)
			reg = C[jt].R[ik]
			pl0 = plot_space(reg,p,jt,fixy = true)
			frame(anim,pl0)
		end
		gif(anim, joinpath(dbpath,"x-section-$ik.gif"),fps=0.5)
		push!(anims,anim)
	end
	anims
end


"single spatial cross sections region"
function plot_rho_space(m::Region,p::Param,it::Int; fixy = false)

	lvec = collect(range(0.,area(m),length=100))  # space for region ik
	setperiod!(p,it)  # set time on parameter!
	plts = Any[]
	if fixy
		push!(plts, plot(lvec, [ϵ(i,m.ϕ,p) for i in lvec],title = L"\epsilon(l)", ylims = (2,4.1), linewidth = 2, leg = false))
		# push!(plts, plot(lvec, [D(i,p,m) for i in lvec],title = L"D(l)", ylims = (-3,60), linewidth = 2, leg = false))
		# push!(plts, plot(lvec, [H(i,p,m) for i in lvec],title = L"H(l)", ylims = (-0.1,15), linewidth = 2, leg = false))
		push!(plts, plot(lvec, [ρ(i,p,m) for i in lvec],title = L"\rho(l)", ylims = (-0.1,10), linewidth = 2, leg = false))
		# push!(plts, plot(lvec, [w(m.Lu,i,m.ϕ,p) for i in lvec],title = L"w(l)", ylims = (-0.1,5.5), linewidth = 2, leg = false))
		# push!(plts, plot(lvec, [h(i,p,m) for i in lvec],title = L"h(l)", ylims = (-0.1,1.5), linewidth = 2, leg = false))
		# push!(plts, plot(lvec, [q(i,p,m) for i in lvec],title = L"q(l)", ylims = (0.5,5.5), linewidth = 2, leg = false))
	else
		push!(plts, plot(lvec, [ϵ(i,m.ϕ,p) for i in lvec],title = L"\epsilon(l)",linewidth = 2, leg = false))
		# push!(plts, plot(lvec, [D(i,p,m) for i in lvec],title = L"D(l)",linewidth = 2, leg = false))
		# push!(plts, plot(lvec, [H(i,p,m) for i in lvec],title = L"H(l)",linewidth = 2, leg = false))
		push!(plts, plot(lvec, [ρ(i,p,m) for i in lvec],title = L"\rho(l)",linewidth = 2, leg = false))
		# push!(plts, plot(lvec, [w(m.Lu,i,m.ϕ,p) for i in lvec],title = L"w(l)",linewidth = 2, leg = false))
		# push!(plts, plot(lvec, [h(i,p,m) for i in lvec],title = L"h(l)",linewidth = 2, leg = false))
		# push!(plts, plot(lvec, [q(i,p,m) for i in lvec],title = L"q(l)", linewidth = 2, leg = false))

	end

	ti = plot(title = "X-Section", grid = false, showaxis = false, bottom_margin = -30Plots.px)
	pl = plot(ti,plot(plts...,layout = (2,1), link = :x), layout = @layout([A{0.05h}; B]))
	pl
end

function TS_impl(s::DataFrame; year = nothing, xlim = nothing,ylim = nothing,tstring = nothing)

	# sims = [[:Lr,:Lu], [:Sr, :ϕ, :Srh], [:qr, :r], [:wr , :wu0]]
	sims = [["Lr","Lu"], ["ϕ"], ["r"], ["wr" , "wu0"]]
	# sims = [[:Lr,:Lu], [:ϕ], [:r], [:q0, :ρ0] , [:qr, :ρr], [:wr , :wu0],[:Hr , :hr]]
	# titles = ["Labor"; "Land"; "Rents"; "Central prices"; "Fringe prices"; "wages"; "Housing at phi"]
	nms = [[L"L_r" L"L_u"], [L"S_r" L"\phi" L"S_{rh}"] , [L"q" L"r"], [L"w_r" L"w_u"]]
	# nms = [[L"L_r" L"L_u"], L"\phi" , L"r", [L"q_0" L"\rho_0"], [L"q_r" L"\rho_r"], [L"w_r" L"w_u = \theta_u"], ["Supply" "Demand"]]
	nms = [[L"L_r" L"L_u"], L"\phi" , L"r", [L"q_0" L"q_r"]]

	plt = Any[]
	if isnothing(year)
		for i in 1:length(sims)
			x = @linq s |>
				where((:variable .∈ Ref(sims[i])))
			px = @df x plot(:year, :value, group = :variable,
							# title = titles[i],
							titlefontsize=10,
							label=nms[i],
							legend = i < 3 ? :bottomright : :topleft,
							linewidth=2,marker = (:circle,3), legendfontsize = 8)
			push!(plt, px)
		end
		# push!(plt,plot())
		# push!(plt,plot())
		# ti = plot(title = "Time Series", grid = false, showaxis = false, bottom_margin = -30Plots.px)
		# pl = plot(ti,plot(plt...,layout = (2,3), link = :x), layout = @layout([A{0.05h}; B]))
		pl = plot(plt...,layout = (2,2), link = :x)
		# return extrema of each subplot
		xlims = [plt[i].subplots[1][:xaxis][:extrema] for i in 1:length(plt)]
		ylims = [plt[i].subplots[1][:yaxis][:extrema] for i in 1:length(plt)]
		xlims = [[xlims[i].emin,xlims[i].emax] for i in 1:length(plt)]
		ylims = [[ylims[i].emin,ylims[i].emax] for i in 1:length(plt)]
		xlims = [(xlims[i][1] .- 0.05*diff(xlims[i])[1],xlims[i][2] .+ 0.05*diff(xlims[i])[1]) for i in 1:length(plt)]
		ylims = [(ylims[i][1] .- 0.05*diff(ylims[i])[1],ylims[i][2] .+ 0.05*diff(ylims[i])[1]) for i in 1:length(plt)]
		return pl,xlims,ylims
	else
		@assert isa(xlim, AbstractArray) && isa(ylim, AbstractArray)
		for i in 1:length(sims)
			x = @linq s |>
				where((:variable .∈ Ref(sims[i])) .& (:year .<= year))
			legpos = (sims[i][1] == :Lr) ? :bottomright : :topleft
			px = @df x plot(:year, :value, group = :variable,
							title = titles[i],
							titlefontsize=10,
							label=nms[i],
							xlims = xlim[i],
							ylims = ylim[i],
							legend = legpos,
							linewidth=2,marker = (:circle,3))
			push!(plt, px)
		end
		if isnothing(tstring)
			ti = "Time Series until $year"
		else
			ti = tstring
		end
		# ti = plot(title = ti, grid = false, showaxis = false, bottom_margin = -30Plots.px)
		# pl = plot(ti,plot(plt...,layout = (2,3), link = :x), layout = @layout([A{0.05h}; B]))
		pl = plot(plt...,layout = (2,2), link = :x)

		return pl
	end

end

function doit()
	p = Param()
	x,M,p = run(p)
	plot_ts_xsect(M,p,1)
end

function plot_ts(M::Vector{Region},p::Param,it::Int)
	df = dataframe(M,p.T)
	s = stack(df, Not([:year]))
	# prepare a TS

	ts0,xlims,ylims = TS_impl(s)

	yr = p.T[it]

	LandUse.setperiod!(p,1)
	θ1 = p.θu
	LandUse.setperiod!(p,14)
	θ14 = p.θu

	# make a TS up to period jt: keeping extrema fixed, however, so the plot "grows" nicely
	pl = TS_impl(s, year = yr, xlim = xlims, ylim = ylims , tstring = "first and last period prod: [$θ1,$(round(θ14,digits = 2))]")
	plot(pl)
end

function plot_ts(M::Vector{Region},p::Param)
	df = dataframe(M,p.T)
	s = stack(df, Not([:year]))
	# prepare a TS
	LandUse.setperiod!(p,1)
	θ1 = p.θu
	LandUse.setperiod!(p,14)
	θ14 = p.θu

	# make a TS up to period jt: keeping extrema fixed, however, so the plot "grows" nicely
	ts0,xlims,ylims = TS_impl(s,tstring = "first and last period prod: [$θ1,$(round(θ14,digits = 2))]")
	ts0
end

function plot_ts_xsect(M::Vector{Region},p::Param,it::Int)
	df = dataframe(M,p.T)
	s = stack(df, Not([:year]))
	# prepare a TS

	ts0,xlims,ylims = TS_impl(s)

	yr = p.T[it]

	LandUse.setperiod!(p,1)
	θ1 = p.θu
	LandUse.setperiod!(p,14)
	θ14 = p.θu

	# make a TS up to period jt: keeping extrema fixed, however, so the plot "grows" nicely
	ts = TS_impl(s, year = yr, xlim = xlims, ylim = ylims , tstring = "first and last period prod: [$θ1,$(round(θ14,digits = 2))]")

	# make a x-section
	xs = plot_rho_space(M[it],p,it, fixy = true)

	# combine both
	# o = plot(ts, xs, layout = (1,2), size = (1300,500))
	o = plot(ts, xs, layout = @layout([A{0.75w} B]))
	o

end






function single_TS(s::DataFrame, ik::Int ; year = nothing, xlim = nothing,ylim = nothing)
	# sims = [[:Lr,:Lu], [:Sr, :ϕ, :Srh], [:qr, :r], [:wr , :wu0]]
	sims = [[:Lr,:Lu], [:Sr, :ϕ, :Srh], [:r], [:wr , :wu0]]
	titles = ["Labor"; "Land"; "Rents"; "Wages"]
	# nms = [[L"L_r" L"L_u"], [L"S_r" L"\phi" L"S_{rh}"] , [L"q" L"r"], [L"w_r" L"w_u"]]
	nms = [[L"L_r" L"L_u"], [L"S_r" L"\phi" L"S_{rh}"] , [L"r"], [L"w_r" L"w_u"]]

	plt = Any[]
	if isnothing(year)
		for i in 1:length(sims)
			x = @linq s |>
				where((:variable .∈ Ref(sims[i])) .& (:region .== ik))
			px = @df x plot(:year, :value, group = :variable,
							title = titles[i],
							titlefontsize=10,
							label=nms[i],
							legend = :topleft,
							linewidth=2,marker = (:circle,3))
			push!(plt, px)
		end
		ti = plot(title = "Time Series Region $ik", grid = false, showaxis = false, bottom_margin = -30Plots.px)
		pl = plot(ti,plot(plt...,layout = (2,2), link = :x), layout = @layout([A{0.05h}; B]))
		# return extrema of each subplot
		xlims = [plt[i].subplots[1][:xaxis][:extrema] for i in 1:length(plt)]
		ylims = [plt[i].subplots[1][:yaxis][:extrema] for i in 1:length(plt)]
		xlims = [[xlims[i].emin,xlims[i].emax] for i in 1:length(plt)]
		ylims = [[ylims[i].emin,ylims[i].emax] for i in 1:length(plt)]
		xlims = [(xlims[i][1] .- 0.05*diff(xlims[i])[1],xlims[i][2] .+ 0.05*diff(xlims[i])[1]) for i in 1:length(plt)]
		ylims = [(ylims[i][1] .- 0.05*diff(ylims[i])[1],ylims[i][2] .+ 0.05*diff(ylims[i])[1]) for i in 1:length(plt)]
		return pl,xlims,ylims
	else
		@assert isa(xlim, AbstractArray) && isa(ylim, AbstractArray)
		for i in 1:length(sims)
			x = @linq s |>
				where((:variable .∈ Ref(sims[i])) .& (:region .== ik) .& (:year .<= year))
			px = @df x plot(:year, :value, group = :variable,
							title = titles[i],
							titlefontsize=10,
							label=nms[i],
							xlims = xlim[i],
							ylims = ylim[i],
							legend = :topleft,
							linewidth=2,marker = (:circle,3))
			push!(plt, px)
		end
		ti = plot(title = "Time Series Region $ik in $year", grid = false, showaxis = false, bottom_margin = -30Plots.px)
		pl = plot(ti,plot(plt...,layout = (2,2), link = :x), layout = @layout([A{0.05h}; B]))

		return pl
	end

end


function plot_ts_xsect(C::Vector{Country},pp::Vector{Param},ik::Int)

	# prepare a TS
	df = dataframe(C)
	s = stack(df, Not([:year, :region]))
	ts0,xlims,ylims = single_TS(s,ik)  # to get extrema on each subplot right



	# loop over time
	anim = Animation()
	for (jt,it) in enumerate(C[1].T)

		# make a TS up to period jt: keeping extrema fixed, however, so the plot "grows" nicely
		ts = single_TS(s,ik, year = it, xlim = xlims, ylim = ylims )

		# make a x-section
		xs = plot_space(C[jt].R[ik],pp[ik],jt, fixy = true)

		# combine both
		o = plot(ts, xs, layout = (1,2), size = (1300,500))

		# push onto an anim
		frame(anim)

	end
	g = gif(anim, joinpath(dbpath,"TS-xsection-$ik.gif"),fps=0.5)
	# gif(anim, joinpath(dbpath,"TS-xsection-$ik.gif"),fps=0.5, variable_palette = true)
	return anim
end


function plot_ts_impl(mm::Vector{Country},pp::Vector{Param},ik::Int,it::Int)

	df = dataframe(mm)
	s = stack(df, Not([:year, :region]))
	ts0,xlims,ylims = LandUse.single_TS(s,ik)  # to get extrema on each subplot right

	single_TS(s,ik, year = it, xlim = xlims, ylim = ylims )

end

function plot_dyn(mm::Vector{Model},p::Param,it::Int)
	plot(plot_dyn_impl(mm,p,it)..., layout = (2,3), link = :x)
end

function plot_dynxs(mm::Vector{Model},p::Param,it::Int)
	dyns = plot_dyn_impl(mm,p,it)
	xs = plot_xs_impl(mm[it],p)
	plot([dyns;xs]..., size = (800,400), layout = @layout [grid(2,3){0.8w} grid(2,1){0.2w}])
end


# plotter for optimization tracking

function traceplot(it)
	ft = hcat(Ftrace...)'
	xt = hcat(Xtrace...)'
	K = Int((length(Ftrace[1]) - 3)/2)  # 2K+3 equations
	# nrows = min(size(xt)[1],1000)
	# p1 = plot(ft[1:nrows,:],title = "Ftrace $it",
	p1 = plot(ft,title = "Ftrace $it",
	         label = hcat(["Labor"],reshape(["Land_$i" for i in 1:K],1,K),["rents"],["Urban Good"],reshape(["citysize_$i" for i in 1:K],1,K)),
			 xlabel = "iteration")
	p2 = plot(xt,title = "xtrace",label = hcat(["LS" "r" "pr"],reshape(["SR_$i" for i in 1:K],1,K),reshape(["Lu_$i" for i in 1:K],1,K)),xlabel = "iteration")
	# p2 = plot(xt[1:nrows,:],title = "xtrace",label = hcat(["LS" "r" "pr"],reshape(["SR_$i" for i in 1:K],1,K)),xlabel = "iteration")
	pl = plot(p1,p2,layout = (1,2))
	savefig(pl,joinpath(@__DIR__,"..","images","country_trace$it.pdf"))
end


function plotsol(x)
	K = Int((length(x[1])-3)/2)
	plot(hcat(x...)',lab = hcat(["LS" "r" "pr"],reshape(["SR_$i" for i in 1:K],1,K),reshape(["Lu_$i" for i in 1:K],1,K)),title="solution paths")
end


# nicolas question about rho vs y
function plot_ts0(M::Vector{Region},p::Param)
	df = dataframe(M,p.T)
	# vars = (:ρr, :qr, :Lr, :Lu, :wu0, :wr, :Sr, :Srh, :r, :pr, :ϕ, :icu_input, :iDensity, :icu, :icr, :iτ, :iq, :iy)
	# @df df plot(:year, cols(2:size(df,2)))
	s = stack(df, Not(:year))
	sims = [(:Lr,:Lu,:Sr,:pr,:qr,:r); (:wr, :wu0,:U, :icu , :iy); (:ϕ, :iτ , :ρr ,:Srh)]
	nms = [(L"L_r",L"L_u",L"S_r",L"p_r",L"q_r",L"r"); (L"w_r", L"w_u",:U, L"\int c_u(l) D(l) dl" , L"\int w(l) D(l) dl"); (L"\phi",L"i\tau",L"\rho_r",L"S_{rh}")]
	plts = Any[]
	for i in 1:length(sims)
		x = @linq s |>
			where(:variable .∈ Ref(sims[i]))
		px = @df x plot(:year, :value, group = :variable, layout = length(sims[i]),title = reshape(collect(nms[i]),1,length(nms[i])),titlefontsize=10,leg=false,linewidth=2)
		savefig(px,joinpath(dbpath,"ts_$i.pdf"))
		push!(plts, px)
	end
	return plts

end



function plot_ρ_ts(M::Vector{Region},p::Param)
	plts = Any[]
	push!(plts,plot(p.T,[ M[i].ρr for i in 1:length(p.T)], label="rho_r"))
	push!(plts,plot(p.T,[ pcy(M[i],p) for i in 1:length(p.T)], label="y"))
	push!(plts,plot(p.T,[M[i].ρr / pcy(M[i],p) for i in 1:length(p.T)], label="rho_r/y"))
	pl = plot(plts...,size=(900,450))
	savefig(pl,joinpath(dbpath,"rho_r_y.pdf"))
	pl
end

function plot_y_ts(M::Vector{Region},p::Param)
	x = hcat([ LandUse.pcy(M[i],p) for i in 1:length(p.T)],
	[ M[i].r for i in 1:length(p.T)],
	[ M[i].iy / p.L for i in 1:length(p.T)],
	[ M[i].wr * M[i].Lr / p.L for i in 1:length(p.T)])

	pl = plot(p.T,x,label = ["agg y" "r" "urban y" "rural y"],legend=:topleft)
	savefig(pl,joinpath(dbpath,"y-components.pdf"))
	pl
end



# archive

function plot_static2(m::Region,p::Param)
	lvec = collect(range(0.,area(m),length=100))
	xt = round.(collect(range(0.,area(m),length=3)),digits=2)
	plts = Any[]
	push!(plts,plot(lvec,[ϵ(i,m.ϕ,p) for i in lvec], title = "epsilon(l)",leg = false, xticks = xt,titlefontsize=10,guidefontsize=6))
	push!(plts,plot(lvec,[D(i,p,m) for i in lvec], title = "D(l)",leg = false, xticks = xt,titlefontsize=10,guidefontsize=6))
	push!(plts,plot(lvec,[H(i,p,m) for i in lvec], title = "H(l)",leg = false, xticks = xt,titlefontsize=10,guidefontsize=6))
	push!(plts,plot(lvec,[ρ(i,p,m) for i in lvec], title = "rho(l)",leg = false, xticks = xt,titlefontsize=10,guidefontsize=6))
	push!(plts,plot(lvec,[w(m.Lu,i,m.ϕ,p) for i in lvec], title = "w(l)",leg = false, xticks = xt,titlefontsize=10,tickfont=font(6)))
	push!(plts,plot(lvec,[h(i,p,m) for i in lvec], title = "h(l)",leg = false, xticks = xt,titlefontsize=10,guidefontsize=6))
	plot(plts...)
end

function plot_static(m::Model,p::Param)
	lvec = collect(range(0.,1.0,length=100))
	gr()
	xt = [0.0; 0.5;1.0]
	plts = Any[]
	push!(plts,plot(lvec,[τ(i,m.ϕ,p) for i in lvec], title = "tau",leg = false, xticks = xt,titlefontsize=10,tickfont=font(6)))
	push!(plts,plot(lvec,[w(m.Lu,i,m.ϕ,p) for i in lvec], title = "w(l)",leg = false, xticks = xt,titlefontsize=10,tickfont=font(6)))
	push!(plts,plot(lvec,[cu(i,p,m) for i in lvec], title = "cu(l)",leg = false, xticks = xt,titlefontsize=10,xguidefontsize=6))
	# push!(plts,plot(lvec,[cr(i,p,m) for i in lvec], title = "cr(l)",leg = false, xticks = xt,titlefontsize=10,guidefontsize=6))
	push!(plts,plot(lvec,[cost(i,m.ϕ,p) for i in lvec], title = "cost(l)",leg = false, xticks = xt,titlefontsize=10,guidefontsize=6))
	push!(plts,plot(lvec,[q(i,p,m) for i in lvec], title = "q(l)",leg = false, xticks = xt,titlefontsize=10,guidefontsize=6))
	push!(plts,plot(lvec,[ρ(i,p,m) for i in lvec], title = "rho(l)",leg = false, xticks = xt,titlefontsize=10,guidefontsize=6))
	push!(plts,plot(lvec,[h(i,p,m) for i in lvec], title = "h(l)",leg = false, xticks = xt,titlefontsize=10,guidefontsize=6))
	push!(plts,plot(lvec,[H(i,p,m) for i in lvec], title = "H(l)",leg = false, xticks = xt,titlefontsize=10,guidefontsize=6))
	push!(plts,plot(lvec,[D(i,p,m) for i in lvec], title = "D(l)",leg = false, xticks = xt,titlefontsize=10,guidefontsize=6))
	push!(plts,plot(lvec,[utility(i,p,m) for i in lvec], title = "U(l)",leg = false, xticks = xt,titlefontsize=10,guidefontsize=6))
	push!(plts,plot(lvec,[ϵ(i,m.ϕ,p) for i in lvec], title = "epsilon(l)",leg = false, xticks = xt,titlefontsize=10,guidefontsize=6))
	push!(plts,plot(lvec,[χ(i,m.ϕ,p) for i in lvec], title = "chi(l)",leg = false, xticks = xt,titlefontsize=10,guidefontsize=6))
	plot(plts...)
end

function tplot()
	x,C,cpar,par = LandUse.runk()
	a  = LandUse.anim_cross(C,par)
	a
end
