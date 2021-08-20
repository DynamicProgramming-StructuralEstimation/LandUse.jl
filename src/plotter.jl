
function setup_color(l::Int)
    # setup the colors
    grad   = ColorSchemes.inferno
    colors = [grad[k] for k in range(0,stop=1,length=l)]
    # colors = repeat(colors, 1, length(colors))
    return colors
end

"plot densities for all years"
function plot_dens_years(M::Vector{Region},p::Param)
	x = plot(M[1].ϕmids, M[1].iDensities,leg = false)
	for i in 2:length(M)
		plot!(x, M[i].ϕmids, M[i].iDensities)
	end
	x
end

"""
	dashboard(M::Vector{Region},p::Param,i::Int; objvalue = nothing)

Model dashboard for a single [`Region`](@ref).
"""
function dashboard(M::Vector{Region},p::Param,i::Int; objvalue = nothing)
	pl = LandUse.ts_plots(M,p,fixy = false)
	pc = LandUse.cs_plots(M[i], p, i, objvalue = objvalue)
	po = plot(pl[:Lr_data],pl[:spending],pl[:pr_data],pl[:productivity],
			pl[:n_densities], pl[:densities], pl[:mode], pl[:ctime],
			pl[:phi] , plot(pl[:phi_Lu_data], title = "",ylab = ""), pl[:r_y], pl[:r_rho],
			pc[:ϵ] , pc[:D], pc[:mode] , pl[:speedshares],
			layout = (4,4),size = (1200,800))
	po
end


"""
	dashboard(C::Vector{Country},it::Int)

Dashboard for a [`Country`](@ref) in period `it`
"""
function dashboard(C::Vector{Country},it::Int)
	K = C[1].K
	pl = Any[]

	for ik in 1:K
		cs = cs_plots(C[it].R[ik],C[it].pp[ik], it)
		ts = ts_plots([C[ij].R[ik] for ij in eachindex(C[1].T)],C[it].pp[ik])
		push!(pl, ts[:phi], ts[:pop], ts[:densities], ts[:productivity])
	end
	plot(pl..., layout = (K,4),size = (1100,700))

end

"""
	dashk20(;save = false)

Produce and optionally save model output for 20 city version.
"""
function dashk20(;save = false)

	# compare model to data
	# =====================
	# get a single city to compare
	x,M,p = runm()

	di = Dict()

	ddd = k20_dfexport()
	K = maximum(ddd[1].region)

	ddd[1][!,:scenario] .= "baseline"
	ddd[2][!,:scenario] .= "d1-d2"
	dd = vcat(ddd[1],ddd[2])
	d = ddd[1]
	d2 = ddd[2]


	# baseline vs policy
	# ==================

	# exponential coefficient
	for yr in p.T[[1,5,10,15,19]]
		x = select(subset(dd, :year => x -> x.== yr), :exp_coef => (x -> abs.(x)) => :exp_coef, :region, :scenario)
		di[Symbol("d1d2_coefs_$yr")] = @df x plot(:region, :exp_coef, group = :scenario, marker = (:auto, 3), leg = :bottomright,
									ylab = "exponential exponent (abs)", xlab = "city rank", title = "Exponential coefs in $yr")

		if save
			savefig(di[Symbol("d1d2_coefs_$yr")], joinpath(LandUse.dbplots,"k20-exp-coefs-$yr.pdf"))
		end
	end

	densyears = [1870, 1950,1970,1990,2010]

	dy = subset(d, :year => x -> x .∈ Ref(densyears))
	d2y = subset(d2, :year => x -> x .∈ Ref(densyears))

	meandata = mean(dy.density_data)
	meanmodel = mean(dy.citydensity)
	transform!(dy, :citydensity => (x -> x .* (meandata / meanmodel)) => :ncitydensity)

	meandata = mean(d2y.density_data)
	meanmodel = mean(d2y.citydensity)
	transform!(d2y, :citydensity => (x -> x .* (meandata / meanmodel)) => :ncitydensity)

	dly = transform(dy, :density_data => (x -> log.(x)) => :ldensity_data, 
	                    :citydensity => (x -> log.(x)) => :lcitydensity,
						:ncitydensity => (x -> log.(x)) => :nlcitydensity)

	d2ly = transform(d2y, :density_data => (x -> log.(x)) => :ldensity_data, 
	                    :citydensity => (x -> log.(x)) => :lcitydensity,
						:ncitydensity => (x -> log.(x)) => :nlcitydensity)

	r0 = lm(@formula( ldensity_data ~ nlcitydensity), dly)
	r1 = lm(@formula( ldensity_data ~ nlcitydensity), d2ly)


	di[:dens_mod_data_base] = @df dly plot(:nlcitydensity, :ldensity_data, group = :LIBGEO, leg = :outerright, 
		xlab = "log normalized model density", ylab = "log data density", title = "Urban Density X-section over time",
		color = reshape(palette(:tab20)[1:20], 1,20),
		arrow = 2, 
		markershape = :circle,
		markersize = 4)
	plot!(di[:dens_mod_data_base], x -> coef(r0)[1] + coef(r0)[2] * x, lab = "", linewidth = 2, color = "black")
	annotate!(di[:dens_mod_data_base], [(11.5, 11.8, Plots.text("slope=$(round(coef(r0)[2],digits = 3))", 12))])

	di[:dens_mod_data_d1d2] = @df d2ly plot(:nlcitydensity, :ldensity_data, group = :LIBGEO, leg = :outerright, 
		xlab = "log normalized model density", ylab = "log data density", title = "X-section with d1-d2",
		color = reshape(palette(:tab20)[1:20], 1,20),
		arrow = 2, 
		markershape = :circle,
		markersize = 4)
	plot!(di[:dens_mod_data_d1d2], x -> coef(r1)[1] + coef(r1)[2] * x, lab = "", linewidth = 2, color = "black")
	annotate!(di[:dens_mod_data_d1d2], [(11.5, 11.8, Plots.text("slope=$(round(coef(r1)[2],digits = 3))", 12))])

	
	if save
		savefig(di[:dens_mod_data_base], joinpath(LandUse.dbplots,"k20-xsect-time.pdf"))
		savefig(di[:dens_mod_data_d1d2], joinpath(LandUse.dbplots,"k20-xsect-time-d1d2.pdf"))
	end



	# add rel pop and rel dens
	d = transform(groupby(d,:year), [:Lu, :region]          => ((a,b) -> a ./ a[b .== 1]) => :rel_Lu,
	                                [:citydensity, :region] => ((a,b) -> a ./ a[b .== 1]) => :rel_density)

	# no paris
	nop = subset(d, :region => x -> x .> 1)

	# relative population
	c20 = palette(:tab20)
	pl = @df subset(nop, :region => x -> x.==2) plot(:year, :rel_Lu, leg = :outertopright, linecolor = c20[2], lab = "Lyon",
	ylab = "Population rel to Paris" )
	@df subset(nop, :region => x -> x.==2) scatter!(:datayears, :relpop_data, markercolor = c20[2], lab = "" )
	for ik in 3:K
		dk = subset(nop, :region => x -> x.==ik)
		plot!(pl,dk.year, dk.rel_Lu,  linecolor = c20[ik], lab = dk[1,:LIBGEO] )
		scatter!(pl,dk.datayears, dk.relpop_data,  markercolor = c20[ik], lab = "" )
	end
	title!(pl, "Model (Lines) vs Data (Points)")

	di[:relpop_data] = pl
	

	di[:relpop_nop] = @df nop plot(:year, :rel_Lu, group = :LIBGEO, leg = :outertopright, col = palette(:tab20))
	@df nop scatter!(di[:relpop_nop],:datayears, :relpop_data, group = :LIBGEO, lab = false, col = palette(:tab20))

	di[:rpop] = @df d plot(:year,:Lr, group = :region, ylab = "Share of Rural Population")
	di[:pop] = @df d plot(:year,:Lu, group = :region, ylab = "Urban Population")
	di[:popn] = @df d plot(:year,:Lu_n, group = :region, ylab = "Share of Urban Population")
	di[:relpop] = @df nop plot(:year,:rel_Lu, group = :region,
	                   ylab = "Population relative to Biggest City")
	di[:dens_pop] = @df d plot(:Lu, :citydensity, group = :region, 
	                    xaxis = ("log urban pop", :log10), yaxis = ("log density", :log10), leg = false)

	dparis = subset(d, :region => x-> x.==1)
	di[:dens_paris] = @df dparis plot(:year, :citydensity, label = "Average", lw = 2, title = "Paris Density", ylab = "Density")
	plot!(di[:dens_paris], dparis.year, dparis.d0, label = "Central", lw = 2)

	di[:densn_paris] = @df dparis plot(:year, :avgd_n, label = "Average", lw = 2, title = "Paris Normalized Density", ylab = "Density")
	plot!(di[:densn_paris], dparis.year, dparis.d0_n, label = "Central", lw = 2)
	plot!(di[:densn_paris], dparis.year, dparis.dr_n, label = "Fringe", lw = 2)

	dparis = subset(d, :region => x-> x.==2)
	di[:dens_lyon] = @df dparis plot(:year, :citydensity, label = "Average", lw = 2, title = "Lyon Density", ylab = "Density")
	plot!(di[:dens_lyon], dparis.year, dparis.d0, label = "Central", lw = 2)

	di[:densn_lyon] = @df dparis plot(:year, :avgd_n, label = "Average", lw = 2, title = "Lyon Normalized Density", ylab = "Density")
	plot!(di[:densn_lyon], dparis.year, dparis.d0_n, label = "Central", lw = 2)
	plot!(di[:densn_lyon], dparis.year, dparis.dr_n, label = "Fringe", lw = 2)

	g3 = combine(groupby(d, :year), :citydensity => mean => :avgdensity)

	
	d1 = dataframe(M,p)
	di[:avg_density] = @df g3 plot(:year, :avgdensity, label = "Avg over $K" , title = "Average Densities", lw = 2)
	plot!(di[:avg_density], d1.year, d1.citydensity, label = "Single city", lw = 2)

	di[:density_kvs1] = @df d plot(:year, :citydensity, group = :LIBGEO, label = "", color = :lightgray, leg = false)
	plot!(di[:density_kvs1],d1.year, d1.citydensity, color = :red ,lw = 3, title = "Urban Density")

	di[:relarea_kvs1] = @df d plot(:year, :rel_cityarea, group = :LIBGEO, label = "", color = :lightgray, leg = false)
	plot!(di[:relarea_kvs1],d1.year, d1.rel_cityarea, color = :red ,lw = 3, title = "Urban Area")
	
	di[:Lu_kvs1] = @df d plot(:year, :Lu, group = :LIBGEO, label = "", color = :lightgray)
	plot!(di[:Lu_kvs1],d1.year, d1.Lu, label = "Single city", color = :red ,lw = 3, title = "Urban Population", leg = :topleft)

	# cross section of densities
	sort!(d, [:region, :year])
	di[:density_3d] = wireframe(unique(d.region), unique(d.year), 
	                            reshape(d.citydensity,19,20), 
								camera = (70,40),xlab = "City", 
								ylab = "year",zlab = "Average Density", 
								title = "Density Cross Section over Time")


	

	# within density in year 2020
	de2020 = hcat( Matrix(select(subset(d , :year => x -> x .== 2020 ),  :iDensities))... )
	di[:dens_2020] = wireframe(1:K,1:p.int_bins,de2020, 
	                          xlab = "city rank", ylab = "distance bin",
							  zlab = "density",camera = (70,40),
							  title = "year 2020 within city gradients")

	# 1 by 3 plot of 1 vs k
	di[:p1vsk] = plot(di[:Lu_kvs1],di[:relarea_kvs1],di[:density_kvs1], layout = (1,3), size = (900,300))

	if save
		savefig(di[:p1vsk], joinpath(LandUse.dbplots,"k20-density-vs1.pdf"))
		savefig(di[:dens_2020], joinpath(LandUse.dbplots,"k20-density3D-2020.pdf"))
		savefig(di[:relpop_data], joinpath(LandUse.dbplots,"k20-relpop-model-data.pdf"))
		savefig(di[:density_3d], joinpath(LandUse.dbplots,"k20-density3D.pdf"))
	end


	# compare model to data
	# =====================

	di
end

"single spatial cross sections region"
function cs_plots(m::Region,p::Param,it::Int; fixy = false, objvalue = nothing)

	lvec = collect(range(0.,m.ϕ,length=100))  # space for region ik
	lvec0 = collect(range(0.001,m.ϕ,length=100))  # space for region ik
	setperiod!(p,it)  # set time on parameter!
	ti = p.T[it]  # acutal year for printing
	d = Dict()

	# get data
	ϵd = [ϵ(i,m.ϕ,p) for i in lvec]
	Dd = [D(i,p,m) for i in lvec]
	Hd = [H(i,p,m) for i in lvec]
	ρd = [ρ(i,p,m) for i in lvec]
	qd = [q(i,p,m) for i in lvec]	
	md = [mode(i,p,m.Lu) for i in lvec0]

	# get ratio first over last point
	ϵg = round(ϵd[1]/ϵd[end],digits =1)
	Dg = round(Dd[1]/Dd[end],digits =1)
	Hg = round(Hd[1]/Hd[end],digits =1)
	ρg = round(ρd[1]/ρd[end],digits =1)
	qg = round(qd[1]/qd[end],digits =1)

	# get 90/10 ratio of density
	d1090 = round(m.iDensity_q10 / m.iDensity_q90, digits = 1)

	# run exponential decay model
	ndensities = m.iDensities ./ m.iDensities[1]
	gradient,emod = expmodel(1:p.int_bins, ndensities)
	MSE = round(1000 * mse(emod),digits = 3)

	# get elasticity of speed wrt distance
	# it's a straight line
	spela = diff(log.(md[[1,end]])) ./ diff(log.(lvec0[[1,end]]))

	if isnothing(objvalue)
		manno = (m.ϕ*0.2 , 0.5*maximum(md) , "Elasticity: $(round(spela[1],digits=2))")
	else
		manno = ([m.ϕ*0.2,m.ϕ*0.2 ] , [0.5*maximum(md),0.3*maximum(md)] , ["Elasticity: $(round(spela[1],digits=2))", "Objective: $(round(objvalue,digits = 2))"])
	end


	if fixy
		d[:ϵ] = plot(lvec , ϵd , title = latexstring("\\epsilon(l,$(ti))") , ylims = (2    , 4.1) , linewidth = 2 , leg = false , xlab = "distance" , annotations = (m.ϕ*0.8 , 0.9*maximum(ϵd) , "$(ϵg)x"))                   

		d[:D] = Plots.scatter(1:p.int_bins, ndensities, m = (:circle, :red, 4), leg = false,title = latexstring("D(l,$(ti))")  )
		plot!(d[:D],1:p.int_bins, x -> gradient[1] .* exp.(gradient[2] * x), linewidth = 2, xlab = "distance", 
		            annotations = ([p.int_bins*0.7 ] , [0.9], ["exp.coef=$(round(gradient[2],digits=2))\n10/90=$(d1090)\nMSE=$MSE"]))

		# d[:D] = plot(lvec , Dd , title = latexstring("D(l,$(ti))")         , ylims = (-3   , 60)  , linewidth = 2 , leg = false , xlab = "distance" , annotations = ([m.ϕ*0.7 ] , [0.9*maximum(Dd)], ["10/90=$(round(m.iDensity_q10,digits=1))/$(round(m.iDensity_q90,digits=1))\n=$(d1090)"]))
		# vline!(d[:D],[m.ϕ10, m.ϕ90], color = :red,leg = false)
		d[:H] = plot(lvec , Hd , title = latexstring("H(l,$(ti))")         , ylims = (-0.1 , 15)  , linewidth = 2 , leg = false , xlab = "distance" , annotations = (m.ϕ*0.8 , 0.9*maximum(Hd) , "$(Hg)x"))
		d[:ρ] = plot(lvec , ρd , title = latexstring("\\rho(l,$(ti))")     , ylims = (-0.1 , 10)  , linewidth = 2 , leg = false , xlab = "distance" , annotations = (m.ϕ*0.8 , 0.9*maximum(ρd) , "$(ρg)x"))
		d[:q] = plot(lvec , qd , title = latexstring("q(l,$(ti))")         , ylims = (0.5  , 5.5) , linewidth = 2 , leg = false , xlab = "distance" , annotations = (m.ϕ*0.8 , 0.9*maximum(qd) , "$(qg)x"))
		d[:mode] = plot(lvec0 , md , title = latexstring("mode(l,$(ti))")       , linewidth = 2 , leg = false , xlab = "distance (log scale)" , xscale = :log10, yscale = :log10, annotations = manno)
	else
		d[:ϵ] = plot(lvec , ϵd , title = latexstring("\\epsilon(l,$(ti))")     , linewidth = 2 , leg = false , xlab = "distance" , annotations = (m.ϕ*0.8 , 0.9*maximum(ϵd) , "$(ϵg)x"))
		d[:D] = Plots.scatter(1:p.int_bins, ndensities, m = (:circle, :red, 4), leg = false,title = latexstring("D(l,$(ti))")  )
		plot!(d[:D],1:p.int_bins, x -> gradient[1] .* exp.(gradient[2] * x), linewidth = 2, xlab = "distance", 
		            annotations =  ([p.int_bins*0.7 ] , [0.9], ["exp.coef=$(round(gradient[2],digits=2))\n10/90=$(d1090)\nMSE=$MSE"]))
		d[:H] = plot(lvec , Hd , title = latexstring("H(l,$(ti))")             , linewidth = 2 , leg = false , xlab = "distance" , annotations = (m.ϕ*0.8 , 0.9*maximum(Hd) , "$(Hg)x"))
		d[:ρ] = plot(lvec , ρd , title = latexstring("\\rho(l,$(ti))")         , linewidth = 2 , leg = false , xlab = "distance" , annotations = (m.ϕ*0.8 , 0.9*maximum(ρd) , "$(ρg)x"))
		d[:q] = plot(lvec , qd , title = latexstring("q(l,$(ti))")             , linewidth = 2 , leg = false , xlab = "distance" , annotations = (m.ϕ*0.8 , 0.9*maximum(qd) , "$(qg)x"))
		d[:mode] = plot(lvec0 , md , title = latexstring("mode(l,$(ti))")       , linewidth = 2 , leg = false , xlab = "distance (log scale)" , xscale = :log10, yscale = :log10, annotations = manno)
	end
	d
end

"""
	ts_plots(M,p::Param;fixy = false)

Time series plots for a single [`Region`](@ref) `M`.
"""
function ts_plots(M,p::Param;fixy = false)
	d = dataframe(M,p)
	dd = Dict()
	df = select(d, :year, :hshare => :h, :ushare => :u, :rshare => :r)
	ds = stack(df, Not(:year))
	
	# marker
	mmark = (:circle, 4)

	# colors
	rg = ["red" "green"]
	brg = ["blue" "red" "green"]

	# subsets
	d1880 = subset(d, :year => x -> x .>= 1880)

	# normalizations
	transform!(d1880, :cityarea => firstnorm => :cityarea_n,
	                  :rel_cityarea => firstnorm => :rel_cityarea_n,
	                  :Lu => firstnorm => :Lu_n)

	i1900 = argmin( abs.(p.moments.year .- 1900) )
	i2015 = argmin( abs.(p.moments.year .- 2015) )
	i2010 = argmin( abs.(p.moments.year .- 2015) )
	i1870 = argmin( abs.(p.moments.year .- 1870) )

	h1900 = df[i1900,:h]
	hend = df[i2010,:h]
	dd[:spending] = @df ds plot(:year,:value, group = :variable,
			   linewidth = 2, title = "Spending Shares",
			   ylims = (0.0,0.83), marker = mmark, legend = :right, color = brg,
			   annotations = ([p.T[i1900],p.T.stop],[h1900-0.1, hend-0.1],["$(round(h1900,digits = 3))","$(round(hend,digits = 3))"]))

    # dd[:spending_data] = plot!(dd[:spending], p.moments.year, p.moments[!,[:SpendingShare_Housing, :SpendingShare_Urban,:SpendingShare_Rural]], color = brg)

	ds2 = stack(select(d,:year,:Lu, :Lr), Not(:year))
	dd[:pop] = @df ds2 plot(:year, :value, group = :variable,
					 linewidth = 2, title = "Population", marker = mmark,
					 legend = :right, color = rg)

	dd[:Lr_data] = @df d plot(:year, :Lr ./ p.Lt, label = "model",marker = mmark, color = "blue",linewidth = 2, title = "Rural Employment")
	plot!(dd[:Lr_data], p.moments.year, p.moments.Employment_rural, label = "data", color = "red",linewidth = 2)

	dd[:pr_data] = @df d plot(:year, :pr, label = "model",marker = mmark, color = "blue",linewidth = 2, title = "Rural Price")
	plot!(dd[:pr_data], p.moments.year, p.moments.P_rural, label = "data", color = "red",linewidth = 2)


	dd[:ρ0_real] = @df d plot(:year, :ρ0_real, linewidth = 2, color = "black", leg =false, title = "Central Land Values", marker = mmark, ylims = fixy ? (0.0,29.0) : false )
	dd[:ρ0] = @df d plot(:year, :ρ0, linewidth = 2, color = "black", leg =false, title = "Central Land Values", marker = mmark, ylims = fixy ? (0.0,29.0) : false )
	dd[:qbar] = @df d plot(:year, :qbar, linewidth = 2, color = "black", leg =false, title = "Average House Prices", marker = mmark, ylims = fixy ? (0.0,29.0) : false )
	dd[:qbar_real] = @df d plot(:year, :qbar_real, linewidth = 2, color = "black", leg =false, title = "Average House Prices", marker = mmark, ylims = fixy ? (0.0,29.0) : false )
	dd[:ρ0_y] = @df d plot(:year, :ρ0_y, linewidth = 2, color = "black", leg =false, title = "Central Land Rents over Income", marker = mmark, ylims = fixy ? (0.0,200.0) : false )
	dd[:ρr_y] = @df d plot(:year, :ρr_y, linewidth = 2, color = "black", leg =false, title = "Fringe Land Rents over Income", marker = mmark, ylims = fixy ? (0.0,200.0) : false )
 	dd[:ρr] = @df d plot(:year, :ρr, linewidth = 2, color = "black", leg =false, title = "Fringe Land Values", marker = mmark, ylims = fixy ? (0.19,0.36) : false )
 	dd[:q0] = @df d plot(:year, :q0, linewidth = 2, color = "black", leg =false, title = "Central House Prices", marker = mmark, ylims = fixy ? (0.0,5.0) : false )
 	dd[:qr] = @df d plot(:year, :qr, linewidth = 2, color = "black", leg =false, title = "Fringe House Prices", marker = mmark, ylims = fixy ? (0.5,1.5) : false )
 	dd[:qr_real] = @df d plot(:year, :qr_real, linewidth = 2, color = "black", leg =false, title = "Fringe House Prices", marker = mmark, ylims = fixy ? (0.5,1.5) : false )
	dd[:hr] = @df d plot(:year, :hr, linewidth = 2, color = "black", leg =false, title = "Fringe Housing Demand", marker = mmark, ylims = fixy ? (0.0,5.0) : false )

	df = @linq d |>
		 select(:year,:rr,:r, :y ,:iq, :rr_real, :ru_real, :pop) |>
		 transform(r_y = 100*(:r .* :pop ./ :y), rr_y = 100*(:rr ./ :y), ru_y = 100*(:iq ./ :y))
	dd[:r_y] = @df df plot(:year, [:r_y, :ru_y,:rr_y], labels = ["Total" "Urban" "Rural"],
	            linewidth = 2, title = "Rents as % of Income", marker = mmark,
				ylims = fixy ? (0,50) : false, color = brg,legend = :bottomright)

	df = @linq d |>
		 select(:year,:rr_real, :ru_real) |>
		 transform(rr_real_1 = :rr_real ./ :rr_real[1], ru_real_1 = :ru_real ./ :ru_real[1])
	dd[:rents_real] = @df df plot(:year, [:ru_real, :rr_real ], labels = ["Rural Rents" "Urban Rents"],
	            linewidth = 2, title = "Real Rents", marker = mmark,
				ylims = fixy ? (0,50) : false, legend = :left, color = rg)
	dd[:rents_real_1] = @df df plot(:year, [:ru_real_1,:rr_real_1], labels = ["Rural Rents" "Urban Rents"],
				linewidth = 2, title = "Real Rents", marker = mmark, ylims = fixy ? (0,50) : false,
				legend = :left, color = rg)


	y1950 = findmin(abs.(d.year .- 1950))[2]
	dtemp = transform(select(d, :year, :r, :r_real, :ρ0, :ρr, :q0, :qr, :hr, :Hr, :qbar, :qq1),
	                  :r => ((x)  -> 100*(x ./ x[y1950])) => :r,
	                  :r_real => ((x)  -> 100*(x ./ x[y1950])) => :r_real,
	                  :ρ0 => ((x) -> 100*(x ./ x[y1950])) => :ρ0,
	                  :ρr => ((x) -> 100*(x ./ x[y1950])) => :ρr,
	                  :q0 => ((x) -> 100*(x ./ x[y1950])) => :q0,
	                  :qq1 => ((x) -> 100*(x ./ x[y1950])) => :qq1,
	                  :hr => ((x) -> 100*(x ./ x[y1950])) => :hr,
	                  :Hr => ((x) -> 100*(x ./ x[y1950])) => :Hr,
	                  :qbar => ((x) -> 100*(x ./ x[y1950])) => :qbar,
	                  :qr => ((x) -> 100*(x ./ x[y1950])) => :qr
					  )
	dd[:r] = @df dtemp plot(:year, :r, linewidth = 2, color = "black", leg =false, title = "Land Rents",ylab = "1950=100", marker = mmark, ylims = fixy ? (95,500) : false)
	dd[:r_real] = @df dtemp plot(:year, :r_real, linewidth = 2, color = "black", leg =false, title = "Real Land Rents",ylab = "1950=100", marker = mmark, ylims = fixy ? (95,500) : false)
	dd[:r_rho] = @df dtemp plot(:year, [:r, :ρ0],linewidth = 2, lab = [L"r" L"\rho_0" ],
	                           ylab = "1950=100",
							   title = "Rents and Land Value",
							   leg = :topleft,
							   marker = mmark, ylims = fixy ? (50,500) : false)

    dd[:q0_qr] = @df dtemp plot(:year, [:q0, :qr],linewidth = 2, lab = [L"q_0" L"q_r"],
	                           title = "Central and Fringe House Prices (1945=100)",
							   leg = :topleft,
							   marker = mmark, ylims = fixy ? (85,150) : false)
   dd[:q0_qr_qq1] = @df dtemp plot(:year, [:q0, :qr , :qq1],linewidth = 2, lab = [L"q_0" L"q_r" L"q_q1"],
	                           title = "Central and Fringe House Prices (1945=100)",
							   leg = :topleft,
							   marker = mmark, ylims = fixy ? (85,150) : false)
   dd[:qbar_100] = @df dtemp plot(:year, :qbar, linewidth = 2, color = "black", leg =false, title = "Average House Prices", marker = mmark, ylims = fixy ? (0.0,29.0) : false )

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
					  linewidth = 2, title = "Log Productivity", ylims = fixy ? (0,20) : false, marker = mmark,
					  legend = :left, yscale = :log10)
	# ds4 = stack(select(d,:year,:ϕ), Not(:year))
	ds4 = select(d,:year,:cityarea)
	incphi = round(d1880.rel_cityarea[end] / d1880.rel_cityarea[1],digits = 1)

	dd[:phi] = @df d1880 plot(:year, :rel_cityarea,
					 linewidth = 2, title = "Rel City Area. 2010=$(round(d.rel_cityarea[i2010],digits=2))", color = "black",
					 leg = fixy ? :topleft : false, marker = mmark, annotate = (p.T[end],0.2*maximum(d.rel_cityarea),"$(round(incphi,digits=1))x"),
					 ylims = fixy ? (0.0,0.15) : false)


	dd[:phi_Lu] = @df stack(select(d1880, :year, :Lu_n => :population, :rel_cityarea_n => :rel_area), Not(:year)) plot(:year, :value, 
	                       group = :variable, leg = :left, linewidth = 2,
						   linecolor = ["orange" "red" ], marker = (:circle, ["orange" "red" ]),
						   yscale = :identity, formatter = y->string(round(Int,y)),
						   title = "Urban Population and Area (rel to rural)",
						   ylab = "data (1876) = model (1880) = 1")

	pad = copy(p.poparea_data)
	transform!(pad, :population => firstnorm => :population, 
	                :area => firstnorm => :area)
	dd[:phi_Lu_data] = scatter!(dd[:phi_Lu], pad.year, [pad.population, pad.area],
	                            marker = (:star5, ["orange" "red" ]), label = ["pop data" "area data"])					   
    
	dd[:Sr] = @df d plot(:year, :Sr,
				  linewidth = 2, color = "black",title = "Agricultural Land",
				  leg = false, marker = mmark, ylims = fixy ? (0.6,1.0) : false)
	
    dens = select(d,:year,:d0, :dr,:citydensity)
	df4 = dens
	ds4 = stack(dens, Not(:year))
	# ds4 = stack(select(df4,:year, :avgd), Not(:year))
	dd[:n_densities] = @df stack(select(d,:year,:d0_n, :dr_n, :avgd_n), Not(:year)) plot(:year, :value, group = :variable,
 					 linewidth = 2, title = "Normalized Densities", color = brg,
					 leg = :topright, ylims = fixy ? (0,300) : false)


    incdens = d.citydensity[i1870] / d.citydensity[i2010]
    diffdens = d.citydensity[i1870] - d.citydensity[i2010]
	ancdens = d.citydensity[end] + 0.2 * diffdens
	dd[:avdensity] = @df d plot(:year, :citydensity,
					 linewidth = 2, title = "Avg Density", marker = mmark,
					 legend = false,color = "black", ylims = fixy ? (0,200) : false, annotations = (p.T.stop,ancdens,"$(round(incdens,digits = 1))x"))

	dd[:densities] = @df stack(select(d,:year,:d0, :dr, :citydensity), Not(:year)) plot(:year, :value, group = :variable,
					 linewidth = 2, title = "Densities", leg = :topright, ylims = fixy ? (0,300) : false,
					 annotations = (p.T[2],0.2*maximum(d.d0),"$(round(incdens,digits = 2))x"))


    dss = stack(select(d,:year,:mode0,:modeϕ,:imode), Not(:year))
	facs = combine(groupby(dss,:variable),:value => (x -> x[end]/x[1]) => :factor)
	labs = String.(facs.variable) .* " : " .* string.(round.(facs.factor,digits=1))  .* "x"

	dss = [d[!,:year] mapcols(x -> x ./ x[1],select(d,:mode0,:modeϕ,:imode,:ctime0,:ctimeϕ,:ictime))]
	rename!(dss, :x1 => :year)
	dmode = stack(select(dss,:year,:mode0,:modeϕ,:imode) , Not(:year))

	dd[:mode] = @df dmode plot(:year, :value, group = :variable,
					 linewidth = 2, title = "mode increase", leg = :left, ylims = fixy ? (0,1.4) : false , marker = mmark)
    # dss = stack(select(d,:year,:ctime0,:ctimeϕ,:ictime), Not(:year))
	# facs = combine(groupby(dss,:variable),:value => (x -> x[end]/x[1]) => :factor)
	# labs = String.(facs.variable) .* " : " .* string.(round.(facs.factor,digits=1))  .* "x"

	dtime = stack(select(dss,:year,:ctime0,:ctimeϕ,:ictime) , Not(:year))
	dd[:ctime] = @df dtime plot(:year, :value, group = :variable,
					 # linewidth = 2, title = "commute time", leg = :topleft, labels = reshape(labs,1,3) )
					 linewidth = 2, title = "ctime increase", leg = :topleft, marker = mmark)

	# distribution of modes chosen in each year to get a view of how many people choose which mode in which year
	z = select(d, :iDensity, :iDensitySpeeds , :year, [:iDensity, :iDensitySpeeds] => ((a,b) -> b ./ a) => :SpeedShares)
	z = hcat(select(d,:year), DataFrame(hcat(Array(select(z, :SpeedShares))...)',[:walk, :bus, :car]))
	sy = stack(z, Not(:year))
	dd[:speedshares] = @df sy groupedbar(:year, :value, group = :variable)

	# plot(pl,pl2,pl3,pl4,pl5, layout = (1,5),size = (700,250))
	# plot(pl2,pl6,pl4,pl7 ,layout = (1,4),size = (900,250))
	dd
end


function plot_i0k(di::Dict;it = 19)

	p0 = LandUse.Param(par = di, use_estimatedθ = false)

	x,M,p = LandUse.run(p0, estimateθ = false)
	pl = LandUse.ts_plots(M,p0,fixy = false)
	pc = LandUse.cs_plots(M[it], p0, it)

						# multi country	
	x,C,p = runk(par = merge(Dict(:K => 2, :kshare => [0.5,0.5], :factors => [1.0,1.05]),di))
	x = impl_plot_slopes(C)
	pl1 = plot(x[3])

	pl2 = plot(pl[:Lr_data],pl[:spending],pl[:pr_data],pl[:productivity],
						     pl[:n_densities], pl[:densities], pl[:mode], pl[:ctime],
							 pl[:phi] , pl[:qbar_real], pl[:r_y], pl[:r_rho],
							 pc[:ϵ] , pc[:D], pc[:q] , pc[:H],
							 layout = (4,4))
	plot(pl2, pl1, size = (1600,700), layout = @layout [a{0.7w} b{0.3w}])	
end
