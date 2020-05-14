

"""
https://github.com/floswald/LandUse.jl/issues/22
"""
function i22()
	sbs = 0.0:0.05:0.5
	cbs = 0.0:0.05:0.5
	sis = 0.1:0.1:0.99
	ags = 1.0:0.1:1.3
	θus = 1.0:0.05:1.2
	θrs = 1.0:0.05:1.2
	p = Param()
	@manipulate for sb in slider(sbs, label = "sbar", value = p.sbar),
		            cb in slider(cbs, label = "cbar", value = p.cbar),
		            sig in slider(sis, label = "σ", value = p.σ),
		            agg in slider(ags, label = "agg_g", value = p.θagg_g),
		            θut in slider(θus, label = "θu", value = 1.0),
		            θrt in slider(θrs, label = "θr", value = 1.0)
		θr = [θrt for i in p.T]
		θu = [θut for i in p.T]
	    x,M,p = run(Region,
		             Param(par = Dict(:ϵs => 0.0, :ϵsmax => 0.0,
					                  :sbar => sb, :cbar => cb, :σ => sig,
									  :θagg_g => agg, :θut => θu, :θrt => θr)))
	    d = dataframe(M,p.T)
	    df = @linq d |>
	         select(:year,:Ch,:Cu,:Cr,:C ) |>
	         transform(h = :Ch ./ :C,u = :Cu ./ :C,r = :Cr ./ :C) |>
	         select(:year, :h, :u , :r)
	    ds = stack(df, Not(:year))
	    pl = @df ds plot(:year,:value, group = :variable,
	               linewidth = 2, title = "Spending Shares")

		ds2 = stack(select(d,:year,:Lu, :Lr), Not(:year))
		pl2 = @df ds2 plot(:year, :value, group = :variable,
		                 linewidth = 2, title = "Population")
	    d3  = @linq d |>
		        select(:year,:θu, :θr) |>
				transform(theta_u = :θu, theta_r = :θr) |>
				select(:year,:theta_u, :theta_r)
	    ds3 = stack(d3, Not(:year))
	    pl3 = @df ds3 plot(:year, :value, group = :variable,
					      linewidth = 2, title = "Productivity")
		plot(pl,pl2,pl3, layout = (1,3),size = (700,250))
    end
end


function i0()
	epsimax = 0.0:10.0
	gfac = 1.00:0.05:2.0
	sigmas = 0.1:0.1:0.99
	sbs = 0.0:0.1:0.5

	@manipulate for growthtype = Dict("orig" => 1,"u-r const" => 4),
					ugrowth in slider(gfac,value = 1.01, label = "u-growthrate"),
					rgrowth in slider(gfac,value = 1.01, label = "r-growthrate"),
					epsim in slider(epsimax, value = 1.0, label = "ϵr"),
					sigma in slider(sigmas, value = 0.9, label = "sigma"),
					sb in slider(sbs, label = "sbar")

					if growthtype == 1
						p0 = LandUse.Param(par = Dict(:ϵr => epsim,:ϵsmax => 0.0,:σ => sigma , :sbar => sb))
						x,M,p = LandUse.run(LandUse.Region,p0)
						LandUse.plot_ts(M,p0)

					elseif growthtype == 2
						println("does not work")
						# p0 = LandUse.Param(par = Dict(:θug => [ugrowth for i in 1:14],
						# 					 :ϵsmax => epsim,
						# 					 :σ => sigma ))
						# x,M,p = LandUse.run(LandUse.Region,p0)
						# LandUse.plot_ts(M,p0)
					#
					elseif growthtype == 3
						println("does not work")
					# 		p0 = LandUse.Param(par = Dict(:θrg => [rgrowth for i in 1:14],
					# 							 :ϵsmax => epsim,
					# 							 :σ => sigma ))
					# 		x,M,p = LandUse.run(LandUse.Region,p0)
					# 		LandUse.plot_ts(M,p0,time)
					#
					elseif growthtype == 4
						p0 = LandUse.Param(par = Dict(:θut => [ugrowth for i in 1:14],
											 :θrt => [rgrowth for i in 1:14],
						                     :ϵr => epsim,:ϵsmax => 0.0,
											 :σ => sigma, :sbar => sb ))
					    x,M,p = LandUse.run(LandUse.Region,p0)
						LandUse.plot_ts(M,p0)
					end
	end
end


function i1()

	ti = 1:14
	epsimax = 0.0:10.0
	gfac = 1.00:0.05:1.25
	sigmas = 0.1:0.1:0.99

	@manipulate for growthtype = Dict("orig" => 1, "orig-u" => 2, "orig-r" => 3, "u-r const" => 4),
		            time in slider(ti, value = 14, label = "time"),
					ugrowth in slider(gfac,value = 1.24, label = "u-growthrate"),
					rgrowth in slider(gfac,value = 1.24, label = "r-growthrate"),
					epsim in slider(epsimax, value = 0.0, label = "ϵ-slope"),
					sigma in slider(sigmas, value = 0.9, label = "sigma")

					if growthtype == 1
						p0 = LandUse.Param(par = Dict(:ϵsmax => epsim,:σ => sigma ))
						x,M,p = LandUse.run(LandUse.Region,p0)
						LandUse.plot_ts(M,p0,time)

					elseif growthtype == 2
						p0 = LandUse.Param(par = Dict(:θut => [ugrowth for i in 1:14],
											 :ϵsmax => epsim,
											 :σ => sigma ))
						x,M,p = LandUse.run(LandUse.Region,p0)
						LandUse.plot_ts(M,p0,time)
					#
					# elseif growthtype == 3
					# 		p0 = LandUse.Param(par = Dict(:θrg => [rgrowth for i in 1:14],
					# 							 :ϵsmax => epsim,
					# 							 :σ => sigma ))
					# 		x,M,p = LandUse.run(LandUse.Region,p0)
					# 		LandUse.plot_ts(M,p0,time)
					#
					# elseif growthtype == 4
					# 	p0 = LandUse.Param(par = Dict(:θug => [ugrowth for i in 1:14],
					# 						 :θrg => [rgrowth for i in 1:14],
					# 	                     :ϵsmax => epsim,
					# 						 :σ => sigma ))
					#     x,M,p = LandUse.run(LandUse.Region,p0)
					# 	LandUse.plot_ts(M,p0,time)


					# elseif growthtype == 3
					# 	gg = (LandUse.originalθ[end] - LandUse.originalθ[1]) / 13
					# 	g0 = LandUse.originalθ[1] .+ gg .* [i for i in 0:13]
					# 	g = g0[2:end] ./ g0[1:end-1]
					# 	p0 = LandUse.Param(par = Dict(:θug => [growth for i in 1:14],
					# 						 :θrg => [growth for i in 1:14],
					# 	                     :ϵsmax => epsim))
					#     x,M,p = LandUse.run(p0)
					# 	LandUse.plot_ts_xsect(M,p0,time)

					end
	end
end

function interact1()
	K = 2
	es = 1.0:0.5:10.0
	g2 = 1.0:0.01:1.06

    # mp = @manipulate for θr in slider(θrs, label = "θr"), θu in slider(θus, label = "θu")
    mp = @manipulate for e in slider(es, label = "ϵr", value = 4.0 ),
		                 g in slider(g2, label = "growth2", value = 1.01)
		x = LandUse.issue12_1(e; gf = [1.0,g])
		x[5]
    end
end

function interact2()
	K = 2
	es = 1.0:0.5:10.0
	g2 = 1.0:0.01:1.06

    # mp = @manipulate for θr in slider(θrs, label = "θr"), θu in slider(θus, label = "θu")
    mp = @manipulate for e in slider(es, label = "ϵr", value = 4.0 ),
		                 g in slider(g2, label = "growth2", value = 1.01)
		x = LandUse.issue12_1(e; gf = [1.0,g])
		x[5]
    end
end
